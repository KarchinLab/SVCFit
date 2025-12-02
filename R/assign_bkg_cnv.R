library(GenomicRanges)
library(dplyr)
library(tidyr)

assign_background_cnv <- function(sv_sum, cnv, flank = 1000) {
  # GRanges for CNV segments
  cnv_gr <- GRanges(
    seqnames = cnv$chrom,
    ranges   = IRanges(cnv$start, cnv$end),
    cnv_idx  = seq_len(nrow(cnv))
  )
  
  sv_idx <- seq_len(nrow(sv_sum))
  
  ## 1) Build left and right flank GRanges around each SV
  # left: [POS - flank, POS - 1]
  left_gr <- GRanges(
    seqnames = sv_sum$CHROM,
    ranges   = IRanges(
      start = pmax(1, sv_sum$POS - flank),
      end   = pmax(0, sv_sum$POS - 1)
    ),
    sv_idx  = sv_idx
  )
  # drop empty ranges (can happen near chromosome start)
  left_gr <- left_gr[width(left_gr) > 0]
  
  # right: [END + 1, END + flank]
  right_gr <- GRanges(
    seqnames = sv_sum$CHROM,
    ranges   = IRanges(
      start = sv_sum$END + 1,
      end   = sv_sum$END + flank
    ),
    sv_idx  = sv_idx
  )
  right_gr <- right_gr[width(right_gr) > 0]
  
  ## 2) Helper to get best-overlapping CNV per flank
  get_best <- function(gr, side_label) {
    if (length(gr) == 0L) {
      return(tibble())  # nothing to do
    }
    
    hits <- findOverlaps(gr, cnv_gr)
    
    if (length(hits) == 0L) {
      return(tibble())
    }
    
    ov_df <- as.data.frame(hits) %>%
      mutate(
        width = width(
          pintersect(
            gr[queryHits],
            cnv_gr[subjectHits]
          )
        ),
        side = side_label
      )
    
    # for each SV flank pick CNV with max overlap
    best <- ov_df %>%
      group_by(queryHits, side) %>%
      slice_max(width, with_ties = FALSE) %>%
      ungroup()
    
    best %>%
      transmute(
        sv_idx = gr$sv_idx[queryHits],
        side,
        cna   = cnv$tcn.em[subjectHits],
        minor = ifelse(is.na(cnv$lcn.em[subjectHits]), 0, cnv$lcn.em[subjectHits]),
        major = cna - minor,
        cncf  = cnv$cf.em[subjectHits]
      )
  }
  
  ## 3) Get best CNV for left and right flanks
  bg_tbl <- bind_rows(
    get_best(left_gr,  "left"),
    get_best(right_gr, "right")
  )
  
  # if nothing overlapped, just return sv_sum with NAs
  if (nrow(bg_tbl) == 0L) {
    return(
      sv_sum %>%
        mutate(
          cna_left = NA_real_, minor_left = NA_real_, major_left = NA_real_, cncf_left = NA_real_,
          cna_right = NA_real_, minor_right = NA_real_, major_right = NA_real_, cncf_right = NA_real_
        )
    )
  }
  
  ## 4) Spread left/right into columns and merge back to SVs
  bg_wide <- bg_tbl %>%
    pivot_wider(
      id_cols    = sv_idx,
      names_from = side,
      values_from = c(cna, minor, major, cncf),
      names_glue = "{.value}_{side}"
    )
  
  sv_sum %>%
    mutate(sv_idx = row_number()) %>%
    left_join(bg_wide, by = "sv_idx") %>%
    select(-sv_idx)%>%
    mutate(bkg_cna=rowMeans(.[, c("cna_left", "cna_right")], na.rm = TRUE),
           bkg_cnv=case_when(
             bkg_cna > 2 ~ 'DUP',
             bkg_cna ==2 ~ 'norm',
             bkg_cna < 2 ~ 'DEL'
           ))
}

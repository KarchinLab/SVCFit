#' Assign background CNV segments to structural variants
#'
#' For each SV in \code{sv_sum}, finds the best-overlapping CNV segment within
#' a flanking window on both sides of the breakpoint and annotates the SV with
#' total copy number, minor/major allele counts, and clonal cell fraction from
#' those flanking segments.  The background CNV state is summarised as the mean
#' of the left and right total copy numbers and classified as duplication,
#' normal, or deletion.
#'
#' @param sv_sum data.frame. SV summary table with at minimum the columns
#'   \code{CHROM}, \code{POS}, and \code{END}.
#' @param cnv data.frame. CNV segment table in FACETS format; must contain the
#'   columns \code{chrom}, \code{start}, \code{end}, \code{tcn.em}
#'   (total copy number), \code{lcn.em} (minor allele copy number), and
#'   \code{cf.em} (clonal cell fraction).
#' @param flank Integer. Half-width (bp) of the flanking window searched on
#'   each side of the SV breakpoints for overlapping CNV segments.
#'   Default \code{1000}.
#'
#' @return The input \code{sv_sum} data.frame with the following additional
#'   columns:
#' \describe{
#'   \item{\code{cna_left}, \code{minor_left}, \code{major_left}, \code{cncf_left}}{
#'     Total CN, minor allele CN, major allele CN, and clonal fraction for the
#'     best-overlapping CNV segment in the left flank.}
#'   \item{\code{cna_right}, \code{minor_right}, \code{major_right}, \code{cncf_right}}{
#'     Same as above for the right flank.}
#'   \item{\code{bkg_cna}}{Mean total copy number across left and right flanks
#'     (\code{NA} values ignored).}
#'   \item{\code{bkg_cnv}}{Background CNV classification: \code{"DUP"}
#'     (\code{bkg_cna > 2}), \code{"norm"} (\code{bkg_cna == 2}), or
#'     \code{"DEL"} (\code{bkg_cna < 2}).}
#' }
#' If no CNV segments overlap the flanks the function still returns all rows of
#' \code{sv_sum} with the eight CNV columns set to \code{NA}.
#'
#' @export
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

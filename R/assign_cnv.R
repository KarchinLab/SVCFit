#' assign an SV to a CNV based on genomic region overlap
#'
#' @param sv_sum an object of class 'data frame'. This object stores the output from `sum_sv_info`.
#' @param cnv an object of class 'data frame'. This object stores the cnv data loaded from `load_data`.
#'
#' @returns sth
#' @export
#'
assign_cnv <- function(sv_sum, cnv) {
  # GRanges for SV segments
  sv_seg <- GRanges(
    seqnames = sv_sum$CHROM,
    ranges   = IRanges(
      pmin(sv_sum$POS, sv_sum$END),
      pmax(sv_sum$POS, sv_sum$END)
    ),
    sv_idx   = seq_len(nrow(sv_sum))
  )
  # GRanges for CNV segments
  cnv_gr <- GRanges(
    seqnames = cnv$chrom,
    ranges   = IRanges(cnv$start, cnv$end),
    cnv_idx  = seq_len(nrow(cnv))
  )
  # overlaps
  hits <- findOverlaps(sv_seg, cnv_gr)
  ov_df <- as.data.frame(hits) %>%
    mutate(
      width = width(
        pintersect(
          sv_seg[queryHits],
          cnv_gr[subjectHits]
        )
      )
    )
  # pick best overlap per SV
  best <- ov_df %>%
    group_by(queryHits) %>%
    slice_max(width, with_ties = FALSE) %>%
    ungroup()
  # merge back
  sv_cnv <- best %>%
    transmute(
      row_num = queryHits,
      cna     = cnv$tcn.em[subjectHits],
      minor   = ifelse(is.na(cnv$lcn.em[subjectHits]), 0, cnv$lcn.em[subjectHits]),
      major   = cna - minor,
      cncf    = cnv$cf.em[subjectHits]
    ) %>%
    inner_join(
      sv_sum %>% mutate(row_num = row_number()),
      by = 'row_num'
    )
  return(sv_cnv)
}

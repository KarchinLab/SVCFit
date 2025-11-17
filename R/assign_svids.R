#' assign each SNP to an SV
#'
#' @param sv_phase an object of class 'data frame'. This object stores the output from `parse_snp_on_sv`.
#' @param sv_info an object of class 'data frame'. This object stores the output from `parse_sv_info`.
#' @param flank an object of class 'integer'. This object describes the upper limit of location difference
#' for an SNP to be assigned to an SV
#'
#' @returns sth
#' @export
#'
assign_svids <- function(sv_phase, sv_info, flank = 500) {
  # GRanges for SNPs (sv_phase) and SVs (sv_info)
  snp_gr <- GRanges(
    seqnames = sv_phase$CHROM,
    ranges   = IRanges(sv_phase$POS, sv_phase$POS+1),
    snp_id  = seq_len(nrow(sv_phase))
  )
  sv_gr <- GRanges(
    seqnames = sv_info$CHROM,
    ranges   = IRanges(sv_info$nPOS, sv_info$mPOS),
    sv_id    = sv_info$ID
  )
  
  #Create ±500 bp “start” and “end” windows
  starts  <- pmax(start(sv_gr) - flank, 1)
  ends_st <- start(sv_gr) + flank
  sv_start_win <- GRanges(
    seqnames = seqnames(sv_gr),
    ranges   = IRanges(start=starts, end=ends_st),
    sv_id    = sv_gr$sv_id
  )
  
  ends    <- end(sv_gr) + flank
  starts_en <- pmax(end(sv_gr) - flank, 1)
  sv_end_win <- GRanges(
    seqnames = seqnames(sv_gr),
    ranges   = IRanges(start=starts_en, end=ends),
    sv_id    = sv_gr$sv_id
  )
  
  #Find overlaps for each
  hits_start <- findOverlaps(snp_gr, sv_start_win, ignore.strand=TRUE)
  hits_end   <- findOverlaps(snp_gr, sv_end_win,   ignore.strand=TRUE)
  
  #Turn into a table, tagging region=“start” or “end”
  ## some trans has start bigger than end, so i need to filter
  start_id=sv_info %>% filter(POS < END)
  end_id=sv_info %>% filter(POS > END)
  
  assign_start <- data.frame(
    snp_id = mcols(snp_gr)$snp_id[ queryHits(hits_start) ],
    ID  = mcols(sv_start_win)$sv_id[ subjectHits(hits_start) ],
    region = "start",
    stringsAsFactors = FALSE
  )%>%
    filter(ID %in% start_id$ID)
  
  assign_end <- data.frame(
    snp_id = mcols(snp_gr)$snp_id[ queryHits(hits_end) ],
    ID  = mcols(sv_end_win)$sv_id[ subjectHits(hits_end) ],
    region = "end",
    stringsAsFactors = FALSE
  )%>%
    filter(ID %in% end_id$ID)
  
  assign_id <- rbind(assign_start, assign_end)
  return(assign_id)
}

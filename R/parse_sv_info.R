#' parse SV information
#'
#' @param sv an object of class 'data frame'. This object stores the sv data loaded from `load_data`.
#' @param bnd an object of class 'data frame'. This object stores the output from `proc_bnd`.
#' @param del an object of class 'data frame'. This object stores the output from `proc_bnd`.
#' @param QUAL_tresh an object of class 'integer'. This object describes the minimum QUAL score
#' to include for analysis.
#' @param min_alt an object of class 'integer'. This object describes the minimum alternative read
#' count for a structural varinats to have to be included for anlaysis.
#'
#' @returns sth
#' @export
#'
parse_sv_info <- function(sv, bnd, del, QUAL_thresh=100, min_alt=2) {
  sv_info <- sv %>%
    filter(grepl('INS', INFO),
           FILTER %in% c("PASS") | QUAL > QUAL_thresh,
           !ID %in% del$ID)%>%
    rowwise()%>%
    mutate(sv_ref = setNames(strsplit(tumor,  ":")[[1]], strsplit(FORMAT, ":")[[1]])["RO"],
           sv_alt = setNames(strsplit(tumor,  ":")[[1]], strsplit(FORMAT, ":")[[1]])["AO"],
           sv_ref=as.integer(sv_ref),
           sv_alt=as.integer(sv_alt),
           #sv_ref=as.integer(gsub(".*:(\\d+)(:([^:]*)){9}$","\\1", tumor)),
           #sv_alt=as.integer(gsub(".*:(\\d+)(:([^:]*)){8}$","\\1", tumor)),
           classification=gsub(".*SVTYPE=(\\w+).*","\\1", INFO),
           SVID=gsub("((:[^:\\D:]){4}$)","", ID))%>%
    group_by(SVID, CHROM)%>%
    mutate(
      #inv has two records, reduce them to one
      sv_ref=ifelse(classification=='INV', mean(as.integer(sv_ref)), sv_ref),
      sv_alt=ifelse(classification=='INV', mean(as.integer(sv_alt)), sv_alt),
      # half the ref reads for dup and del to unify equations
      sv_ref=ifelse(classification %in% c('DUP','DEL'), sv_ref/2, sv_ref))%>%
    ungroup()%>%
    left_join(bnd)%>%
    filter(sv_alt > min_alt)%>%
    mutate(classification=case_when(
      class=='rtl' ~ 'RBND',
      class=='cuttl' ~ 'CUTBND',
      class=='coptl' ~ 'COPBND',
      is.na(class) ~ classification),
      END=ifelse(grepl("BND", ID),mPOS, gsub(".*END=(\\d+).*", "\\1", INFO)),
      END=as.numeric(END),
      END=ifelse(abs(POS-END)<50, nPOS, END),
      nPOS=ifelse(!grepl('BND',ID), POS, nPOS),
      mPOS=ifelse(!grepl('BND',ID), END, mPOS),
      chr2=ifelse(grepl('BND', ID), receiver, CHROM))%>%
    filter(!is.na(END))%>%
    select(CHROM, chr2,POS, END,nPOS, mPOS, ID, sv_ref, sv_alt, class, donor, receiver, classification) %>%
    mutate(CHROM=ifelse(grepl('BND', classification), donor, CHROM),
           chr2=ifelse(grepl('BND', classification), receiver, chr2),
           POS=ifelse(grepl('BND', classification), nPOS, POS),
           END=ifelse(grepl('BND', classification), mPOS, END))

  return(sv_info)
}

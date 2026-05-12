#' parse SV information
#'
#' @param sv an object of class 'data frame'. This object stores the sv data loaded from `load_data`.
#' @param bnd an object of class 'data frame'. This object stores the output from `proc_bnd`.
#' @param del an object of class 'data frame'. This object stores the output from `proc_bnd`.
#' @param QUAL_thresh an object of class 'integer'. This object describes the minimum QUAL score
#' to include for analysis.
#' @param min_alt an object of class 'integer'. This object describes the minimum alternative read
#' count for a structural varinats to have to be included for anlaysis.
#'
#' @return A data.frame of filtered SVs with columns \code{CHROM}, \code{chr2},
#'   \code{pos2}, \code{POS}, \code{END}, \code{nPOS}, \code{mPOS}, \code{ID},
#'   \code{sv_ref}, \code{sv_alt}, \code{class}, \code{donor}, \code{receiver},
#'   \code{classification}, \code{zygosity_tmp}, and \code{mate}.  BND records
#'   are annotated with translocation class; INV records share read counts across
#'   mate pairs.
#' @export
#'
# parse_sv_info <- function(sv, bnd, del) {
#   sv_info <- sv %>%
#     filter(!grepl('SVTYPE=INS', INFO),
#            FILTER %in% c("PASS") | QUAL >100,
#            !ID %in% del$ID)%>%
#     rowwise()%>%
#     mutate(sv_ref = setNames(strsplit(tumor,  ":")[[1]], strsplit(FORMAT, ":")[[1]])["RO"],
#            sv_alt = setNames(strsplit(tumor,  ":")[[1]], strsplit(FORMAT, ":")[[1]])["AO"],
#            sv_ref=as.integer(sv_ref),
#            sv_alt=as.integer(sv_alt),
#            # sv_ref=as.integer(gsub(".*:(\\d+)(:([^:]*)){9}$","\\1", tumor)),
#            # sv_alt=as.integer(gsub(".*:(\\d+)(:([^:]*)){8}$","\\1", tumor)),
#            classification=gsub(".*SVTYPE=(\\w+).*","\\1", INFO),
#            SVID=gsub("((:[^:\\D:]){4}$)","", ID))%>%
#     group_by(SVID, CHROM)%>%
#     mutate(
#       #inv has two records, reduce them to one
#       sv_ref=ifelse(classification=='INV', mean(as.integer(sv_ref)), sv_ref),
#       sv_alt=ifelse(classification=='INV', mean(as.integer(sv_alt)), sv_alt),
#       # half the ref reads for dup and del to unify equations
#       sv_ref=ifelse(classification %in% c('DUP','DEL'), sv_ref/2, sv_ref))%>%
#     ungroup()%>%
#     left_join(bnd)%>%
#     filter(sv_alt > 2)%>%
#     mutate(classification=case_when(
#       class=='rtl' ~ 'RBND',
#       class=='cuttl' ~ 'CUTBND',
#       class=='coptl' ~ 'COPBND',
#       is.na(class) ~ classification),
#       END=ifelse(grepl("BND", INFO),mPOS, gsub(".*END=(\\d+).*", "\\1", INFO)),
#       END=as.numeric(END),
#       END=ifelse(abs(POS-END)<50, nPOS, END),
#       nPOS=ifelse(!grepl('BND',INFO), POS, nPOS),
#       mPOS=ifelse(!grepl('BND',INFO), END, mPOS),
#       chr2=ifelse(grepl('BND', INFO), receiver, CHROM))%>%
#     filter(!is.na(END))%>%
#     select(CHROM, chr2,POS, END,nPOS, mPOS, ID, sv_ref, sv_alt, class, donor, receiver, classification)%>%
#     mutate(CHROM=ifelse(grepl('BND', classification), donor, CHROM),
#            chr2=ifelse(grepl('BND', classification), receiver, chr2),
#            POS=ifelse(grepl('BND', classification), nPOS, POS),
#            END=ifelse(grepl('BND', classification), mPOS, END))
#   
#   return(sv_info)
# }
parse_sv_info <- function(sv, bnd, del, QUAL_thresh=100, min_alt=2) {
  sv_info <- sv %>%
    filter(!grepl('SVTYPE=INS', INFO),
           FILTER %in% c("PASS") | QUAL > QUAL_thresh,
           !ID %in% del$ID)%>%
    left_join(bnd)%>%
    rowwise()%>%
    #collect read count
    mutate(sv_ref = setNames(strsplit(tumor,  ":")[[1]], strsplit(FORMAT, ":")[[1]])["RO"],
           sv_alt = setNames(strsplit(tumor,  ":")[[1]], strsplit(FORMAT, ":")[[1]])["AO"],
           sv_ref=as.integer(sv_ref),
           sv_alt=as.integer(sv_alt),
           classification=gsub(".*SVTYPE=(\\w+).*","\\1", INFO),
           mate={idx <- which(abs(POS-sv$POS)<50); if(grepl('INV', classification) && length(idx)>0) paste0('inv', idx[1L]) else mate})%>%
    group_by(mate, CHROM)%>%
    mutate(
      #inv has two records, reduce them to one
      sv_ref=ifelse(classification=='INV', mean(as.integer(sv_ref)), sv_ref),
      sv_alt=ifelse(classification=='INV', mean(as.integer(sv_alt)), sv_alt),
      # half the ref reads for dup and del to unify equations
      sv_ref=ifelse(classification %in% c('DUP','DEL'), sv_ref/2, sv_ref))%>%
    ungroup()%>%
    filter(sv_alt > min_alt)%>%
    mutate(classification=case_when(
      class=='rtl' ~ 'RBND',
      class=='cuttl' ~ 'CUTBND',
      class=='coptl' ~ 'COPBND',
      is.na(class) ~ classification),
      END=ifelse(grepl("BND", INFO),mPOS, gsub(".*END=(\\d+);.*", "\\1", INFO)),
      END=as.numeric(END),
      END=ifelse(abs(POS-END)<50, nPOS, END),
      nPOS=ifelse(!grepl('BND',INFO), POS, nPOS),
      mPOS=ifelse(!grepl('BND',INFO), END, mPOS),
      chr2=ifelse(grepl('BND', INFO), receiver, CHROM),
      zygosity_tmp=gsub('^(\\d/\\d):.*','\\1', tumor),
      z1=gsub('(\\d)/\\d','\\1', zygosity_tmp),
      z2=gsub('\\d/(\\d)','\\1', zygosity_tmp),
      zygosity_tmp=ifelse(z1==z2, 'hom','het'))%>%
    filter(!is.na(END))%>%
    select(CHROM, chr2,pos2,POS, END,nPOS, mPOS, ID, sv_ref, sv_alt, class, donor, receiver, classification, zygosity_tmp,mate)%>%
    mutate(CHROM=ifelse(grepl('BND', classification), donor, CHROM),
           chr2=ifelse(grepl('BND', classification), receiver, chr2),
           POS=ifelse(grepl('BND', classification), nPOS, POS),
           END=ifelse(grepl('BND', classification), mPOS, END))%>%                                                                   
    mutate(bnd_group = case_when(
      grepl('BND', classification) ~ sub("(MantaBND:\\d+):.*", "\\1", ID),
      TRUE                         ~ ID
    ))%>%
    group_by(bnd_group)%>%
    mutate(
    best_sv_alt = max(sv_alt, na.rm = TRUE),
    best_sv_ref = sv_ref[which.max(sv_alt)]
  ) %>%
  slice_max(sv_ref, n=1, with_ties=FALSE) %>%
  mutate(
    use_mate = grepl("BND", classification) & sv_alt < best_sv_alt,
    sv_ref   = ifelse(use_mate, best_sv_ref, sv_ref),
    sv_alt   = ifelse(use_mate, best_sv_alt, sv_alt)
  ) %>%
  select(-best_sv_alt, -best_sv_ref, -use_mate) %>%
  ungroup()
    #mutate(sv_ref_pool = sum(as.integer(sv_ref)),
    #         sv_alt_pool = sum(as.integer(sv_alt)))%>%
    #  slice_max(sv_ref, n=1, with_ties=FALSE)%>%
    #  mutate(sv_ref = sv_ref_pool, sv_alt = sv_alt_pool)%>%
    #  select(-sv_ref_pool, -sv_alt_pool)%>%
    #  ungroup()
    ##slice_max(sv_alt, n=1, with_ties=FALSE)%>%
    ##ungroup()
  
  return(sv_info)
}

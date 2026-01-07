#' sumarise SV information (zygosity, phasing, reads)
#'
#' @param sv_phase an object of class 'data frame'. This object stores the output from `parse_snp_on_sv`.
#' @param assign_id an object of class 'data frame'. This object stores the output from `assign_svids`.
#' @param sv_info an object of class 'data frame'. This object stores the output from `parse_sv_info`.
#'
#' @returns sth
#' @export
#'
sum_sv_info <- function(sv_phase, assign_id, sv_info){
  sv_sum    <- sv_phase %>%
    left_join(assign_id, by=c('snp_id'))%>%
    filter(!is.na(ID))%>%
    # use het snp to decide zygosity and phasing
    mutate(sv_zy = ifelse(a_count == 1, "het", "hom"),
           sv_phase = ifelse(allele == REF, "mat", "pat"))%>%
    group_by(ID, CHROM) %>%
    summarise(
      POS       = mean(POS),
      snp_alt   = mean(snp_alt, na.rm=TRUE),
      snp_ref   = mean(snp_ref, na.rm=TRUE),
      dep       = mean(dep, na.rm=TRUE),
      onsv_alt  = mean(onsv_alt, na.rm=TRUE),
      onsv_ref  = mean(onsv_ref, na.rm=TRUE),
      zygosity_snp = ifelse(sum(sv_zy == "hom") > sum(sv_zy == "het"), "hom", "het"),
      sv_phase = ifelse(sum(sv_phase == "mat") > sum(sv_phase == "pat"), "mat", "pat"),
      .groups="drop"
    )%>%
    select(-POS)%>%
    right_join(sv_info, by=c('ID','CHROM'))%>%
    #assign info to SV with no het snps
    mutate(no_snp = ifelse(is.na(snp_alt), 1, 0),
           sv_phase=ifelse(no_snp==1, 'pat', sv_phase),
    zygosity=ifelse(no_snp==1, 'het', zygosity_snp))%>%
    filter(!is.na(POS))
  
  return(sv_sum)
}

#' parse information for germline het snp that are close to SVs
#'
#' @param het_on_sv an object of class 'data frame'. This object stores the het snp on sv data loaded from `load_data`.
#' @param snp_df an object of class 'data frame'. This object stores the output from `parse_het_snps`.
#'
#' @returns sth
#' @export
#'
parse_snp_on_sv <- function(het_on_sv, snp_df) {
  sv_phase <- het_on_sv %>%
    select(-ALT)%>%
    # get allele count for SNPs
    mutate(snp_id=row_number(),
           AD=gsub(".*:","",bulk),
           len=nchar(gsub("\\d","",AD)),
           DEP=gsub("DP=(\\d+);.*","\\1", INFO),
           onsv_ref=as.integer(gsub(",.*","", AD)),
           onsv_alt=as.integer(ifelse(len==2, gsub(".*,(\\d+),.*","\\1", AD),0)),
           onsv_alt2=as.integer(ifelse(len==2, gsub(".*,","\\1", AD), 0)), # check if this SNP has three alleles, which is filtered out
           a_count=(onsv_ref!=0) + (onsv_alt!=0),
           allele=ifelse(a_count==1 & onsv_ref!=0, REF, "other"))%>% # other means the allele is not reference allele (used to tell mat or pat origin)
    filter(DEP > 0,
           onsv_alt2==0)%>%
    # join in all het SNP
    left_join(snp_df, by = c("CHROM", "POS", "REF"))%>%
    select(-ID)
  
  return(sv_phase)
}

#' extract information
#'
#' @param p_het object of class 'character'. This object stores the file path 
#' to vcf file for heterozygous SNPs.
#' @param p_onsv object of class 'character'. This object stores the file path 
#' to vcf files for heterozygous SNPs that are on SV supporting reads.
#' @param p_sv object of class 'character'. This object stores the file path 
#' to vcf files for SVs.
#' @param p_cnv object of class 'character'. This object stores the file path 
#' to CNV file. 
#' @param flank_del object of class 'numeric'. This object describes the maximum allowed
#' differences in genomic locations for a deletion to be considered as overlapping to a translocation.  
#' @param QUAL_tresh object of class 'numeric'. This object describes the minimum quality score 
#'  allowed to include an SV.
#' @param min_alt object of class 'numeric'. This object describes the minimum 
#' amount of SV supporting read counts to include an SV.
#' @param tumor_only object of class 'Boolean'. This object describes if the SVs were called only using tumor BAM
#'
#' @returns sth
#' @export
#'

extract_info <- function(p_het, p_onsv, p_sv, p_cnv, chr_lst=NULL, flank_del=50, QUAL_tresh=100, min_alt=2, tum_only){
  data      <- load_data(p_het, p_onsv, p_sv, p_cnv, chr=chr_lst, tumor_only=tum_only)
  bnd_info   <- proc_bnd(data$sv, flank_del=50)
  sv_info   <- parse_sv_info(data$sv, bnd_info[[1]], bnd_info[[2]], QUAL_tresh=100, min_alt=2)
  snp_df    <- parse_het_snps(data$het_snp)
  sv_phase  <- parse_snp_on_sv(data$het_on_sv, snp_df)
  return(list(data, sv_info, snp_df, sv_phase))
}

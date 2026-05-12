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
#' @param QUAL_thresh object of class 'numeric'. This object describes the minimum quality score 
#'  allowed to include an SV.
#' @param chr_lst Character vector or \code{NULL}. Restrict analysis to these
#'   chromosomes (e.g. \code{"chr1"}). Default \code{NULL} (all chromosomes).
#' @param min_alt object of class 'numeric'. This object describes the minimum
#' amount of SV supporting read counts to include an SV.
#' @param tum_only Logical. Whether SVs were called from a tumor-only BAM.
#'   Default \code{FALSE}.
#'
#' @return A list of four elements: \code{[[1]]} the raw data list from
#'   \code{\link{load_data}} (named elements \code{het_snp}, \code{het_on_sv},
#'   \code{sv}, \code{cnv}); \code{[[2]]} \code{sv_info} from
#'   \code{\link{parse_sv_info}}; \code{[[3]]} \code{snp_df} from
#'   \code{\link{parse_het_snps}}; \code{[[4]]} \code{sv_phase} from
#'   \code{\link{parse_snp_on_sv}}.
#' @export
#'

extract_info <- function(p_het, p_onsv, p_sv, p_cnv, chr_lst=NULL, flank_del=50, QUAL_thresh=100, min_alt=2, tum_only){
  data      <- load_data(p_het, p_onsv, p_sv, p_cnv, chr=chr_lst, tumor_only=tum_only)
  bnd_info   <- proc_bnd(data$sv, flank_del=flank_del)
  sv_info   <- parse_sv_info(data$sv, bnd_info$bnd, bnd_info$del, QUAL_thresh=QUAL_thresh, min_alt=min_alt)
  snp_df    <- parse_het_snps(data$het_snp)
  sv_phase  <- parse_snp_on_sv(data$het_on_sv, snp_df)
  return(list(data, sv_info, snp_df, sv_phase))
}

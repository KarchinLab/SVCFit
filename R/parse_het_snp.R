#' parse information for germline heterozygous snps
#'
#' @param het_snp an object of class 'data frame'. This object stores the het snp data loaded from `load_data`.
#'
#' @return A data.frame with columns \code{CHROM}, \code{POS}, \code{REF},
#'   \code{ALT}, \code{snp_ref} (reference-allele read count), \code{snp_alt}
#'   (alternate-allele read count), and \code{dep} (total sequencing depth).
#'   Restricted to single-nucleotide SNPs (length-1 REF and ALT alleles).
#' @export
#'
parse_het_snps <- function(het_snp) {
  snp_df = het_snp %>%
    filter(nchar(ALT)==1,
           nchar(REF)==1)%>%
    mutate(snp_ref=as.integer(gsub(".*:(\\d+),(\\d+):.*", "\\1", bulk)),
           snp_alt=as.integer(gsub(".*:(\\d+),(\\d+):.*", "\\2", bulk)),
           dep=as.integer(gsub("^(?:[^:]+:){2}([^:]+).*", "\\1",bulk)))%>%
    select(CHROM, POS, REF,ALT,snp_ref, snp_alt, dep)
  
  return(snp_df)
}

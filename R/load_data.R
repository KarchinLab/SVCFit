#' load input data
#'
#' @param p_het object of class 'character'. This object stores the file path 
#' to vcf file for heterozygous SNPs.
#' @param p_onsv object of class 'character'. This object stores the file path 
#' to vcf files for heterozygous SNPs that are on SV supporting reads.
#' @param p_sv object of class 'character'. This object stores the file path 
#' to vcf files for SVs.
#' @param p_cnv object of class 'character'. This object stores the file path 
#' to CNV file. 
#' @param chr an object of class 'character'. This object describes what chromosome to include for analysis.
#' @param tumor_only object of class 'Boolean'. This object describes if the SVs were called only using tumor BAM
#' 
#' @returns sth
#' @export
#'
load_data <- function(p_het, p_onsv, p_sv, p_cnv, chr=NULL, tumor_only) {

  het_snp <- read.table(p_het, quote="\"")
  colnames(het_snp) = c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','bulk')
  if(!is.null(chr)){
    het_snp <- het_snp %>% filter(CHROM %in% chr)
  }
  het_on_sv <- read.table(p_onsv, quote = '"')
  colnames(het_on_sv) = c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','bulk')
  if(!is.null(chr)){
    het_on_sv <- het_on_sv %>% filter(CHROM %in% chr)
  }
  sv <- read.table(p_sv, quote = '"')
  if(tumor_only){
    colnames(sv) = c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','tumor')
  }else{
    colnames(sv) = c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','normal','tumor')
  }
  if(!is.null(chr)){
    sv <- sv %>% filter(CHROM %in% chr)
  }
  cnv <- read.delim(p_cnv, header = TRUE) %>%
    mutate(lcn.em=ifelse(is.na(lcn.em), 1, lcn.em),
           chrom=paste0("chr",chrom))
  if(!is.null(chr)){
    cnv <- cnv %>% filter(chrom %in% chr)
  }

  return(list(het_snp = het_snp,
              het_on_sv = het_on_sv,
              sv = sv,
              cnv = cnv))
}

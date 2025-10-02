#' load input data
#'
#' @param samp an object of class 'character'. This object describes the sample name.
#' @param exper an object of class 'character'. This object describes the experiment number.
#' @param base_dir an object of class 'character'. This object describes the path where file is stored.
#' @param chr an object of class 'character'. This object describes what chromosome to include for analysis.
#'
#' @returns sth
#' @export
#'
load_data <- function(samp, exper, base_dir = "~/Documents/Karchin_lab/visor/m_benchmark", chr=NULL) {

  het_snp <- read.table(file.path(base_dir, 'SNP', samp, exper, paste0("het_near_sv_",exper, ".vcf")), quote="\"")
  colnames(het_snp) = c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','bulk')
  if(!is.null(chr)){
    het_snp <- het_snp %>% filter(CHROM %in% chr)
  }
  het_on_sv <- read.table(file.path(base_dir, 'SNP', samp, exper, paste0("het_on_sv_",exper,".vcf")), quote = '"')
  colnames(het_on_sv) = c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','bulk')
  if(!is.null(chr)){
    het_on_sv <- het_on_sv %>% filter(CHROM %in% chr)
  }
  sv <- read.table(file.path(base_dir, 'svtyp', samp, exper, paste0('svt_', exper, '.vcf')), quote = '"')
  colnames(sv) = c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','normal','tumor')
  if(!is.null(chr)){
    sv <- sv %>% filter(CHROM %in% chr)
  }
  cnv <- read.delim(file.path(base_dir, 'facet', samp, exper, paste0(exper, '.bed')), header = TRUE) %>%
    mutate(lcn.em=ifelse(is.na(lcn.em), 1, lcn.em),
           chrom=paste0("chr",chrom))
  if(!is.null(chr)){
    cnv <- het_snp %>% filter(chrom %in% chr)
  }

  return(list(het_snp = het_snp,
              het_on_sv = het_on_sv,
              sv = sv,
              cnv = cnv))
}

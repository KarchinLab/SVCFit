#' Load simulation truth data
#'
#' @param truth_path Character. Path to the directory containing BED files of
#'   true SV clone assignments.  For non-overlapping simulations the files must
#'   be named \code{c1.bed}, \code{c2.bed}, \code{c3.bed}; for overlapping
#'   simulations \code{c11.bed}, \code{c22.bed}, \code{c33.bed}.
#' @param overlap Logical. Whether the simulation has SV-CNV overlap.
#'   Default \code{FALSE}.
#'
#' @return A data.frame with columns \code{CHROM}, \code{start}, \code{end},
#'   \code{type}, \code{info}, \code{flank}, \code{samp} (clone label such as
#'   \code{"c1"}), \code{exp}, \code{clone} (\code{"clonal"}, \code{"sub"}, or
#'   \code{"2"}), \code{chr2}, and \code{pos2}.  One row per unique genomic
#'   start position across all simulated clones.
#' @export
#'
load_truth <- function(truth_path, overlap=FALSE){
  
  if(!overlap){
    c1=read.delim(paste0(truth_path,"/c1.bed"), header=FALSE) %>%
      mutate(samp='c1')
    c2=read.delim(paste0(truth_path,"/c2.bed"), header=FALSE) %>%
      mutate(samp='c2')
    c3=read.delim(paste0(truth_path,"/c3.bed"), header=FALSE) %>%
      mutate(samp='c3')
    tmp=rbind(c1,c2,c3)
    colnames(tmp)=c('CHROM','start','end','type','info','flank', 'samp')
  }else{
    c1=read.delim(paste0(truth_path,"/c11.bed"), header=FALSE) %>%
      mutate(samp='c11')
    c2=read.delim(paste0(truth_path,"/c22.bed"), header=FALSE) %>%
      mutate(samp='c22')
    c3=read.delim(paste0(truth_path,"/c33.bed"), header=FALSE) %>%
      mutate(samp='c33')
    tmp=rbind(c1,c2,c3)
    colnames(tmp)=c('CHROM','start','end','type','info','flank','samp')
  }
  truth=tmp %>%
    mutate(CHROM = paste0("chr", sub("^chr", "", CHROM)),
           exp=ifelse(samp %in% c('c1','c2','c3'), 'exp1','other'))%>%
    group_by(exp, start, end)%>%
    mutate(num=n(),
           clone=ifelse(num==1, 'sub','clonal'),
           clone=ifelse(grepl('2', samp) & num==1, "2",clone))%>%
    ungroup()%>%
    distinct(start, .keep_all = T)%>%
    mutate(chr2=gsub('.*(chr\\d+).*', "\\1",info),
           pos2=ifelse(grepl('trans', type),gsub('.*:(\\d+).*', '\\1',info), end),
           pos2=as.integer(pos2))
  
  return(truth)
}

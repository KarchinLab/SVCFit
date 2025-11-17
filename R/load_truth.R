#' load simulation truth data
#'
#' @param exper an object of class 'character'. This object describes the experiment number.
#'
#' @returns sth
#' @export
#'
load_truth <- function(exper){
  #path_lst=list.files("~/Documents/Karchin_lab/visor/m_benchmark/truth", full.names = T)
  if(exper=='exp1'){
    c1=read.delim("~/Documents/Karchin_lab/visor/m_benchmark/truth/c1.bed", header=FALSE) %>%
      mutate(samp='c1')
    c2=read.delim("~/Documents/Karchin_lab/visor/m_benchmark/truth/c2.bed", header=FALSE) %>%
      mutate(samp='c2')
    c3=read.delim("~/Documents/Karchin_lab/visor/m_benchmark/truth/c3.bed", header=FALSE) %>%
      mutate(samp='c3')
    tmp=rbind(c1,c2,c3)
    colnames(tmp)=c('CHROM','start','end','type','info','flank', 'samp')
  }else{
    c1=read.delim("~/Documents/Karchin_lab/visor/m_benchmark/truth/c11.bed", header=FALSE) %>%
      mutate(samp='c11')
    c2=read.delim("~/Documents/Karchin_lab/visor/m_benchmark/truth/c22.bed", header=FALSE) %>%
      mutate(samp='c22')
    c3=read.delim("~/Documents/Karchin_lab/visor/m_benchmark/truth/c33.bed", header=FALSE) %>%
      mutate(samp='c33')
    tmp=rbind(c1,c2,c3)
    colnames(tmp)=c('CHROM','start','end','type','info','flank','samp')
  }
  truth=tmp %>%
    mutate(exp=ifelse(samp %in% c('c1','c2','c3'), 'exp1','other'))%>%
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

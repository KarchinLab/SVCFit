#' attach truth of simulation to the output of SVCFit
#'
#' @param svcf_out an object of class 'data frame'. This object stores the output from `calc_svcf`.
#' @param truth an object of class 'data frame'. This object stores the output from `load_truth`.
#'
#' @returns sth
#' @export
#'
attach_truth <- function(svcf_out, truth){
  svcf_truth=svcf_out %>%
    rowwise()%>%
    mutate(type=any((CHROM==truth$CHROM & abs(POS-truth$start)<100) | (chr2==truth$chr2 & abs(END-truth$pos2)<100)))%>%
    filter(type)%>%
    mutate(row=which((CHROM==truth$CHROM & abs(POS-truth$start)<100) | (chr2==truth$chr2 & abs(END-truth$pos2)<100)),
           clone=truth$clone[row])
  return(svcf_truth)
}
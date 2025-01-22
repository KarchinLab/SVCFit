
#' Calculate structural variant cellular fraction from vcf files
#'
#' @param vcf_path an object of class "Character". This variable is the path to vcf files
#' @param simulation an object of class "Boolean". This variable indicates if the vcf files are from a simulation
#' @param overlap an object of class "Boolean". This variable indicates if a structural variants should be filtered based on coordinates overlap
#' @param tolerance an object of class "Integer". This variable sets the threshold coordinates overlaps for two SVs
#' @param window an object of class "Integer". This variable sets how many SVs should be compared at the same time
#' @param multiple an object of class "Boolean". This variable indicates whether the sample has multiple clone (only used for simulated data)
#' @param tumor_only an object of class "Boolean". This variable indicates whether the vcf files were created under tumor-only mode
#' @param truth_path an object of class "Character". This variable is the path to the ground truth designed for simulation
#' @param length_filter an object of class "integer". This variable set lowest threshold on the size of a structural variants
#'
#' @return an annotated vcf file
#' @export
#' @import tidyverse
#' @examples example
SVCFit <- function(vcf_path, simulation=FALSE, overlap=TRUE, tolerance=6, window=100, multiple=FALSE, tumor_only=FALSE, truth_path=NULL, length_filter=0){

  # load vcf files
  vcf=extract_info(vcf_path, tumor_only)

  if(simulation){
    ground_truth=ifelse(is.null(truth_path), "/Users/lyz928/Karchin Lab Dropbox/YunZhou Liu/SVCFit-2024/script/sv/visor/input/", truth_path)
    p <- list.files(ground_truth, pattern = "^c.*.bed", full.names = T)
    truth = read_clone(p)
  }else{
    truth=data.frame()
  }

  # filer SV based on simulation ground truth and overlaps
  if(simulation & overlap){
    over=check_overlap(vcf, vcf)
    checked=check_overlap(over, truth)
  }else if(simulation & !overlap){
    checked=check_overlap(vcf, truth)
  }else if((!simulation) & overlap){
    checked=check_overlap(vcf, vcf)
  }else{
    checked=vcf
  }

  #calculate svcf
  if(tumor_only==FALSE){
    output <- checked %>%
      dplyr::filter(!classification%in%c("INS","BND"))%>%
      dplyr::mutate(vaf = alt/(alt+0.5*ref),
             tcn= ifelse(classification =="DUP", (4*alt+2*ref)/ref, 2),
             inferred_icn = ifelse(classification=="DUP", round(tcn*vaf+2),2),
             svcf = ifelse(inferred_icn<4, tcn*vaf, tcn*vaf/(inferred_icn-2)))%>%
      dplyr::select(sample, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,normal,tumor,classification,pos2,vaf, tcn, inferred_icn, svcf)
  }else{
    output <- checked %>%
      dplyr::filter(!classification%in%c("INS","BND"))%>%
      dplyr::mutate(vaf = alt/(alt+0.5*ref),
             tcn= ifelse(classification =="DUP", (4*alt+2*ref)/ref, 2),
             inferred_icn = ifelse(classification=="DUP", round(tcn*vaf+2),2),
             svcf = ifelse(inferred_icn<4, tcn*vaf, tcn*vaf/(inferred_icn-2)))%>%
      dplyr::select(sample, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,tumor,classification,pos2,vaf, tcn, inferred_icn, svcf)
  }

  if(multiple&simulation&nrow(truth)>0){ # make sure truth is loaded
    output = attach_clone(output, truth)
  }

  return(output)
}

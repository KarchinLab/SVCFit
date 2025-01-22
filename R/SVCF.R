
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
#' @import magrittr
#' @import dplyr
#' @examples example
SVCF <- function(vcf_path, simulation=FALSE, overlap=TRUE, tolerance=6, window=100, multiple=FALSE, tumor_only=FALSE, truth_path=NULL, length_filter=0){

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
      dplyr::filter(!.data$classification%in%c("INS","BND"))%>%
      dplyr::mutate(vaf = .data$alt/(.data$alt+0.5*.data$ref),
             tcn= ifelse(.data$classification =="DUP", (4*.data$alt+2*.data$ref)/.data$ref, 2),
             inferred_icn = ifelse(.data$classification=="DUP", round(.data$tcn*.data$vaf+2),2),
             svcf = ifelse(.data$inferred_icn<4, .data$tcn*.data$vaf, .data$tcn*.data$vaf/(.data$inferred_icn-2)))%>%
      dplyr::select(.data$sample, .data$CHROM,.data$POS,.data$ID,.data$REF,.data$ALT,.data$QUAL,.data$FILTER,.data$INFO,.data$FORMAT,.data$normal,.data$tumor,.data$classification,.data$pos2,.data$vaf, .data$tcn, .data$inferred_icn, .data$svcf)
  }else{
    output <- checked %>%
      dplyr::filter(!.data$classification%in%c("INS","BND"))%>%
      dplyr::mutate(vaf = .data$alt/(.data$alt+0.5*.data$ref),
             tcn= ifelse(.data$classification =="DUP", (4*.data$alt+2*.data$ref)/.data$ref, 2),
             inferred_icn = ifelse(.data$classification=="DUP", round(.data$tcn*.data$vaf+2),2),
             svcf = ifelse(.data$inferred_icn<4, .data$tcn*.data$vaf, .data$tcn*.data$vaf/(.data$inferred_icn-2)))%>%
      dplyr::select(.data$sample, .data$CHROM,.data$POS,.data$ID,.data$REF,.data$ALT,.data$QUAL,.data$FILTER,.data$INFO,.data$FORMAT,.data$tumor,.data$classification,.data$pos2,.data$vaf, .data$tcn, .data$inferred_icn, .data$svcf)
  }

  if(multiple&simulation&nrow(truth)>0){ # make sure truth is loaded
    output = attach_clone(output, truth)
  }

  return(output)
}

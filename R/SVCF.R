
#' Calculate structural variant cellular fraction from vcf files
#'
#' @param vcf_path an object of class "Character". This variable is the path to
#' vcf files
#' @param overlap an object of class "Boolean". This variable indicates if a
#' structural variants should be filtered based on coordinates overlap
#' @param tolerance an object of class "Integer". This variable sets the threshold
#' coordinates overlaps for two SVs
#' @param window an object of class "Integer". This variable sets how many SVs
#' should be compared at the same time
#' @param multiple an object of class "Boolean". This variable indicates whether
#' the sample has multiple clone (only used for simulated data)
#' @param tumor_only an object of class "Boolean". This variable indicates whether
#' the vcf files were created under tumor-only mode
#' @param truth_path an onject of class "Character". This variable is a path to
#' bed files storing true structural variants information with clonal assignment.
#' Each bed file should be named like "c1.bed, c2.bed" etc. Structural variants
#' should be saved in a seperate bed file if they belong to different (sub)clone.
#' @param mode an onject of class "Character". This variable describe how true
#' clonal information is saved. In "inherited" mode, bed files for all children
#' clone contains all ancestral structural variants of their parents. In "distinct"
#' mode, children clones don't contain any ancestral structural variants.
#' @param length_threshold an object of class "integer". This variable set lowest
#' threshold on the size of a structural variants
#'
#' @return an annotated vcf file
#' @export
#' @import tidyverse

SVCF <- function(vcf_path=NULL, overlap=TRUE, tolerance=6, window=100, multiple=FALSE,
                 tumor_only=FALSE, truth_path=NULL, mode="inherited", length_threshold=0){

  if(is.null(vcf_path)){
    stop("No vcf file path is provided.\n")
  }
  # load vcf files
  vcf=extract_info(vcf_path, tumor_only, length_threshold)

  if(!is.null(truth_path)){
    #"/Users/lyz928/Karchin Lab Dropbox/YunZhou Liu/SVCFit-2024/script/sv/visor/input/"
    p <- list.files(truth_path, pattern = "^c.*.bed", full.names = T)
    truth = read_clone(p, mode)
  }else{
    truth=data.frame()
  }

  # filer SV based on simulation ground truth and overlaps
  if(overlap & !is.null(truth_path)){

    over=check_overlap(vcf, vcf, tolerance, window)
    checked=check_overlap(over, truth, tolerance, window)

  }else if(!overlap & !is.null(truth_path)){

    checked=check_overlap(vcf, truth, tolerance, window)

  }else if(overlap & is.null(truth_path)){

    checked=check_overlap(vcf, vcf, tolerance, window)

  }else{
    checked=vcf
  }

  #calculate svcf
  output=calculate_svcf(checked, tumor_only)

  if(multiple & nrow(truth)>0){ # make sure truth is loaded
    output = attach_clone(output, truth, tolerance)
  }

  return(output)
}

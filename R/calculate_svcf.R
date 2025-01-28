#' Title
#'
#' @param input an object of class "dataframe". This object stores the set of information for structural variants used to calculate SVCF
#' @param tumor_only an object of class "Boolean". This variable indicates whether the vcf files were created under tumor-only mode
#'
#' @return an annotated vcf file
#' @export
#'
calculate_svcf <- function(input, tumor_only=FALSE){
  #calculate svcf
  if(tumor_only==FALSE){
    output <- input %>%
      dplyr::filter(!.data$classification%in%c("INS","BND"))%>%
      dplyr::mutate(vaf = .data$alt/(.data$alt+0.5*.data$ref),
                    tcn= ifelse(.data$classification =="DUP", (4*.data$alt+2*.data$ref)/.data$ref, 2),
                    inferred_icn = ifelse(.data$classification=="DUP", round(.data$tcn*.data$vaf+2),2),
                    svcf = ifelse(.data$inferred_icn<4, .data$tcn*.data$vaf, .data$tcn*.data$vaf/(.data$inferred_icn-2)))%>%
      dplyr::select(sample, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,normal,tumor,classification,pos2,vaf,tcn,inferred_icn,svcf)
  }else{
    output <- input %>%
      dplyr::filter(!.data$classification%in%c("INS","BND"))%>%
      dplyr::mutate(vaf = .data$alt/(.data$alt+0.5*.data$ref),
                    tcn= ifelse(.data$classification =="DUP", (4*.data$alt+2*.data$ref)/.data$ref, 2),
                    inferred_icn = ifelse(.data$classification=="DUP", round(.data$tcn*.data$vaf+2),2),
                    svcf = ifelse(.data$inferred_icn<4, .data$tcn*.data$vaf, .data$tcn*.data$vaf/(.data$inferred_icn-2)))%>%
      dplyr::select(sample, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,tumor,classification,pos2,vaf,tcn,inferred_icn,svcf)
  }

  return(output)
}

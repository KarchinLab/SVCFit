#' Title
#'
#' @param input an object of class "dataframe". This object stores the set of
#' information for structural variants used to calculate SVCF
#' @param tumor_only an object of class "Boolean". This variable indicates whether
#' the vcf files were created under tumor-only mode
#'
#' @return an annotated vcf file
#'
#' @import tidyverse
#' @export
#'
calculate_svcf <- function(input, tumor_only=FALSE){
  #calculate svcf
  if(tumor_only==FALSE){
    output <- input %>%
      filter(!classification%in%c("INS","BND"))%>%
      mutate(vaf = alt/(alt+0.5*ref),
                    Rbar= ifelse(classification =="DUP",
                                (4*alt+2*ref)/ref, 2),
                    r = ifelse(classification=="DUP",
                                          round(Rbar*vaf+2),2),
                    svcf = ifelse(r<4,
                                  Rbar*vaf, Rbar*vaf/(r-2)))%>%
      select(sample, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,normal,
                    tumor,classification,pos2,vaf,Rbar,r,svcf)
  }else{
    output <- input %>%
      filter(!classification%in%c("INS","BND"))%>%
      mutate(vaf = alt/(alt+0.5*ref),
                    Rbar= ifelse(classification =="DUP",
                                (4*alt+2*ref)/ref, 2),
                    r = ifelse(classification=="DUP",
                                          round(Rbar*vaf+2),2),
                    svcf = ifelse(r<4,
                                  Rbar*vaf, Rbar*vaf/(r-2)))%>%
      select(sample, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,tumor,
             classification,pos2,vaf,Rbar,r,svcf)
  }

  return(output)
}

#' extract_info is used to get the information from vcf files to calculate structural variant cellular fraction
#'
#' @param path an object of class "Character". This variable is the path to vcf files
#' @param tumor_only an object of class "Boolean". This variable indicates whether
#' the vcf files were created under tumor-only mode
#' @param length_threshold an object of class "integer". This variable set lowest
#' threshold on the size of a structural variants
#'
#' @return a dataframe with read information for each structural variant
#' @export
#' @import tidyverse
#' @importFrom utils read.table
#' @examples example
extract_info <- function(path, tumor_only=FALSE, length_threshold=0){
  svt <- read.table(path, quote="\"")
  sample_id <- gsub(".*(c.*).vcf","\\1",path)

  if(tumor_only){
    value=ncol(svt)!=length(c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","tumor"))
    if(value){
      stop("Number of vcf input column doesn't match the column name length, consider change tumor_only option.\n")
    }
    colnames(svt) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","tumor")
  }else{
    value=ncol(svt)!=length(c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","normal","tumor"))
    if(value){
      stop("Number of vcf input column doesn't match the column name length, consider change tumor_only option.\n")
    }
    colnames(svt) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","normal","tumor")
  }


  if(tumor_only){
    out <- svt %>%
      mutate(chr2 = gsub(".*(chr.).*","\\1", ALT),
             chr2 = ifelse(grepl("chr", chr2), chr2, CHROM),
             pos2 = gsub(".*:(\\d+).*","\\1", ALT),
             pos2 = ifelse(grepl("\\d+", pos2), pos2, gsub("END=(\\d+).*","\\1",INFO)),
             pos2 = as.integer(pos2),
             classification = gsub(".*SVTYPE=(\\w+).*","\\1", INFO),
             length = abs(POS-pos2),
             ID = gsub(":\\d$","",ID),
             ## here, sum the spanning and split read
             PR = gsub(":.*$", "", tumor),
             SR = ifelse(FORMAT == "PR", "0,0", gsub("^.*:", "", tumor)),
             rpr = gsub(",\\d+","", PR),
             apr = gsub("\\d+,","", PR),
             rsr = gsub(",\\d+","", SR),
             asr = gsub("\\d+,","", SR),
             ref = as.numeric(rpr)+as.numeric(rsr),
             alt = as.numeric(apr)+as.numeric(asr),
             sample = sample_id,
             row = as.integer(row_number()),
             iid=paste0(CHROM,"_",POS))%>%
      filter(FILTER == "PASS",
                    alt > 2,
                    length>length_threshold)
  }else{
  out <- svt %>%
    mutate(chr2 = gsub(".*(chr.).*","\\1", ALT),
           chr2 = ifelse(grepl("chr", chr2), chr2, CHROM),
           pos2 = gsub(".*:(\\d+).*","\\1", ALT),
           pos2 = ifelse(grepl("\\d+", pos2), pos2, gsub("END=(\\d+).*","\\1",INFO)),
           pos2 = as.integer(pos2),
           length = abs(POS-pos2),
           classification = gsub(".*SVTYPE=(\\w+).*","\\1", INFO),
           ID = gsub(":\\d$","",ID),
           ## here, sum the spanning and split read
           PR = gsub(":.*$", "", tumor),
           SR = ifelse(FORMAT == "PR", "0,0", gsub("^.*:", "", tumor)),
           NPR = gsub(":.*$", "", normal),
           NSR = ifelse(FORMAT == "PR", "0,0", gsub("^.*:", "", normal)),
           anpr = gsub("\\d+,","", NPR),
           ansr = gsub("\\d+,","", NSR),
           rpr = gsub(",\\d+","", PR),
           apr = gsub("\\d+,","", PR),
           rsr = gsub(",\\d+","", SR),
           asr = gsub("\\d+,","", SR),
           ref = as.numeric(rpr)+as.numeric(rsr),
           alt = as.numeric(apr)+as.numeric(asr),
           sample = sample_id,
           row = as.integer(row_number()),
           iid=paste0(CHROM,"_",POS))%>%
    filter(FILTER == "PASS",
                  alt > 2,
                  ansr == 0,
                  anpr == 0,
                  length>length_threshold)
  }
  return(out)
}

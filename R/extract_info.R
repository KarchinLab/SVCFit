#' extract_info is used to get the information from vcf files to calculate structural variant cellular fraction
#'
#' @param path an object of class "Character". This variable is the path to vcf files
#' @param tumor_only an object of class "Boolean". This variable indicates whether the vcf files were created under tumor-only mode
#' @param length_threshold an object of class "integer". This variable set lowest threshold on the size of a structural variants
#'
#' @return a dataframe with read information for each structural variant
#' @export
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @importFrom utils read.table
#' @examples example
extract_info <- function(path, tumor_only=FALSE, length_threshold=0){
  svt <- read.table(path, quote="\"")
  sample_id <- gsub(".*(c.*).vcf","\\1",path)

  if(tumor_only){
    value=ncol(svt)!=length(c("CHROM",  "POS",     "ID",      "REF",     "ALT",     "QUAL",    "FILTER",  "INFO",    "FORMAT", "tumor"))
    if(value){
      stop("Number of vcf input column doesn't match the column name length, consider change tumor_only option.\n")
    }
    colnames(svt) <- c("CHROM",  "POS",     "ID",      "REF",     "ALT",     "QUAL",    "FILTER",  "INFO",    "FORMAT", "tumor")
  }else{
    value=ncol(svt)!=length(c("CHROM",  "POS",     "ID",      "REF",     "ALT",     "QUAL",    "FILTER",  "INFO",    "FORMAT","normal", "tumor"))
    if(value){
      stop("Number of vcf input column doesn't match the column name length, consider change tumor_only option.\n")
    }
    colnames(svt) <- c("CHROM",  "POS",     "ID",      "REF",     "ALT",     "QUAL",    "FILTER",  "INFO",    "FORMAT", "normal", "tumor")
  }


  if(tumor_only){
    out <- svt %>%
      dplyr::mutate(chr2 = gsub(".*(chr.).*","\\1", .data$ALT),
             chr2 = ifelse(grepl("chr", .data$chr2), .data$chr2, .data$CHROM),
             pos2 = gsub(".*:(\\d+).*","\\1", .data$ALT),
             pos2 = ifelse(grepl("\\d+", .data$pos2), .data$pos2, gsub("END=(\\d+).*","\\1",.data$INFO)),
             pos2 = as.integer(.data$pos2),
             classification = gsub(".*SVTYPE=(\\w+).*","\\1", .data$INFO),
             length = abs(.data$POS-.data$pos2),
             ID = gsub(":\\d$","",.data$ID),
             ## here, sum the spanning and split read
             PR = gsub(":.*$", "", .data$tumor),
             SR = ifelse(.data$FORMAT == "PR", "0,0", gsub("^.*:", "", .data$tumor)),
             rpr = gsub(",\\d+","", .data$PR),
             apr = gsub("\\d+,","", .data$PR),
             rsr = gsub(",\\d+","", .data$SR),
             asr = gsub("\\d+,","", .data$SR),
             ref = as.numeric(.data$rpr)+as.numeric(.data$rsr),
             alt = as.numeric(.data$apr)+as.numeric(.data$asr),
             sample = sample_id,
             row = as.integer(row_number()),
             iid=paste0(.data$CHROM,"_",.data$POS))%>%
      dplyr::filter(.data$FILTER == "PASS",
                    .data$alt > 2,
                    .data$length>length_threshold)
  }else{
  out <- svt %>%
    dplyr::mutate(chr2 = gsub(".*(chr.).*","\\1", .data$ALT),
           chr2 = ifelse(grepl("chr", .data$chr2), .data$chr2, .data$CHROM),
           pos2 = gsub(".*:(\\d+).*","\\1", .data$ALT),
           pos2 = ifelse(grepl("\\d+", .data$pos2), .data$pos2, gsub("END=(\\d+).*","\\1",.data$INFO)),
           pos2 = as.integer(.data$pos2),
           length = abs(.data$POS-.data$pos2),
           classification = gsub(".*SVTYPE=(\\w+).*","\\1", .data$INFO),
           ID = gsub(":\\d$","",.data$ID),
           ## here, sum the spanning and split read
           PR = gsub(":.*$", "", .data$tumor),
           SR = ifelse(.data$FORMAT == "PR", "0,0", gsub("^.*:", "", .data$tumor)),
           NPR = gsub(":.*$", "", .data$normal),
           NSR = ifelse(.data$FORMAT == "PR", "0,0", gsub("^.*:", "", .data$normal)),
           anpr = gsub("\\d+,","", .data$NPR),
           ansr = gsub("\\d+,","", .data$NSR),
           rpr = gsub(",\\d+","", .data$PR),
           apr = gsub("\\d+,","", .data$PR),
           rsr = gsub(",\\d+","", .data$SR),
           asr = gsub("\\d+,","", .data$SR),
           ref = as.numeric(.data$rpr)+as.numeric(.data$rsr),
           alt = as.numeric(.data$apr)+as.numeric(.data$asr),
           sample = sample_id,
           row = as.integer(row_number()),
           iid=paste0(.data$CHROM,"_",.data$POS))%>%
    dplyr::filter(.data$FILTER == "PASS",
                  .data$alt > 2,
                  .data$ansr == 0,
                  .data$anpr == 0,
                  .data$length>length_threshold)
  }
  return(out)
}

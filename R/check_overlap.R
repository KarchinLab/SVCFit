#' check_overlap compare two sets of coordinates of structural variants and merge the structural variants with coordinates within tolerance.
#'
#' @param dat an object of class "dataframe". This object stores the first set of structural variants to be compared
#' @param compare an object of class "dataframe". This object stores the second set of structural variants for comparison
#' @param tolerance an object of class "integer". This variable sets the threshold coordinates overlaps for two SVs
#' @param window an object of class "Integer". This variable sets how many SVs should be compared at the same time
#'
#' @return an filtered annotated vcf file
#' @export
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @examples example
check_overlap <- function(dat, compare, tolerance=6, window=1000){

  result = list()

  for (i in seq(1, nrow(dat), by = window)) {
    # first start with a chunk of tmp
    chunk1 <- dat$iid[i:min(i + window - 1, nrow(dat))]
    chr1=gsub("_.*$","",chunk1)
    pos1=as.numeric(gsub("^.*_","",chunk1))
    for (j in seq(1, nrow(compare), by = window)) {
      # the compare the chunk of tmp with every tmp_clone
      chunk2 <- compare$iid[j:min(j + window - 1, nrow(compare))]
      chr2=gsub("_.*$","",chunk2)
      pos2=as.numeric(gsub("^.*_","",chunk2))
      # find pos difference
      diff_pos <- abs(outer(pos1, pos2, "-"))
      diff_pos[diff_pos > tolerance]=-10
      # check whether on the same chromosome
      idx <- which(diff_pos != -10, arr.ind = TRUE) #find which positions overlap
      filter <- chr1[idx[,1]]==chr2[idx[,2]]
      sub_idx = idx[filter, ,drop=FALSE] # index of overlapping entry
      # fix the row index, to better align back to the data
      sub_idx[,1]=idx[,1]+(i-1)
      sub_idx[,2]=idx[,2]+(j-1)
      result = append(result, list(sub_idx))
    }
  }
  rr = as.data.frame(do.call(rbind,result))
  out <- dat %>%
    dplyr::mutate(row=as.integer(row_number()))%>%
    dplyr::filter(.data$row %in% rr$row)%>% # only kept ones that overlaps: self-compare will retain everything, truth-compare will keep overlap
    dplyr::left_join(rr)%>%
    dplyr::group_by(col)%>%
    ## use mean for alt and ref
    dplyr::mutate(alt=round(mean(.data$alt)),
           ref=round(mean(.data$ref)))%>%
    dplyr::ungroup()%>%
    dplyr::distinct(.data$iid, .keep_all = T)%>%
    dplyr::select(-c("col", "row"))

  return(out)
}

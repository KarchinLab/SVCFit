#' check_overlap compare two sets of coordinates of structural variants and merge the structural variants with coordinates within tolerance.
#'
#' @param dat an object of class "dataframe". This object stores the first set of
#' structural variants to be compared
#' @param compare an object of class "dataframe". This object stores the second
#' set of structural variants for comparison
#' @param tolerance an object of class "integer". This variable sets the threshold
#' coordinates overlaps for two SVs
#' @param window an object of class "Integer". This variable sets how many SVs
#' should be compared at the same time
#'
#' @return an filtered annotated vcf file
#' @export
#' @import tidyverse
#' @import igraph

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
  ## in rr, when row and column doesn't match there is an overlap
  rr = as.data.frame(do.call(rbind,result))%>%
    mutate(Row=ifelse(row>col, col, row),
           Col=ifelse(col>row, col, row),
           type=paste0(Row,"-",Col),
           row=Row,
           col=Col)%>%
    distinct(type,.keep_all = T)%>%
    arrange(Row)%>%
    select(row, col)
  ## here, group overlapping with common variable
  g <- graph_from_edgelist(as.matrix(rr), directed = FALSE)
  comp <- components(g)
  component_df <- data.frame(node = unique(c(rr$row, rr$col)), group = comp$membership)
  rr= rr %>% left_join(component_df, by = c("row" = "node"))

  out <- dat %>%
    mutate(row=as.integer(row_number()))%>%
    # only kept ones that overlaps: self-compare will retain everything, truth-compare will keep overlap
    filter(row %in% rr$row)%>%
    left_join(rr)%>%
    group_by(group)%>%
    ## use mean for alt and ref
    mutate(alt=round(mean(alt)),
           ref=round(mean(ref)))%>%
    ungroup()%>%
    distinct(group, .keep_all = T)%>%
    select(-c("col", "row","group"))

  return(out)
}

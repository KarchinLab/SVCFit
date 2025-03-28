#' attach_clone is used to assign each structural variant to a clone designed in simulation. Only applicable if the ground truth is known and the simulation has multiple clones.
#'
#' @param dat an object of class "dataframe". This object is a annotated vcf file
#' storing structural variants to be compared
#' @param truth an object of class "dataframe". This object stores the clone assignment
#' for each structural variants designed in a simulation
#' @param tolerance an object of class "integer". This variable sets the threshold
#' of coordinates overlaps for two SVs
#'
#' @return add one more column "clone_num" to the dat
#' @export
#' @import tidyverse

attach_clone <- function(dat, truth, tolerance=6){
  out <- dat %>%
    rowwise()%>%
    mutate(type = any(abs(POS-truth$pos1)<=tolerance))%>%
    filter(type == TRUE)%>%
    mutate(clone_id = as.integer(which(abs(POS-truth$pos1)<=tolerance)),
           clone_num = truth$id[which(abs(POS-truth$pos1)<=tolerance)])%>%
    ungroup()%>%
    select(-type, -clone_id)
  return(out)
}

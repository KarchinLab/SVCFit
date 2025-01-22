
#' This function reads the ground truth designed for a simulation
#'
#' @param fp an onject of class "Character". This variable is a path to bed files storing structural variants information
#'
#' @return a dataframe
#' @export
#' @import tidyverse
#' @examples example
read_clone <- function(fp){
  # functions to read in true SV locations from bed
  read <- function(path){
    clone <- paste0("clone",gsub(".*/\\w(\\d+).bed","\\1",path))
    tmp <- read.delim(path, header=FALSE)%>%
      dplyr::mutate(id=clone)
    colnames(tmp) <- c("chr","pos1","pos2","sv","info","flank","id")
    return(tmp)
  }

  # process clone
clone <- lapply(fp,function(x) read(x))%>%
  do.call(rbind,.)%>%
  dplyr::group_by(pos1) %>%
  dplyr::mutate(n = n(),
         id = ifelse(n == 3, "clone3", id))%>%
  dplyr::ungroup()%>%
  dplyr::group_by(id, pos1)%>%
  dplyr::distinct(pos1,.keep_all = T)%>%
  dplyr::arrange(pos1)%>%
  dplyr::ungroup()%>%
  as.data.frame() %>%
  dplyr::mutate(iid = paste0(chr,"_",pos1))

return(clone)
}

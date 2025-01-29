
#' This function reads the ground truth designed for a simulation
#'
#' @param truth_path an onject of class "Character". This variable is a path to
#' bed files storing true structural variants information with clonal assignment.
#' Each bed file should be named as "clone" + "number". Structural variants should
#' be saved in a seperate bed file if they belong to different (sub)clone.
#' @param mode an onject of class "Character". This variable describe how true
#' clonal information is saved. In "heritage" mode, bed files for all children
#' clone contains all ancestral structural variants of their parents. In "separate"
#' mode, children clones don't contain any ancestral structural variants.
#'
#' @return a dataframe
#' @export
#' @import tidyverse
#' @importFrom utils read.delim

read_clone <- function(truth_path, mode="heritage"){
  # functions to read in true SV locations from bed
  read <- function(path){
    clone <- paste0("clone",gsub(".*/\\w(\\d+).bed","\\1",path))
    tmp <- read.delim(path, header=FALSE)%>%
      dplyr::mutate(id=clone)
    colnames(tmp) <- c("chr","pos1","pos2","sv","info","flank","id")
    return(tmp)
  }

  # process clone
  if(mode =="separate"){
    clone <- lapply(truth_path,function(x) read(x))
    out = do.call(rbind, clone)%>%
      arrange(.data$pos1)

  }else if (mode == "heritage"){
    clone <- lapply(truth_path,function(x) read(x))
    out = do.call(rbind, clone)%>%
      group_by(.data$pos1) %>%
      mutate(num = n(),
             id = ifelse(.data$num != 1, paste0("clone", num), .data$id))%>%
      ungroup()%>%
      group_by(.data$id, .data$pos1)%>%
      distinct(.data$pos1,.keep_all = T)%>%
      arrange(.data$pos1)%>%
      ungroup()%>%
      as.data.frame() %>%
      mutate(iid = paste0(.data$chr,"_",.data$pos1))
  }
  return(out)
}

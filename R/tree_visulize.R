library(igraph)
library(dplyr)
library(viridis)

# -------------------------------------------------------------------------
# VISUALIZATION
# -------------------------------------------------------------------------

plotTree <- function(edges, palette=viridis::viridis) {
  plotGraph(edgesToAmLong(edges), colorScheme(edges, palette))
}

colorScheme <- function(edges, palette=viridis::viridis) {
  v_sorted = sort(unique(c(edges$parent, edges$child)))
  v_sorted = c(sort(as.integer(v_sorted[!v_sorted=='G'])), "G")
  colors <- c(palette(length(v_sorted)-1), "white")
  v_color <- dplyr::tibble(v_sorted, colors)
  return(v_color)
}

plotGraph <- function(am.long, v_color){
  am.long <- dplyr::mutate(am.long, child = as.numeric(am.long$child)) %>%
    dplyr::arrange(parent, child)
  am.long <- dplyr::mutate(am.long, child = as.character(am.long$child))
  
  am <- toWide(am.long)
  rownames(am) <- c("G", colnames(am))
  am <- cbind(G=0, am) 
  colnames(am) <- rownames(am)
  am[is.na(am)] <- 0
  
  ig <- igraph::graph_from_adjacency_matrix(am, mode = "directed", weighted = TRUE,
                                            diag = FALSE, add.row = TRUE) 
  igraph::V(ig)$color <- as.list(v_color %>% dplyr::arrange(match(v_sorted, names(igraph::V(ig)))) %>% dplyr::select(colors))$colors
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
                      vertex.size=24, vertex.frame.color = "#000000", vertex.label.cex = 1.5,
                      vertex.label.family = "Helvetica", vertex.label.color = "#000000",
                      edge.arrow.size = 0.5, edge.arrow.width = 2)
}


#' Calculate proportions of subclones in each sample (assumes CCFs comply with lineage precedence and sum condition)
#'
#' @param w_mat Matrix of CCF estimates (from \code{estimateCCFs})
#' @param tree_edges Tibble of tree edges with columns edge, parent, and child
#' @export
calcSubcloneProportions <- function(w_mat, tree_edges) {
  K <- nrow(w_mat)
  S <- ncol(w_mat)
  subclone_props <- matrix(NA, nrow = K, ncol = S)
  
  for (i in seq_len(nrow(w_mat))) {
    children <- tree_edges %>%
      filter(parent == as.character(i)) %>%
      pull(child) %>%
      as.numeric()
    
    if (length(children) == 1) {
      children_ccfs <- w_mat[children, ]
    } else if (length(children) > 1) {
      children_ccfs <- w_mat[children, ,drop=FALSE] %>%
        colSums
    } else {
      children_ccfs <- rep(0, ncol(w_mat))
    }
    
    subclone_props[i, ] <- w_mat[i, ] - children_ccfs
  }
  
  # normalize subclone_props matrix so the props add up to 1
  subclone_props[subclone_props < 0] = 0
  subclone_props = round(t(t(subclone_props) / colSums(subclone_props)),digit = 3)
  
  return(subclone_props)
}

#' Plot pie charts for subclone proportions in each sample
#'
#' @param subclone_props matrix of subclone proportions (returned from \code{calcSubcloneProportions})
#' @param sample_names (Optional) Vector of sample names. Should be in the order of columns of subclone_props
#' @export
plotSubclonePie <- function(subclone_props, palette=viridis::viridis, sample_names = NULL, title_size=16, legend_size=10) {
  if (is.null(sample_names)) sample_names <- paste0("Sample ", 1:ncol(subclone_props))
  props_tb <- subclone_props %>%
    magrittr::set_colnames(sample_names) %>%
    as_tibble() %>%
    mutate(Subclone = factor(paste0("Clone ", 1:nrow(.)),
                             levels = paste0("Clone ", 1:nrow(.)))) %>%
    pivot_longer(cols = sample_names,
                 names_to = "Sample",
                 values_to = "Proportion")
  
  clone_colors <- palette(nrow(subclone_props))
  ggplot(props_tb, aes(x="", y=Proportion, fill = Subclone)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = clone_colors, drop = F) +
    theme_void() +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size=legend_size), legend.title = element_text(size=legend_size)) + 
    facet_wrap(~Sample) +
    theme(strip.text.x = element_text(size=title_size))
  
}


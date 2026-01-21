
#' @export
build_tree <- function(clones, lineage_precedence_thresh=0.2, sum_filter_thresh=0.2){
  mcf_mat=clones[,3:4]
  mcf_mat <- as.matrix(mcf_mat)
  rownames(mcf_mat) <- clones$cluster_num
  graph_G_pre <- prepareGraph(mcf_mat, lineage_precedence_thresh)
  graph_G <- filterEdgesBasedOnCCFs(graph_G_pre, mcf_mat, thresh = lineage_precedence_thresh)
  graph_G <- prune(graph_G, mcf_mat) %>% arrange(child)
  graph_G <- assign("graph_G", graph_G, envir = .GlobalEnv)
  enumerateSpanningTreesModified(graph_G, mcf_mat, sum_filter_thresh = sum_filter_thresh)
  scores <- calcTreeScores_n(mcf_mat, all_spanning_trees)
  best_tree <- all_spanning_trees[[which.max(scores)]]
  tree <- plotTree(best_tree)
  return(list(best_tree, mcf_mat, tree))
}

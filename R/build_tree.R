
#build_tree <- function(clones, lineage_precedence_thresh=0.2, sum_filter_thresh=0.2){
#  mcf_mat <- as.matrix(clones[, c("f_pre", "f_day85")])
#  rownames(mcf_mat) <- clones$cluster_num
#  graph_G_pre <- prepareGraph(mcf_mat, lineage_precedence_thresh)
#  graph_G <- filterEdgesBasedOnCCFs(graph_G_pre, mcf_mat, thresh = lineage_precedence_thresh)
#  graph_G <- prune(graph_G, mcf_mat) %>% arrange(child)
#  graph_G <- assign("graph_G", graph_G, envir = .GlobalEnv)
#  enumerateSpanningTreesModified(graph_G, mcf_mat, sum_filter_thresh = sum_filter_thresh)
#  if (length(all_spanning_trees) == 0) return(NULL)
#  scores <- calcTreeScores_n(mcf_mat, all_spanning_trees)
#  best_tree <- all_spanning_trees[[which.max(scores)]]
#  tree <- plotTree(best_tree)
#  return(list(best_tree, mcf_mat, tree))
#}


#' Build the best-scoring phylogenetic spanning tree from clone CCFs
#'
#' Constructs a directed graph of possible parent-child clone relationships
#' based on Cancer Cell Fraction (CCF) values across two time points, enumerates
#' all valid spanning trees via the modified Gabow-Myers algorithm, scores each
#' with \code{\link{calcTreeScores_n}}, and returns the highest-scoring tree
#' together with a visualisation.
#'
#' @param clones data.frame. Clone table (one row per clone) containing at
#'   minimum the columns \code{cluster_num} (clone label), \code{f_pre}
#'   (pre-treatment CCF), and \code{f_day85} (on-treatment CCF).  Typically
#'   the third element (\code{[[3]]}) returned by \code{\link{cluster_data}}.
#' @param lineage_precedence_thresh Numeric. Maximum allowed CCF excess of a
#'   child clone over its parent when building the initial candidate edge set
#'   and filtering edges.  Default \code{0.2}.
#' @param sum_filter_thresh Numeric. Maximum allowed excess of the summed
#'   children CCFs over the parent CCF when pruning spanning trees.
#'   Default \code{0.2}.
#' @param linear_penalty Numeric. Penalty weight applied to linear-chain tree
#'   topologies during scoring.  Default \code{0} (no penalty).
#'
#' @return A list of length 3:
#' \describe{
#'   \item{\code{[[1]]}}{data.frame. Edge list of the best-scoring spanning
#'     tree (columns: \code{parent}, \code{child}).}
#'   \item{\code{[[2]]}}{Numeric matrix. CCF matrix used for scoring (rows =
#'     clones named by \code{cluster_num}, columns = time points
#'     \code{f_pre} and \code{f_day85}).}
#'   \item{\code{[[3]]}}{igraph plot object returned by
#'     \code{\link{plotTree}}.}
#' }
#'
#' @export
build_tree <- function(clones, lineage_precedence_thresh=0.2, sum_filter_thresh=0.2, linear_penalty=0){
    mcf_mat <- as.matrix(clones[, c("f_pre", "f_day85")])
    rownames(mcf_mat) <- clones$cluster_num
    graph_G_pre <- prepareGraph(mcf_mat, lineage_precedence_thresh)
    graph_G <- filterEdgesBasedOnCCFs(graph_G_pre, mcf_mat, thresh = lineage_precedence_thresh)
    graph_G <- prune(graph_G, mcf_mat) %>% arrange(child)

    # Floor the trunk cluster to (1,1) for sum-condition checking and scoring.
    # This prevents silent NULL returns when DP-GMM under-estimates trunk CCF.
    # The original mcf_mat (unmodified) is still returned to callers.
    scoring_mat <- mcf_mat
    scoring_mat[which.max(rowSums(scoring_mat)), ] <- 1

    all_spanning_trees <- enumerateSpanningTreesModified(graph_G, scoring_mat, sum_filter_thresh = sum_filter_thresh)
    if (length(all_spanning_trees) == 0) return(NULL)
    scores <- calcTreeScores_n(scoring_mat, all_spanning_trees, linear_penalty = linear_penalty)
    best_tree <- all_spanning_trees[[which.max(scores)]]
    tree <- plotTree(best_tree)
    return(list(best_tree, mcf_mat, tree))
  }

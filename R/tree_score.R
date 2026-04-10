# -------------------------------------------------------------------------
# TREE SCORING
# -------------------------------------------------------------------------

#' Score a collection of spanning trees against a CCF matrix
#'
#' Computes a fitness score for each spanning tree in \code{trees} using a
#' combined topology cost and mass cost, optionally in parallel.  Higher scores
#' indicate better agreement between the tree structure and the observed CCFs.
#'
#' @param mcf_matrix Numeric matrix. CCF matrix (rows = clones, columns =
#'   samples).
#' @param trees List of data.frames. Each element is a spanning tree edge list
#'   as produced by \code{\link{enumerateSpanningTreesModified}}.
#' @param mc.cores Integer. Number of cores for \code{parallel::mclapply}.
#'   Default \code{1}.
#' @param weight_mass Numeric. Weight applied to the mass cost term.
#'   Default \code{1}.
#' @param weight_topology Numeric. Weight applied to the topology cost term.
#'   Default \code{1}.
#' @param scaling_coeff Numeric. Exponential scaling coefficient; larger values
#'   create a sharper peak around the best tree.  Default \code{5}.
#' @param zero_thresh Numeric. CCF values below this threshold are treated as
#'   zero when computing the CPOV matrix.  Default \code{0.001}.
#' @param restriction_val Numeric. Penalty value placed in the CPOV matrix for
#'   parent-child pairs that violate the CCF ordering constraint.  Default
#'   \code{1}.
#' @param tol Numeric. Tolerance added to the parent CCF before checking for
#'   violations.  Default \code{1e-6}.
#'
#' @return Numeric vector of fitness scores, one per tree in \code{trees}.
#'
#' @export
calcTreeScores_n <- function(mcf_matrix, trees, mc.cores = 1,
                             weight_mass = 1, weight_topology = 1, scaling_coeff = 5,
                             zero_thresh = 0.001, restriction_val = 1, tol = 1e-6) {
  cpov <- create.cpov_n(mcf_matrix, zero.thresh = zero_thresh, restriction.val = restriction_val, tol = tol)
  schism_scores <- unlist(parallel::mclapply(trees,
                                             function(x) calcTreeFitness(x, cpov, mcf_matrix, am_format = "edges",
                                                                         weight_mass = weight_mass,
                                                                         weight_topology = weight_topology,
                                                                         scaling_coeff = scaling_coeff),
                                             mc.cores = mc.cores))
  return(schism_scores)
}
#' Compute the fitness score for a single spanning tree
#'
#' Calculates the combined topology cost and mass cost for one tree and
#' converts the total penalty to a fitness score via \code{exp(-scaling_coeff * Z)}.
#'
#' @param admat data.frame or matrix. Adjacency representation of the tree.
#'   Use \code{am_format = "edges"} when passing the edge-list format returned
#'   by \code{\link{enumerateSpanningTreesModified}}.
#' @param cpov Matrix. Conditional probability-of-violation matrix as returned
#'   by \code{\link{create.cpov_n}}.
#' @param mcf_matrix Numeric matrix. CCF matrix (rows = clones, columns =
#'   samples).
#' @param am_format Character. Format of \code{admat}: \code{"long"} (default)
#'   or \code{"edges"}.
#' @param weight_mass Numeric. Weight for the mass cost.  Default \code{1}.
#' @param weight_topology Numeric. Weight for the topology cost.  Default
#'   \code{1}.
#' @param scaling_coeff Numeric. Exponential scaling coefficient.
#'   Default \code{5}.
#'
#' @return Scalar numeric fitness score in \code{(0, 1]}.
#'
#' @export
calcTreeFitness <- function(admat, cpov, mcf_matrix, am_format = "long", weight_mass = 1, weight_topology = 1, scaling_coeff=5) {
  if (am_format == "edges") {
    admat <- edgesToAmLong(admat)
    am_format <- "long"
  }
  
  TC <- calcTopologyCost(admat, cpov, am_format)
  MC <- calcMassCost(admat, mcf_matrix, am_format)
  Z <- weight_topology * TC + weight_mass * MC
  fitness <- exp(-scaling_coeff * Z)
  fitness
}
#' Build the conditional probability-of-violation (CPOV) matrix
#'
#' Constructs a \eqn{(K+1) \times K} matrix encoding whether each directed
#' parent-child clone pair violates the CCF ordering constraint.  The root node
#' occupies row 1; clone \eqn{i} occupies row \eqn{i+1}.  Self-edges
#' (diagonal) are set to \code{restriction_val} to prevent self-loops.
#'
#' @param mcf_matrix Numeric matrix or data.frame. CCF matrix (rows = clones,
#'   columns = samples).
#' @param zero.thresh Numeric. CCF values below this threshold are treated as
#'   zero.  Default \code{0.001}.
#' @param restriction.val Numeric. Penalty value for violating parent-child
#'   pairs.  Default \code{1}.
#' @param tol Numeric. Tolerance added to parent CCF when checking violations.
#'   Default \code{1e-6}.
#'
#' @return Integer matrix of dimension \eqn{(K+1) \times K} where rows are
#'   labelled \code{"root"}, \code{"1"}, ..., \code{"K"} and columns are
#'   \code{"1"}, ..., \code{"K"}.  Entries are \code{0} (no violation) or
#'   \code{restriction_val} (violation).
#'
#' @export
create.cpov_n <- function(mcf_matrix, zero.thresh = 0.001, restriction.val = 1, tol = 1e-6) {
  if (is.null(mcf_matrix) || !is.matrix(mcf_matrix) && !is.data.frame(mcf_matrix)) {
    stop("mcf_matrix must be a numeric matrix/data.frame with rows=clusters and cols=samples.")
  }
  MCF <- as.matrix(mcf_matrix)
  K <- nrow(MCF); S <- ncol(MCF)
  
  rnames <- c("root", rownames(MCF))
  cnames <- if (!is.null(rownames(MCF))) rownames(MCF) else paste0("C", seq_len(K))
  cpov <- matrix(0L, nrow = K + 1, ncol = K, dimnames = list(rnames, cnames))
  
  cpov[1, ] <- 0L
  diag_idx <- cbind(2:(K + 1), 1:K) 
  cpov[diag_idx] <- restriction.val
  
  for (r in 2:(K + 1)) {
    from <- r - 1
    parent <- MCF[from, ]
    for (c in 1:K) {
      if (cpov[r, c] == restriction.val) next 
      child <- MCF[c, ]
      idx <- child > zero.thresh
      if (!any(idx)) {
        cpov[r, c] <- 0L
      } else {
        violates <- any(child[idx] > parent[idx] + tol)
        cpov[r, c] <- ifelse(violates, restriction.val, 0L)
      }
    }
  }
  cpov
}
#' Calculate the topology cost for a spanning tree
#'
#' Sums the CPOV penalties for all edges present in the adjacency matrix.  A
#' topology cost of zero means no parent-child CCF ordering is violated.
#'
#' @param am data.frame or matrix. Adjacency representation of the tree (long
#'   or wide format).
#' @param cpov Matrix. CPOV matrix from \code{\link{create.cpov_n}}.
#' @param am_format Character. \code{"long"} (default) or \code{"wide"}.
#'
#' @return Scalar numeric topology cost (\eqn{\geq 0}).
#'
#' @export
calcTopologyCost <- function(am, cpov, am_format = "long") {
  TC <- 0
  if (am_format == "long") {
    am <- toWide(am)
  } 
  edges <- which(am == 1, arr.ind=TRUE)
  for (i in seq_len(nrow(edges))) {
    TC <- TC + cpov[edges[i,1], edges[i,2]]
  }
  TC
}
#' Calculate the mass cost for a spanning tree
#'
#' Penalises trees where the sum of children CCFs exceeds the parent CCF.  For
#' the long format the maximum per-sample excess is summed across all parent
#' nodes; for the wide format the Euclidean norm of per-sample excesses is
#' summed.
#'
#' @param am data.frame or matrix. Adjacency representation of the tree.
#' @param mcf_matrix Numeric matrix. CCF matrix (rows = clones, columns =
#'   samples).
#' @param am_format Character. \code{"long"} (default) or \code{"wide"}.
#'
#' @return Scalar numeric mass cost (\eqn{\geq 0}).
#'
#' @export
calcMassCost <- function(am, mcf_matrix, am_format="long") {
  if (am_format == "long") {
    edges <- getEdges(am)
    parent_nodes <- unique(edges$parent)
    mass_cost <- rep(0, length(parent_nodes))
    
    for (i in seq_len(length(parent_nodes))) {
      parent_node <- parent_nodes[i]
      if (parent_node == "root") {
        parent_w <- rep(1, ncol(mcf_matrix))
      } else {
        parent_w <- mcf_matrix[as.numeric(parent_node), ,drop=FALSE]
      }
      
      kids <- getChildren(am, parent_node)
      if (length(kids) > 1) {
        children_w <- colSums(mcf_matrix[as.numeric(kids), ,drop=FALSE])
      } else {
        children_w <- mcf_matrix[as.numeric(kids), ,drop=FALSE]
      }
      
      mc_s <- ifelse(parent_w >= children_w, 0, children_w - parent_w)
      mass_cost[i] <- max(mc_s) 
    }
    return(sum(mass_cost))
    
  } else if (am_format == "wide") {
    num_children <- rowSums(am, na.rm = T)
    nodes <- which(num_children > 0, arr.ind = T) 
    mc_node <- rep(0, length(nodes))
    
    for (i in 1:length(nodes)) {
      node <- nodes[i]
      
      parent_w <- rep(1, ncol(mcf_matrix))
      if (node != 1) { 
        parent_w <- mcf_matrix[node-1,]
      }
      
      kids <- which(am[node,] == 1, arr.ind = T)
      if (num_children[node] > 1) {
        children_w <- colSums(mcf_matrix[kids, ])
      } else {
        children_w <- mcf_matrix[kids, ]
      }
      
      mc_s <- ifelse(parent_w >= children_w, 0, children_w - parent_w)
      mc_node[i] <- sqrt(sum(mc_s^2))
    }
    return(sum(mc_node))
  }
}
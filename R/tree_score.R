library(parallel)

# -------------------------------------------------------------------------
# TREE SCORING
# -------------------------------------------------------------------------
#' @export
calcTreeScores_n <- function(mcf_matrix, trees, mc.cores = 1) {
  cpov <- create.cpov_n(mcf_matrix)
  schism_scores <- unlist(parallel::mclapply(trees, 
                                             function(x) calcTreeFitness(x, cpov, mcf_matrix, am_format = "edges"),
                                             mc.cores = mc.cores))
  return(schism_scores)
}
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
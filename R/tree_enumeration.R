# -------------------------------------------------------------------------
# GRAPH PREPARATION
# -------------------------------------------------------------------------

#' Create the initial candidate edge set for the phylogenetic tree graph
#'
#' Constructs a data.frame of directed edges (parent → child) representing all
#' clone pairs where the child's CCF does not exceed the parent's CCF by more
#' than \code{thresh} in any sample.  Every clone is also given an edge from a
#' virtual \code{"root"} node.
#'
#' @param mcf_mat Numeric matrix. CCF matrix with rows = clones and
#'   columns = time-point samples.  Row names are used as node labels.
#' @param thresh Numeric. Maximum allowed CCF excess of child over parent for
#'   an edge to be included.
#'
#' @return data.frame with columns \code{edge} (e.g. \code{"root->1"}),
#'   \code{parent}, and \code{child}.
#'
#' @export
prepareGraph <- function(mcf_mat, thresh) {
  graph_pre <- data.frame(edge = character(), parent = character(), child = character())
  for (i in seq_len(nrow(mcf_mat))) {
    graph_pre <- graph_pre %>% dplyr::add_row(edge = paste("root->", i, sep = ""), parent = "root", child = as.character(i))
    for (j in seq_len(nrow(mcf_mat))) {
      if (i!=j) {
        if (all(mcf_mat[j, ] - mcf_mat[i, ] >= -thresh)) {
          graph_pre <- graph_pre %>% dplyr::add_row(edge = paste(j, "->", i, sep = ""), parent = as.character(j), child = as.character(i))
        }
      }
    }
  }
  return(graph_pre)
}

#' Filter candidate edges based on lineage precedence constraints
#'
#' Removes edges from a candidate graph where the child clone's CCF exceeds the
#' parent clone's CCF by more than \code{thresh} in any sample, enforcing the
#' biological constraint that a child lineage cannot be more prevalent than its
#' ancestor.
#'
#' @param graph_G data.frame. Candidate edge table as produced by
#'   \code{\link{prepareGraph}} (columns: \code{edge}, \code{parent},
#'   \code{child}).
#' @param mcf Numeric matrix. CCF matrix (rows = clones, columns = samples).
#' @param thresh Numeric. Maximum allowed CCF excess of child over parent.
#'   Default \code{0.1}.
#'
#' @return Filtered data.frame with the same columns as \code{graph_G}.
#'
#' @export
filterEdgesBasedOnCCFs <- function(graph_G, mcf, thresh = 0.1) {
  check_edges_logical <- apply(graph_G, 1, function(edge) checkEdge(edge, mcf, thresh))
  filtered_graph_G <- graph_G[check_edges_logical, ]
  return(filtered_graph_G)
}

#' Prune biologically invalid edges from the phylogenetic graph
#'
#' Removes edges that violate birth-order or continuity constraints: a child
#' clone cannot appear before its parent (birth violation), and neither clone
#' should be absent for an interval and then reappear (continuity violation).
#'
#' @param graph data.frame. Edge table (columns: \code{edge}, \code{parent},
#'   \code{child}).
#' @param mcf_mat Numeric matrix. CCF matrix (rows = clones, columns = samples).
#'
#' @return Pruned data.frame with the same columns as \code{graph}, sorted by
#'   \code{child}.
#'
#' @export
prune <- function(graph, mcf_mat){
  roots <- graph %>% dplyr::filter(parent=='root')
  
  tmp = graph %>%
    dplyr::filter(parent!='root') %>%
    dplyr::rowwise() %>%
    dplyr::mutate(p_birth=min(which(mcf_mat[as.numeric(parent),]!=0)),
                  p_last_seen=max(which(mcf_mat[as.numeric(parent),]!=0)),
                  p_die=max(which(mcf_mat[as.numeric(parent),]==0)),
                  p_die=ifelse(p_die==-Inf, ncol(mcf_mat), p_die),
                  c_birth=min(which(mcf_mat[as.numeric(child),]!=0)),
                  c_last_seen=max(which(mcf_mat[as.numeric(child),]!=0)),
                  c_die=max(which(mcf_mat[as.numeric(child),]==0)),
                  c_die=ifelse(c_die==-Inf, ncol(mcf_mat), c_die),
                  vio_birth=p_birth > c_birth,
                  vio_continue=(p_die < p_last_seen) | (c_die < c_last_seen)) %>%
    dplyr::filter(!(vio_birth)) %>%
    dplyr::select(edge, parent, child) %>%
    rbind(roots) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(child)
  
  return(tmp)
}

# -------------------------------------------------------------------------
# GABOW-MYERS ALGORITHM
# -------------------------------------------------------------------------

#' Enumerate all valid spanning trees using the modified Gabow-Myers algorithm
#'
#' Recursively enumerates all spanning trees of the candidate graph that
#' satisfy the CCF sum condition (the sum of children CCFs does not exceed the
#' parent CCF by more than \code{sum_filter_thresh}).  Results are written to
#' the global variable \code{all_spanning_trees}.
#'
#' @param graph_G data.frame. Pruned edge table (columns: \code{edge},
#'   \code{parent}, \code{child}), typically produced by \code{\link{prune}}.
#' @param mcf Numeric matrix. CCF matrix (rows = clones, columns = samples).
#' @param sum_filter_thresh Numeric. Maximum allowed excess of summed children
#'   CCFs over the parent CCF.  Default \code{0.2}.
#'
#' @return Invisibly \code{NULL}.  All valid spanning trees are stored in the
#'   global variable \code{all_spanning_trees} as a list of edge-list
#'   data.frames.
#'
#' @export
enumerateSpanningTreesModified <- function(graph_G, mcf, sum_filter_thresh=0.2) {
  # Initialize globals required for recursion
  all_spanning_trees <- assign("all_spanning_trees", list(), envir = .GlobalEnv)
  F_tb <- assign("F_tb", dplyr::filter(graph_G, parent == "root"), envir = .GlobalEnv)
  
  all_vertices <- verticesInGraph(graph_G)
  tree_T <- dplyr::tibble(parent = character(), child = character())
  
  growModified(tree_T, all_vertices, mcf, sum_filter_thresh)
}

#' Recursive worker for the modified Gabow-Myers tree enumeration
#'
#' Internal recursive function called by \code{\link{enumerateSpanningTreesModified}}.
#' Grows a partial spanning tree one edge at a time, checking the CCF sum
#' condition at each step, and stores each complete spanning tree in the global
#' \code{all_spanning_trees} list.
#'
#' @param tree_T tibble. Current partial spanning tree (columns: \code{parent},
#'   \code{child}).
#' @param all_vertices Character vector. All vertex labels in the graph.
#' @param w Numeric matrix. CCF matrix (rows = clones, columns = samples).
#' @param sum_thresh Numeric. CCF sum-condition threshold.  Default \code{0.2}.
#'
#' @return Invisibly \code{NULL}; results are accumulated in the global
#'   \code{all_spanning_trees}.
#'
#' @export
growModified <- function(tree_T, all_vertices, w, sum_thresh=0.2) {
  
  if (length(verticesInGraph(tree_T)) == length(all_vertices) & nrow(tree_T) == (length(all_vertices)-1)) {
    assign("all_spanning_trees", c(all_spanning_trees, list(tree_T)), envir = .GlobalEnv)
    
  } else {
    FF <- dplyr::tibble(parent = character(), child = character())
    bridge <- FALSE
    
    while(!bridge) {
      if (nrow(F_tb) == 0) stop("F_tb is empty")
      edge_e <- pop(F_tb, "F_tb")
      v <- edge_e$child
      tree_T <- rbind(tree_T, edge_e)
      
      # Check sum condition constraint
      if (satisfiesSumCondition(tree_T, w, sum_thresh)) {
        # update F: push each edge (v,w), w not in T onto F
        in_T <- verticesInGraph(tree_T)
        temp_add_to_F <- dplyr::filter(graph_G, parent == v, !(child %in% in_T))
        assign("F_tb", rbind(temp_add_to_F, F_tb), envir = .GlobalEnv)
        
        # remove each edge (w,v), w in T from F
        w_in_T <- verticesInGraph(tree_T)
        removed_edges <- dplyr::filter(F_tb, parent %in% w_in_T, child == v)
        assign("F_tb", dplyr::filter(F_tb, !edge %in% removed_edges$edge), envir = .GlobalEnv)
        
        # Recurse
        growModified(tree_T, all_vertices, w, sum_thresh)
        
        # Restore F
        not_in_T <- all_vertices[!all_vertices %in% verticesInGraph(tree_T)]
        if (length(not_in_T) > 0 & nrow(F_tb) > 0) {
          edges_to_remove_9 <- paste0(v, "->", not_in_T)
          assign("F_tb", dplyr::filter(F_tb, !edge %in% edges_to_remove_9), envir = .GlobalEnv)
        }
        assign("F_tb", rbind(removed_edges, F_tb), envir = .GlobalEnv)
      }
      
      # delete e from T and from G, add e to FF
      tree_T <- tree_T[tree_T$edge != edge_e$edge, ]
      assign("graph_G", graph_G[graph_G$edge != edge_e$edge, ], envir = .GlobalEnv)
      FF <- rbind(edge_e, FF)
      
      # bridge test
      bridge <- bridgeTestBFS(graph_G, edge_e)
    }
    
    # Restore edges
    if (nrow(FF) > 0) {
      assign("graph_G", rbind(FF, graph_G), envir = .GlobalEnv)
      while (nrow(FF) > 0) {
        assign("F_tb", rbind(FF[1, ], F_tb), envir = .GlobalEnv)
        FF <- FF[-1, ]
      }
    }
  }
}
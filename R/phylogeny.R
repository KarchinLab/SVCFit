### Gabow Mayor to build phylogeny

#' Create tibble of possible edges from CCF values based on w_mat only
#' 
#' @export
#' @param w matrix of CCF values (rows = clusters, columns = samples)
#' @return graph_G tibble of possible edges with columns edge, parent, child
prepareGraph <- function(mcf_mat, thresh) {
  graph_pre <- data.frame(edge = character(), parent = character(), child = character())
  for (i in seq_len(nrow(mcf_mat))) {
    graph_pre <- graph_pre %>% add_row(edge = paste("root->", i, sep = ""), parent = "root", child = as.character(i))
    for (j in seq_len(nrow(mcf_mat))) {
      if (i!=j) {
        i_row = mcf_mat[i, ]
        j_row = mcf_mat[j, ]
        if (all(j_row-i_row >= -thresh)) {
          graph_pre <- graph_pre %>% add_row(edge = paste(j, "->", i, sep = ""), parent = as.character(j), child = as.character(i))
        }
      }
    }
  }
  return(graph_pre)
}


#' Filter possible edges based on lineage precedence 
#' 
#' @export
#' @param graph_G tibble of possible edges with columns edge, parent, child
#' @param w matrix of CCF values (rows = clusters, columns = samples)
#' @param thresh maximum allowed violation of lineage precedence (default = 0.1)
filterEdgesBasedOnCCFs <- function(graph_G, mcf, thresh = 0.1) {
  check_edges_logical <- apply(graph_G, 1, function(edge) checkEdge(edge, mcf, thresh))
  filtered_graph_G <- graph_G[check_edges_logical, ]
  return(filtered_graph_G)
}

checkEdge <- function(edge, mcf, thresh = 0.2) {
  # returns TRUE if satisfies lineage precedence with given threshold
  # returns FALSE if violates i.e. child_ccf - parent_ccf > thresh in any sample
  # edge is in the format c(edge_name, parent, child)
  
  # in case of factors
  p <- as.character(edge[2])
  c <- as.character(edge[3])
  
  if (p == "root") {
    parent_ccfs <- rep(1, ncol(mcf))
  } else {
    parent_ccfs <- mcf[as.numeric(p), ]
  }
  child_ccfs <- mcf[as.numeric(c), ]
  
  diff <- child_ccfs - parent_ccfs
  if (any(diff > thresh)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


prune <- function(graph, mcf_mat){
  
  roots <- graph %>% filter(parent=='root')
  
  tmp = graph %>%
    filter(parent!='root')%>%
    rowwise()%>%
    mutate(p_birth=min(which(mcf_mat[as.numeric(parent),]!=0)),
           p_last_seen=max(which(mcf_mat[as.numeric(parent),]!=0)),
           p_die=max(which(mcf_mat[as.numeric(parent),]==0)),
           p_die=ifelse(p_die==-Inf, ncol(mcf_mat), p_die),
           c_birth=min(which(mcf_mat[as.numeric(child),]!=0)),
           c_last_seen=max(which(mcf_mat[as.numeric(child),]!=0)),
           c_die=max(which(mcf_mat[as.numeric(child),]==0)),
           c_die=ifelse(c_die==-Inf, ncol(mcf_mat), c_die),
           vio_birth=p_birth > c_birth,
           vio_continue=(p_die < p_last_seen) | (c_die < c_last_seen))%>%
    filter(!(vio_birth))%>%
    select(edge, parent, child)%>%
    rbind(roots)%>%
    ungroup()%>%
    arrange(child)
  
  return(tmp)
}


#' Enumerate all spanning trees using modified Gabow-Myers
#' 
#' @export
#' @param graph_G tibble of possible edges with columns edge, parent, child
#' @param w matrix of CCF values (rows = clusters, columns = samples)
#' @param sum_filter_thresh thresh maximum allowed violation of Sum Condition (default = 0.2)
enumerateSpanningTreesModified <- function(graph_G, mcf, purity, sum_filter_thresh=0.2) {
  # all_spanning_trees must be set as an empty list, global variable, before function is called
  # graph_G must be set as global variable before function is called
  all_spanning_trees <- assign("all_spanning_trees", list(), envir = .GlobalEnv)
  #filtered_trees <- assign("filtered_trees", list(), envir = .GlobalEnv)
  F_tb <- assign("F_tb", filter(graph_G, parent == "root"), envir = .GlobalEnv)
  all_vertices <- verticesInGraph(graph_G)
  tree_T <- tibble(parent = character(), child = character())
  
  growModified(tree_T, all_vertices, mcf, purity, sum_filter_thresh)
}


verticesInGraph <- function(tb) {
  unique(c(tb$parent, tb$child))
}


growModified <- function(tree_T, all_vertices, w, sum_thresh=0.2) {
  
  if (length(verticesInGraph(tree_T)) == length(all_vertices) & nrow(tree_T) == (length(all_vertices)-1)) {
    assign("all_spanning_trees", c(all_spanning_trees, list(tree_T)), envir = .GlobalEnv)
    
  } else {
    FF <- tibble(parent = character(), child = character())
    
    bridge <- FALSE
    while(!bridge) {
      # new tree edge
      if (nrow(F_tb) == 0) stop("F_tb is empty")
      edge_e <- pop(F_tb, "F_tb")
      v <- edge_e$child
      tree_T <- rbind(tree_T, edge_e)
      
      # check if adding this node does not violate the constraint
      if (satisfiesSumCondition(tree_T, w, sum_thresh)) {
        # update F
        ## push each edge (v,w), w not in T onto F
        in_T <- verticesInGraph(tree_T)
        temp_add_to_F <- filter(graph_G, parent == v, !(child %in% in_T))
        # temp_add_to_F
        assign("F_tb", rbind(temp_add_to_F, F_tb), envir = .GlobalEnv)
        
        ## remove each edge (w,v), w in T from F
        w_in_T <- verticesInGraph(tree_T)
        removed_edges <- filter(F_tb, parent %in% w_in_T, child == v)
        assign("F_tb", filter(F_tb, !edge %in% removed_edges$edge), envir = .GlobalEnv)
        
        # recurse
        growModified(tree_T, all_vertices, w, sum_thresh)
        
        # restore F
        # pop each edge (v,w), w not in T, from F
        not_in_T <- all_vertices[!all_vertices %in% verticesInGraph(tree_T)]
        if (length(not_in_T) > 0 & nrow(F_tb) > 0) {
          edges_to_remove_9 <- paste0(v, "->", not_in_T)
          assign("F_tb", filter(F_tb, !edge %in% edges_to_remove_9), envir = .GlobalEnv)
        }
        # restore each edge (w,v), w in T, in F
        assign("F_tb", rbind(removed_edges, F_tb), envir = .GlobalEnv)
        
      }
      # delete e from T and from G, add e to FF
      tree_T <- tree_T[tree_T$edge != edge_e$edge, ]
      assign("graph_G",graph_G[graph_G$edge != edge_e$edge, ], envir = .GlobalEnv)
      FF <- rbind(edge_e, FF)
      
      # bridge test
      bridge <- bridgeTestBFS(graph_G, edge_e)
    }
    
    # pop each edge e from FF, push e onto F,and add e to G
    if (nrow(FF) > 0) {
      
      # pop and push all edges at once (same order)
      # assign("F_tb", rbind(FF, F_tb), envir = .GlobalEnv)
      assign("graph_G", rbind(FF, graph_G), envir = .GlobalEnv)
      
      # pop and push edges one by one (rev order in F)
      while (nrow(FF) > 0) {
        assign("F_tb", rbind(FF[1, ], F_tb), envir = .GlobalEnv)
        FF <- FF[-1, ]
      }
    }
  }
}

pop <- function(edges_tb, tb_name) {
  assign(tb_name, edges_tb[-1, ], envir = .GlobalEnv)
  return(edges_tb[1, ])
}

bfsLong2 <- function(graph_G) {
  # returns vector of nodes in main tree (connected to root) including "root" 
  # starting at root
  # does not stop if there is a cycle in graph 
  graph_G$parent <- as.character(graph_G$parent)
  children <- graph_G[(graph_G$parent == "root"), ]$child
  nodes <- c("root", children)
  
  while(length(children) > 0) {
    c <- children[1]
    temp.children <- graph_G[(graph_G$parent == c), ]$child
    
    # remove children already seen
    if (any(temp.children %in% nodes)) {
      temp.children <- temp.children[! temp.children %in% nodes]
    }
    
    children <- c(children, temp.children)
    
    nodes <- c(nodes, temp.children)
    children <- children[-1]
  }
  return(nodes)
}

toWide <- function(am.long){
  am.long$child <- as.numeric(am.long$child)
  am.long %>% select(parent, child, connected) %>%
    tidyr::spread(child, connected) %>%
    select(-parent) %>%
    as.matrix()
}

getEdges <- function(am.long) {
  am.long %>%
    filter(connected == 1) %>%
    mutate(parent = as.character(parent))
}

satisfiesSumCondition <- function(edges, w, thresh = 0.2) {
  # returns TRUE if sum condition is not violated with given threshold (default 0.2)
  
  edges$parent <- as.character(edges$parent)
  all_parents <- unique(edges$parent)
  
  for (p in all_parents) {
    # get parent CCF
    if (p == "root") {
      parent_ccf <- rep(1, ncol(w))
    } else {
      parent_ccf <- w[as.numeric(p), ]
    }
    
    # get children CCF (sum if more than 1 child)
    children <- as.numeric(filter(edges, parent == p)$child)
    if (length(children) > 1) {
      children_ccf <- colSums(w[children, ,drop=FALSE])
    } else {
      children_ccf <- w[children, ]
    }
    
    diff <- children_ccf - parent_ccf
    if (any(diff > thresh)) return(FALSE)
  }
  
  # sum condition is never violated, return TRUE
  return(TRUE)
}


bridgeTestBFS <- function(graph_G, edge_e) {
  node_to_check <- edge_e$child
  
  nodes_connected_to_root <- bfsLong2(graph_G)
  !(node_to_check %in% nodes_connected_to_root)
}


create.cpov_n <- function(mcf_matrix, zero.thresh = 0.001, restriction.val = 1, tol = 1e-6) {
  if (is.null(mcf_matrix) || !is.matrix(mcf_matrix) && !is.data.frame(mcf_matrix)) {
    stop("mcf_matrix must be a numeric matrix/data.frame with rows=clusters and cols=samples.")
  }
  MCF <- as.matrix(mcf_matrix)
  K <- nrow(MCF); S <- ncol(MCF)
  
  # Initialize CPOV: (K+1) x K. Row 1 is the root.
  rnames <- c("root", rownames(MCF))
  cnames <- if (!is.null(rownames(MCF))) rownames(MCF) else paste0("C", seq_len(k))
  cpov <- matrix(0L, nrow = K + 1, ncol = K, dimnames = list(rnames, cnames))
  
  # Root row: all zeros (root can go to anyone)
  cpov[1, ] <- 0L
  
  # Forbid self-edges (cluster -> itself)
  diag_idx <- cbind(2:(K + 1), 1:K) # row r maps to from = r-1
  cpov[diag_idx] <- restriction.val
  
  # Fill remaining entries by deterministic presence/inequality rule
  for (r in 2:(K + 1)) {
    from <- r - 1
    parent <- MCF[from, ]
    for (c in 1:K) {
      if (cpov[r, c] == restriction.val) next  # skip self-edge already set
      child <- MCF[c, ]
      # consider samples where the child is present
      idx <- child > zero.thresh
      if (!any(idx)) {
        # If the child is never present, we don't restrict the edge
        cpov[r, c] <- 0L
      } else {
        # If child exceeds parent (by more than tol) in any present sample -> forbid
        violates <- any(child[idx] > parent[idx] + tol)
        cpov[r, c] <- ifelse(violates, restriction.val, 0L)
      }
    }
  }
  cpov
}

calcTreeScores_n <- function(mcf_matrix, trees, purity, mc.cores = 1) {
  # calculate mean and sd of mcf for each cluster in each sample
  #mcf_stats <- summarizeWChain(mcf_chain) 
  
  # create cpov matrix
  # first: using sample presence to create a binary matrix; 0 is i can be a ancestor of j; 1 if not
  # second: for i,j pair that pass the sample presence test, calc stats of the difference between mcf of all 
  # samples; return binary matrix
  # problem: is the stats calculation in create.cpov correct? only using cluster so is actually POV instead 
  # of CPOV
  cpov <- create.cpov_n(mcf_matrix)
  #mcf_mat <- estimateMCFs(mcf_chain)
  
  # first calculate topology cost: sum of the cpov matrix over all edges
  # potential problem: mcf[i]=0.4->mcf[j]=0.6 has same weight as mcf[i]=0.5->mcf[j]=0.6
  # second calculate mass cost: take the max mass violation among all samples; is there a better strategy?
  # fitness is exp(-5*(topology cost + mass cost))
  schism_scores <- unlist(parallel:::mclapply(trees, 
                                              function(x) calcTreeFitness(x, cpov, mcf_mat, purity, am_format = "edges"),
                                              mc.cores = mc.cores))
  return(schism_scores)
}

calcTreeFitness <- function(admat, cpov, mcf_matrix, purity, am_format = "long", weight_mass = 1, weight_topology = 1, scaling_coeff=5) {
  # if only edges are given, change into long format
  if (am_format == "edges") {
    admat <- edgesToAmLong(admat)
    am_format <- "long"
  }
  
  TC <- calcTopologyCost(admat, cpov, am_format)
  MC <- calcMassCost(admat, mcf_matrix, purity, am_format)
  Z <- weight_topology * TC + weight_mass * MC
  fitness <- exp(-scaling_coeff * Z)
  fitness
}

edgesToAmLong <- function(edges) {
  am_wide <- initEmptyAdmatFromK(length(unique(edges$child)))
  edges[edges$parent == "root", "parent"] <- "0"
  edges <- edges %>%
    mutate(parent = as.numeric(parent) + 1,
           child = as.numeric(child)) %>%
    select(parent, child)
  edges <- as.matrix(edges)
  for (r in 1:nrow(edges)) {
    am_wide[edges[r,1], edges[r,2]] <- 1
  }
  admat <- toLong(am_wide)
  admat <- reversedEdges(admat) %>%
    mutate(reversed_connected=reversedConnection(.),
           bi_directional=NA,
           root_connected=NA)    
  admat <- updateGraphElements(admat)
  return(admat)
}

initEmptyAdmatFromK <- function(K) {
  admat <- matrix(0, K, K)
  diag(admat) <- NA
  am2 <- rbind(0, admat)
  dimnames(am2) <- list(c("root", 1:K), 1:K)
  return(am2)
}

toLong <- function(am) {
  am.long <- as_tibble(am) %>%
    mutate(parent=rownames(am)) %>%
    pivot_longer(-parent,
                 names_to="child",
                 values_to="connected") %>%
    filter(parent != child) %>%
    unite("edge", c("parent", "child"), sep="->",
          remove=FALSE) %>%
    mutate(parent=factor(parent, levels=unique(parent)))   
  return(am.long)
}

reversedEdges <- function(am) {
  am2 <- am %>%
    filter(!is.na(connected)) %>%
    unite("reverse_edge", c("child", "parent"), sep="->",
          remove=FALSE)
  am2
}

reversedConnection <- function(am) {
  connections <- setNames(am$connected, am$edge)
  reversed_connections <- connections[am$reverse_edge] %>%
    "["(!is.na(.))
  reversed <- setNames(rep(0, nrow(am)), am$reverse_edge)
  reversed[names(reversed_connections)] <- reversed_connections
  reversed
}

isBidirectional <- function(am) {
  am %>%
    mutate(bi_directional=(reverse_edge %in% edge) &
             connected==1 &
             reversed_connected == 1) %>%
    pull(bi_directional)
}

isRootConnected <- function(am) isParentConnected(am)[1]

isParentConnected <- function(am) {
  am %>%
    mutate(parent=factor(parent, levels=unique(parent))) %>%
    group_by(parent) %>%
    summarize(n=sum(connected)) %>%
    pull(n) > 0
}

updateGraphElements <- function(am) {
  am %>%
    mutate(parent=factor(parent, levels=unique(parent))) %>%
    mutate(reversed_connected=reversedConnection(.)) %>%
    mutate(bi_directional=isBidirectional(.)) %>%
    mutate(root_connected=isRootConnected(.))
}

calcTopologyCost <- function(am, cpov, am_format = "long") {
  TC <- 0
  
  if (am_format == "long") {
    am <- toWide(am)
  } 
  
  edges <- which(am == 1, arr.ind=TRUE)
  N <- nrow(edges)
  for (i in seq_len(N)) {
    TC <- TC + cpov[edges[i,1], edges[i,2]]
  }
  
  TC
}

calcMassCost <- function(am, mcf_matrix, purity, am_format="long") {
  num_samples <- ncol(mcf_matrix)
  
  if (am_format == "long") {
    edges <- getEdges(am)
    
    parent_nodes <- unique(edges$parent)
    mass_cost <- rep(0, length(parent_nodes)) # mass cost of each parent node
    
    for (i in seq_len(length(parent_nodes))) {
      parent_node <- parent_nodes[i]
      
      # root MCF is purity instead of 1
      if (parent_node == "root") {
        # parent_w <- rep(1, num_samples) # 1 replaced by purity
        parent_w <- purity
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
      #mass_cost[i] <- sqrt(sum(mc_s^2))
      mass_cost[i] <- max(mc_s) # take max across samples instead of euclidean distance
    }
    return(sum(mass_cost))
    
  } else if (am_format == "wide") {
    num_children <- rowSums(am, na.rm = T)
    nodes <- which(num_children > 0, arr.ind = T) # not leaves
    mc_node <- rep(0, length(nodes))
    
    for (i in 1:length(nodes)) {
      node <- nodes[i]
      
      # root node: MCF = 1
      parent_w <- rep(1, ncol(mcf_matrix))
      # not root node: look up MCF in mcf_matrix
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


getChildren <- function(am.long, node) {
  # returns vector of children nodes
  edges <- am.long %>%
    mutate(parent = as.character(parent)) %>%
    filter(connected == 1) %>%
    filter(parent == node)
  return(edges$child)
}


colorScheme <- function(edges, palette=viridis::viridis) {
  v_sorted = sort(unique(c(edges$parent, edges$child)))
  v_sorted = c(sort(as.integer(v_sorted[!v_sorted=='G'])), "G")
  # root_idx <- which(v_sorted=="root")
  colors <- c(palette(length(v_sorted)-1), "white")
  v_color <- tibble(v_sorted, colors)
  return(v_color)
}

plotGraph <- function(am.long, v_color){
  # make sure am.long is sorted by parent and child
  am.long <- mutate(am.long, child = as.numeric(am.long$child)) %>%
    arrange(parent, child)
  am.long <- mutate(am.long, child = as.character(am.long$child))
  
  # change to wide format and plot
  am <- toWide(am.long)
  rownames(am) <- c("G", colnames(am))
  am <- cbind(G=0, am) ## add column for root
  colnames(am) <- rownames(am)
  
  am[is.na(am)] <- 0
  
  ig <- igraph::graph_from_adjacency_matrix(am, mode = "directed", weighted = TRUE,
                                            diag = FALSE, add.row = TRUE) 
  V(ig)$color <- as.list(v_color %>% arrange(match(v_sorted, names(V(ig)))) %>% select(colors))$colors
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
                      vertex.size=24, vertex.frame.color = "#000000", vertex.label.cex = 1.5,
                      vertex.label.family = "Helvetica", vertex.label.color = "#000000",
                      edge.arrow.size = 0.5, edge.arrow.width = 2)
}

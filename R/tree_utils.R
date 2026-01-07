library(dplyr)
library(tidyr)

# -------------------------------------------------------------------------
# DATA STRUCTURE HELPERS
# -------------------------------------------------------------------------

#' Pop the first row from a tibble and update the global variable
pop <- function(edges_tb, tb_name) {
  assign(tb_name, edges_tb[-1, ], envir = .GlobalEnv)
  return(edges_tb[1, ])
}

#' Convert long format adjacency tibble to wide matrix
toWide <- function(am.long){
  am.long$child <- as.numeric(am.long$child)
  am.long %>% 
    dplyr::select(parent, child, connected) %>%
    tidyr::spread(child, connected) %>%
    dplyr::select(-parent) %>%
    as.matrix()
}

#' Convert wide matrix to long format tibble
toLong <- function(am) {
  am.long <- dplyr::as_tibble(am) %>%
    dplyr::mutate(parent=rownames(am)) %>%
    tidyr::pivot_longer(-parent, names_to="child", values_to="connected") %>%
    dplyr::filter(parent != child) %>%
    tidyr::unite("edge", c("parent", "child"), sep="->", remove=FALSE) %>%
    dplyr::mutate(parent=factor(parent, levels=unique(parent)))   
  return(am.long)
}

#' Convert edge list to long adjacency matrix with annotations
edgesToAmLong <- function(edges) {
  am_wide <- initEmptyAdmatFromK(length(unique(edges$child)))
  edges[edges$parent == "root", "parent"] <- "0"
  edges <- edges %>%
    dplyr::mutate(parent = as.numeric(parent) + 1,
                  child = as.numeric(child)) %>%
    dplyr::select(parent, child)
  edges <- as.matrix(edges)
  for (r in 1:nrow(edges)) {
    am_wide[edges[r,1], edges[r,2]] <- 1
  }
  admat <- toLong(am_wide)
  admat <- reversedEdges(admat) %>%
    dplyr::mutate(reversed_connected=reversedConnection(.),
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

# -------------------------------------------------------------------------
# GRAPH TRAVERSAL & ACCESSORS
# -------------------------------------------------------------------------

bfsLong2 <- function(graph_G) {
  # Returns vector of nodes connected to root
  graph_G$parent <- as.character(graph_G$parent)
  children <- graph_G[(graph_G$parent == "root"), ]$child
  nodes <- c("root", children)
  
  while(length(children) > 0) {
    c <- children[1]
    temp.children <- graph_G[(graph_G$parent == c), ]$child
    
    if (any(temp.children %in% nodes)) {
      temp.children <- temp.children[! temp.children %in% nodes]
    }
    
    children <- c(children, temp.children)
    nodes <- c(nodes, temp.children)
    children <- children[-1]
  }
  return(nodes)
}

bridgeTestBFS <- function(graph_G, edge_e) {
  node_to_check <- edge_e$child
  nodes_connected_to_root <- bfsLong2(graph_G)
  !(node_to_check %in% nodes_connected_to_root)
}

getEdges <- function(am.long) {
  am.long %>%
    dplyr::filter(connected == 1) %>%
    dplyr::mutate(parent = as.character(parent))
}

getChildren <- function(am.long, node) {
  edges <- am.long %>%
    dplyr::mutate(parent = as.character(parent)) %>%
    dplyr::filter(connected == 1) %>%
    dplyr::filter(parent == node)
  return(edges$child)
}

verticesInGraph <- function(tb) {
  unique(c(tb$parent, tb$child))
}

# -------------------------------------------------------------------------
# LOGIC & CONSTRAINT CHECKS
# -------------------------------------------------------------------------

checkEdge <- function(edge, mcf, thresh = 0.2) {
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

satisfiesSumCondition <- function(edges, w, thresh = 0.2) {
  edges$parent <- as.character(edges$parent)
  all_parents <- unique(edges$parent)
  
  for (p in all_parents) {
    if (p == "root") {
      parent_ccf <- rep(1, ncol(w))
    } else {
      parent_ccf <- w[as.numeric(p), ]
    }
    
    children <- as.numeric(dplyr::filter(edges, parent == p)$child)
    if (length(children) > 1) {
      children_ccf <- colSums(w[children, ,drop=FALSE])
    } else {
      children_ccf <- w[children, ]
    }
    
    diff <- children_ccf - parent_ccf
    if (any(diff > thresh)) return(FALSE)
  }
  return(TRUE)
}

# -------------------------------------------------------------------------
# GRAPH ANNOTATIONS (Reversal/Connectivity)
# -------------------------------------------------------------------------

updateGraphElements <- function(am) {
  am %>%
    dplyr::mutate(parent=factor(parent, levels=unique(parent))) %>%
    dplyr::mutate(reversed_connected=reversedConnection(.)) %>%
    dplyr::mutate(bi_directional=isBidirectional(.)) %>%
    dplyr::mutate(root_connected=isRootConnected(.))
}

reversedEdges <- function(am) {
  am %>%
    dplyr::filter(!is.na(connected)) %>%
    tidyr::unite("reverse_edge", c("child", "parent"), sep="->", remove=FALSE)
}

reversedConnection <- function(am) {
  connections <- setNames(am$connected, am$edge)
  reversed_connections <- connections[am$reverse_edge] %>% "["(!is.na(.))
  reversed <- setNames(rep(0, nrow(am)), am$reverse_edge)
  reversed[names(reversed_connections)] <- reversed_connections
  reversed
}

isBidirectional <- function(am) {
  am %>%
    dplyr::mutate(bi_directional=(reverse_edge %in% edge) & connected==1 & reversed_connected == 1) %>%
    dplyr::pull(bi_directional)
}

isRootConnected <- function(am) isParentConnected(am)[1]

isParentConnected <- function(am) {
  am %>%
    dplyr::mutate(parent=factor(parent, levels=unique(parent))) %>%
    dplyr::group_by(parent) %>%
    dplyr::summarize(n=sum(connected)) %>%
    dplyr::pull(n) > 0
}
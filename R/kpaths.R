# library(igraph)
# library(tidyverse)

k_shortest_yen_out <- function(graph, source, dest, k){
  dest %>% 
    map(~{
      k_shortest_yen(graph, source, .x, k)
    })
}
  

#'@return the shortest path as a list of vertices or NULL if there is no path between src and dest
shortest_path <- function(graph, src, dest){
  path <- suppressWarnings(get.shortest.paths(graph, from = src, to = dest))
  path <- map(path$vpath, names)
  if (length(path)==0) NULL else path
} 


#'@return the sum of the weights of all the edges in the given path
path_weight <- function(path, graph) sum(E(graph, path=path)$weight)

#'@description sorts a list of paths based on the weight of the path
sort_paths <- function(graph, paths){ 
  # browser()
  paths[paths %>% sapply(path_weight, graph) %>% order]
}

#'@description creates a list of edges that should be deleted
find_edges_to_delete <- function(A,i,rootPath){
  edgesToDelete <- NULL
  for (p in A){
    rootPath_p <- p[1:i]
    if (identical(rootPath_p, rootPath)){
      edge <- paste(p[i], ifelse(is.na(p[i+1]),p[i],p[i+1]), sep = '|')
      edgesToDelete[length(edgesToDelete)+1] <- edge
    }
  }
  unique(edgesToDelete)
}

#returns the k shortest path from src to dest
#sometimes it will return less than k shortest paths. This occurs when the max possible number of paths are less than k
k_shortest_yen <- function(graph, src, dest, k){
  # browser()
  if (identical(src, dest)) stop('src and dest can not be the same (currently)')
  
  #accepted paths
  # A <- list(shortest_path(graph, src, dest))
  A <- shortest_path(graph, src, dest)
  if (k == 1) return (A)
  
  
  # So I am going to need to rewrite this section probably 
  # Either mapping over the below stuff, 
  # or 
  # I may remove the mode = "out" condition and  do a one to
  # one mapping, then just do a nested loop
  
  
  #potential paths
  B <- list()
  
  for (k_i in 2:k){
    prev_path <- A[[k_i-1]] # Take the path found
    num_nodes_to_loop <- length(prev_path)-1
    for(i in 1:num_nodes_to_loop){ # Loop over each of the nodes, deleting and finding a new path  on the total graph
      spurNode <- prev_path[i] # spur might mean current node of inquiry???
      rootPath <- prev_path[1:i]
      
      edgesToDelete <- find_edges_to_delete(A, i,rootPath) # A  is overall path list, i is the index of the current node, root_paath is the nodes  from the  staart up unti the index
      t_g <- delete.edges(graph, edgesToDelete)
      # for (edge in edgesToDelete) t_g <- delete.edges(t_g, edge)
      
      spurPath <- shortest_path(t_g, spurNode, dest)
      # browser()
      if (!is.null(spurPath)){
        # browser()
        total_path <- unlist(c(rootPath[-i], spurPath))
        # if (!list(total_path) %in% B) B[length(B)+1] <- total_path
        # if (!purrr::has_element(B, total_path)) B[length(B)+1] <- total_path
        if (!purrr::has_element(B, total_path)){
          B <-  append(B, list(total_path))
        }
      }
    }
    # browser()
    if (length(B) == 0) break
    B <- sort_paths(graph, B)
    A[k_i] <- B[1]
    B <- B[-1]
  }
  A
}

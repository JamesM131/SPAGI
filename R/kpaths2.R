assign_scores <- function(paths, graph){
  scores <- purrr::map(paths, ~{
    path_length <- length(.x)
    weight_sum <- sum(E(graph, path=.x)$weight)
    score <- 1000*path_length - weight_sum
    return(score)
  })
  
  purrr::map2(paths, scores, ~{
    tmp_path <- .x
    attr(tmp_path, "score") <- .y
    return(tmp_path)
  })
}

k_shortest_paths <- function(graph, from, to, k){

  if(k == 1) {
    paths <- shortest_paths(graph, from, to, mode = "out") %>% 
      purrr::pluck("vpath") %>% 
      purrr::map(names)
    scored_paaths <- assign_scores(paths, graph)

    return(scored_paaths)
  } 
  # browser()
  # pb_deep <- progress::progress_bar$new(total = length(to))
  a <- to %>% 
    purrr::map(~{
      # pb_deep$tick()
      # browser()
      # browser()
      k.shortest.paths(graph, from, to = .x, k)
    }) %>% 
    purrr::map_depth(2, ~names(pluck(.x, "vert", 1))) %>%
    purrr::flatten() %>% 
    assign_scores(graph)
}

# find k shortest paths
k.shortest.paths <- function(graph, from, to, k){
  # browser()
  # first shortest path
  # browser()
  k0 <- get.shortest.paths(graph,from,to, output='both')
  
  if(length(k0$vpath[[1]]) <= 1){
    return(NULL)
  }
  
  # number of currently found shortest paths
  kk <- 1
  
  # list of alternatives
  variants <- list()
  # browser()
  
  # shortest variants
  shortest.variants <- list(list(g=graph, path=k0$epath, vert=k0$vpath, dist=shortest.paths(graph,from,to)))
  
  # until k shortest paths are found
  while(kk<k){
    # take last found shortest path
    last.variant <- shortest.variants[[length(shortest.variants)]]              
    
    # calculate all alternatives
    variants <- calculate.variants(variants, last.variant, from, to)
    if(length(variants) != 0) {
      # find shortest alternative
      sp <- select.shortest.path(variants)
      
      # add to list, increase kk, remove shortest path from list of alternatives
      shortest.variants[[length(shortest.variants)+1]] <- list(g=variants[[sp]]$g, path=variants[[sp]]$variants$path, vert=variants[[sp]]$variants$vert, dist=variants[[sp]]$variants$dist)
      variants <- variants[-sp]
    }
    kk <- kk+1
  }
  
  return(shortest.variants)
}

# found all alternative routes
calculate.variants <- function(variants, variant, from, to){
  # take graph from current path
  g <- variant$g
  # browser()
  # iterate through edges, removing one each iterations
  for (j in unlist(variant$path)){
    
    newgraph <- delete.edges(g, j) # remove adge
    sp <- get.shortest.paths(newgraph,from,to, output='both') # calculate shortest path
    spd <- shortest.paths(newgraph,from,to) # calculate length
    if (spd != Inf && length(sp$vpath[[1]]) > 2){ # the the path is found
      # add to list, unless it already contains the same path
      if (!contains.path(variants, sp$vpath)) {
        variants[[length(variants)+1]] <- list(g=newgraph, variants=list(path=sp$epath, vert=sp$vpath, dist=spd))
      }
    }
  }
  
  return(variants)
}

# does a list contain this path?
contains.path <- function(variants, variant){
  return( any( unlist( lapply( variants, function(x){ identical(x$variant$vert,variant) } ) ) ) )
}

# which path from the list is the shortest?
select.shortest.path <- function(variants){
  return( which.min( unlist( lapply( variants, function(x){x$variants$dist} ) ) ) )
}
library(spagi)
library(tidyverse)
library(igraph)
library(tidygraph)
library(foreach)
library(doAzureParallel)

## Get PPI data for the protein molecules of species "mmusculus".
# mm.ppi <- get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species = "mmusculus")
# ## Get PPI data for the protein molecules of species "hsapiens".
# hs.ppi <- get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species = "hsapiens")
## Now combine and get the filtered PPI and the RP and TF proteins of the combined filtered PPI
comb.ppi.result <- combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
## Generate the pathway path data using the comb.ppi.result and housekeeping.gene data sets
# pathway.path2 <- generate_pathway_path(ppi.result = comb.ppi.result, housekeeping.gene, num_paths = 2)
# pathway.path3 <- generate_pathway_path(ppi.result = comb.ppi.result, housekeeping.gene, num_paths = 3)
# pathway.path4 <- generate_pathway_path(ppi.result = comb.ppi.result, housekeeping.gene, num_paths = 4)


# Create your cluster if it does not exist; this takes a few minutes
setCredentials("credentials.json")
cluster <- makeCluster("cluster.json")

# Register your parallel backend 
registerDoAzureParallel(cluster) 

# Check that the nodes are running 
getDoParWorkers()


comb.ppi.result <- combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)

all.significant.filtered.ppi<-comb.ppi.result$PPI
RPs<-comb.ppi.result$RPs
TFs<-comb.ppi.result$TFs

##First calculate the weight from the score as (1000-score)
edge.weight<- 1000 - all.significant.filtered.ppi$score
##

##Now, add the weight column to the data frame
all.edges<-data.frame(all.significant.filtered.ppi[,1:2], "weight"=edge.weight)

#make the rownames null
rownames(all.edges)<-NULL

#####Now create a graph and generate the pathway paths
##make the graph data frame from the "all.edges"
g1 <- graph.data.frame(d = all.edges, directed = TRUE)

# big_list <- RPs %>% 
#   map(~{
#     k_shortest_paths(graph = g1, from = .x, to = TFs, k =  num_paths)
#     # k_shortest_yen(graph = g1, src = .x, dest =  TFs, k = `num_paths`)
#   })



big_list <- foreach::foreach(src = RPs, .errorhandling = "remove", .options.azure = list(enableCloudCombine = FALSE, autoDeleteJob = FALSE)) %dopar% {
  library(tidyverse)
  library(igraph)
  assign_scores <- function(paths, graph){
    scores <- map(paths, ~{
      path_length <- length(.x)
      weight_sum <- sum(E(graph, path=.x)$weight)
      score <- 1000*path_length - weight_sum
      return(score)
    })
    
    map2(paths, scores, ~{
      tmp_path <- .x
      attr(tmp_path, "score") <- .y
      return(tmp_path)
    })
  }
  
  k_shortest_paths <- function(graph, from, to, k){
    
    if(k == 1) {
      paths <- shortest_paths(graph, from, to, mode = "out") %>% 
        pluck("vpath") %>% 
        map(names)
      # browser()
      scored_paaths <- assign_scores(paths, graph)
      
      return(scored_paaths)
    } 
    # browser()
    # pb_deep <- progress::progress_bar$new(total = length(to))
    a <- to %>% 
      map(~{
        # pb_deep$tick()
        # browser()
        k.shortest.paths(graph, from, to = .x, k)
      }) %>% 
      map_depth(2, ~names(pluck(.x, "vert", 1))) %>%
      flatten() %>% 
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
  
  k_shortest_paths(graph = g1, from = src, to = TFs, k =4)
}
# Get back the main  job and continue processing
big_list <- doAzureParallel::getJobResult("job20190603072855")
# Flatten the list (removing two levels) and discard the 
pb2 <- progress::progress_bar$new(total = length(flatten((big_list))))
# browser()
medium_list <- big_list %>% 
  flatten() %>% 
  # flatten() %>% 
  discard(~{
    pb2$tick()
    path_length <- length(.x)
    path_length < 3 
  })
pb3 <- progress::progress_bar$new(total = length(medium_list))
small_list <- medium_list %>% 
  discard(~{
    pb3$tick()
    all(data.table::`%chin%`(.x, housekeeping.gene))
  })



### cleaning the returned lists
library(tidyverse)


groups <- small_list %>% 
  map_chr(~{
    paste(head(.x, 1), "-", tail(.x, 1))
  })

head(groups)

small_list %>% 
  discard(~{
    
  })

mapping_list <- split(small_list, factor(groups))

first <- mapping_list %>% map(1) %>% compact()
second <- mapping_list %>% map(2) %>% compact()
third <- mapping_list %>% map(3) %>% compact()
fourth <- mapping_list %>% map(4) %>% compact()

spagi_universe <- list(first, second, third, fourth)



## Interesting
# a <- tibble(groups)
# a %>% count(groups) %>% count(n)

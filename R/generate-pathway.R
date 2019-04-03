#' @title generate_pathway_path
#'
#' @description This function generates the background pathway path data using
#'   the list object returned by the function combine_mm_hs_ppi.
#'
#' @rdname generate_pathway_path
#' @name generate_pathway_path
#'
#' @details This function generates the background pathway path data using the
#'   list object returned by the function combine_mm_hs_ppi.
#'   
#'   Here maybe you will get some warnings, but these are only for not reachable paths for
#'   source to destination
#'
#' @return This function returns a list consisting of the pathway path data that
#'   will be used as background data for the SPAGI package.
#'
#' @param ppi.result The list object returned by the function combine_mm_hs_ppi.
#' @param housekeeping.gene A vector consisting of housekeeping genes. This data
#'   is loaded automatically with the package. This data was generated using the
#'   gene expression profiles of different cell types and/or tissues from the
#'   ENCODE human and mouse project.
#' @param max.path.length An integer number indicating the maximum path length.
#'   Default is 7.
#'
#' @importFrom igraph graph.data.frame shortest_paths
#'
#' @export
#'
#' @examples
#' ## Need two folder at working directory for downloading stringdb files for
#' ## each species - stringdb_mouse, stringdb_human. ## It takes some time to
#' ## download the data, and then can reuse the downloaded data. ## Here we will
#' ## use RP.protein, KN.protein, TF.protein protein parameters. These data are
#' ## automatically loaded with the package. You can modify these parameters. ## We
#' ## will generate PPI data for two species - "mmusculus" and "hsapiens" by
#' ## calling the function get_ppi_for_molecules two times. ## Then we will combine
#' ## these two PPI data sets by using the combine_mm_hs_ppi function that will be
#' ## used to generate the pathway path data. ## Finally we will generate the
#' ## pathway path data using the combined PPI data
#'
#' ## Get PPI data for the protein molecules of species "mmusculus".
#' mm.ppi <- get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species = "mmusculus")
#' ## Get PPI data for the protein molecules of species "hsapiens".
#' hs.ppi <- get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species = "hsapiens")
#' ## Now combine and get the filtered PPI and the RP and TF proteins of the combined filtered PPI
#' comb.ppi.result <- combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
#' ## Generate the pathway path data using the comb.ppi.result and housekeeping.gene data sets
#' pathway.path <- generate_pathway_path(ppi.result = comb.ppi.result, housekeeping.gene)
#' head(summary(pathway.path))
generate_pathway_path <- function(ppi.result, housekeeping.gene, max.path.length = 7) {
  ##### preprocess the ppi.result data
  ## Assign the result data to the objects
  all.significant.filtered.ppi <- ppi.result$PPI
  RPs <- ppi.result$RPs
  TFs <- ppi.result$TFs
  ##
  
  ## NOTE
  ## First, it is good to know that when looking up paths, igraph understands weights as costs,
  ## i.e. on edges with higher weight it costs more to travel,
  ## so it will consider shorter the paths with lower sum weight.
  ## It is easy to turn this into the opposite, here we will do the score (range 0 to 999) as weight by 1000-score.
  
  ## First calculate the weight from the score as (1000-score)
  edge.weight <- 1000 - all.significant.filtered.ppi$score
  ##
  
  ## Now, add the weight column to the data frame
  all.edges <- data.frame(all.significant.filtered.ppi[, 1:2], "weight" = edge.weight)
  # make the rownames null
  rownames(all.edges) <- NULL
  ##
  #####
  
  ##### Now create a graph and generate the pathway paths
  ## make the graph data frame from the "all.edges"
  g1 <- graph.data.frame(d = all.edges, directed = TRUE)
  ##
  
  ## Find all shortest paths for all source nodes ("RPs") to destination ("TFs")
  ## Here by default uses the Dijkstra's algorithm for weighted directed graph
  ll.all.path <- list()
  for (i in 1:length(RPs)) {
    browser()
    ll.all.path[[RPs[i]]] <- shortest_paths(g1, from = RPs[i], to = TFs, mode = "out")
    print(i)
  }
  ##
  
  ## Finding all the complete (RP-KN-...-KN-TF) paths
  ## Here we considered for at least 3 layers and maximum 7 layers by default
  ll.all.path.complete <- list()
  for (i in 1:length(ll.all.path)) {
    # this tmp variable is used for all filtered paths for a pathway
    tmp.individual.pathway.paths.clean.ind <- NULL
    # this loop to check each individual path for a pathway
    for (j in 1:length(ll.all.path[[i]]$vpath)) {
      if ((length(unlist(ll.all.path[[i]]$vpath[j])) >= 3) & (length(unlist(ll.all.path[[i]]$vpath[j])) <= max.path.length)) {
        tmp.individual.pathway.paths.clean.ind <- c(tmp.individual.pathway.paths.clean.ind, j)
      }
    }
    # this combines all the filtered paths for a pathway to the respective pathway
    ll.all.path.complete[[names(ll.all.path)[i]]] <- ll.all.path[[i]]$vpath[tmp.individual.pathway.paths.clean.ind]
  }
  ##
  
  ## take only the pathways that have at least one complete path
  ll.all.path.complete.exist <- list()
  for (i in 1:length(ll.all.path.complete)) {
    if (length(ll.all.path.complete[[i]]) != 0) {
      ll.all.path.complete.exist[[names(ll.all.path.complete)[i]]] <- ll.all.path.complete[[i]]
    }
  }
  ##
  
  ## get only the pathway path elements, i.e., the names
  pathway.path.all <- list()
  for (i in 1:length(ll.all.path.complete.exist)) {
    tmp.pathway.path <- lapply(ll.all.path.complete.exist[[i]], function(x) {
      return(names(x))
    })
    pathway.path.all[[names(ll.all.path.complete.exist)[i]]] <- tmp.pathway.path
  }
  ##
  #####
  
  
  
  ##### This section is for removing the paths where all elements are housekeeping genes
  ## get the pathway paths in which all elements are not hk genes
  pathway.path.specific <- list()
  for (i in 1:length(pathway.path.all)) {
    tmp.path.spec <- lapply(pathway.path.all[[i]], function(x) {
      if (!(all(x %in% housekeeping.gene == "TRUE"))) {
        return(x)
      }
    })
    pathway.path.specific[[names(pathway.path.all)[i]]] <- tmp.path.spec
  }
  ##
  
  ## take only the existing pathway paths without null paths
  pathway.path.specific.clean <- lapply(pathway.path.specific, function(x) {
    return(x[!(sapply(x, is.null))])
  })
  ##
  
  ## take only the pathways that have at least one path
  pathway.path.specific.clean.2 <- list()
  for (i in 1:length(pathway.path.specific.clean)) {
    if (length(pathway.path.specific.clean[[i]]) != 0) {
      pathway.path.specific.clean.2[[names(pathway.path.specific.clean)[i]]] <- pathway.path.specific.clean[[i]]
    }
  }
  ##
  
  ## Finally return the pathway.path.specific.clean.2 data
  ## This data will be used as background pathway path data
  return(pathway.path.specific.clean.2)
  ##
  #####
}

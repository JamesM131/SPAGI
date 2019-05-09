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
#' @param num_paths An integer representing how many paths to return for each R and TF 
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
generate_pathway_path<-function(ppi.result, housekeeping.gene, max.path.length=7, num_paths = 1){
  #####preprocess the ppi.result data
  ##Assign the result data to the objects
  all.significant.filtered.ppi<-ppi.result$PPI
  RPs<-ppi.result$RPs
  TFs<-ppi.result$TFs
  ##
  
  ##NOTE
  ##First, it is good to know that when looking up paths, igraph understands weights as costs,
  ##i.e. on edges with higher weight it costs more to travel,
  ##so it will consider shorter the paths with lower sum weight.
  ##It is easy to turn this into the opposite, here we will do the score (range 0 to 999) as weight by 1000-score.
  
  ##First calculate the weight from the score as (1000-score)
  edge.weight<- 1000 - all.significant.filtered.ppi$score
  ##
  
  ##Now, add the weight column to the data frame
  all.edges<-data.frame(all.significant.filtered.ppi[,1:2], "weight"=edge.weight)
  #make the rownames null
  rownames(all.edges)<-NULL
  ##
  #####
  #####Now create a graph and generate the pathway paths
  ##make the graph data frame from the "all.edges"
  g1 <- graph.data.frame(d = all.edges, directed = TRUE)
  ##
  pb <- progress::progress_bar$new(total = length(RPs))
  big_list <- RPs %>% 
    map(~{
      pb$tick()
      k_shortest_yen(g1, src = .x, dest =  TFs, k = num_paths)
    })
  
  # Flatten the list (removing two levels) and discard the 
  big_list %>% 
    flatten() %>% 
    flatten() %>% 
    discard(~(length(.x) < 3 || length(.x) > max.path.length)) %>% 
    discard(~{
      all(.x %in% housekeeping.gene)
    })
}
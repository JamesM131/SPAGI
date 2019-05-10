#' @title identify_active_pathway_path
#'
#' @description
#' This function identifies active pathway paths for RNA-seq gene expression profile. It utilizes background pathway path data to identify the active pathway paths.
#'
#' @rdname identify_active_pathway_path
#' @name identify_active_pathway_path
#'
#' @details
#' This function identifies active pathway paths for RNA-seq gene expression profile. It utilizes background pathway path data to identify the active pathway paths.
#'
#' @return This function returns a list of pathways where each sublist denote a path of the pathway.
#'
#' @param pathway.path A list with pathway path data where each sublist denotes a path of the pathway. This is used as background data.
#' @param processed.query.data A list with expressed query data where each sublist corresponds for each cell/tissue type.
#'
#' @importFrom data.table chmatch
#'
#' @export
#'
#' @examples
#' ## Here we will use "pathway.path" as background data from the SPAGI repository.
#' ## Also we will use "ROR1.data" as query RNA-seq gene expression data. This data is for ocular lens epithelial cell differentiated from human pluripotent stem cells.
#' ## These data sets are loaded automatically with the package.
#' 
#' 
#' comb.ppi.result <- combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
#' # Generate the pathway path data using the comb.ppi.result and housekeeping.gene data sets
#' new_background_paths <- generate_pathway_path(ppi.result = comb.ppi.result, housekeeping.gene)
#' 
#' 
#' ## Pre-process the query data (ROR1.data), the data has already been made in CPM and log2 normalized format. Also we have already made the replicate names same for the data.
#' ROR1.processed.data <- preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
#' ## Identify active pathway paths of the processed query data
#' ROR1.active.pathway <- identify_active_pathway_path(pathway.path = new_background_paths, processed.query.data = ROR1.processed.data)
#' head(ROR1.active.pathway$ROR1_LEC$FGFR1)

identify_active_pathway_path <- function(pathway.path, processed.query.data) {
  # A note on my (James) change here: 
  # The previous version of this function exhibited an unexpected/undocumented
  # behaviour in that it would return 1 list of active proteins for each column
  # in the original (non-processed RNA data). An exception to this behaviour
  # was when the columns were given the same name (as with the ROR1.data) in
  # which only one list of proteins was returned. This new version of the
  # function does not allow this silent combination (although as a result of
  # this, some of the behaviours regarding the ROR1 example may change).

  gene_names <- processed.query.data %>% 
    map(~sort(stringr::str_to_upper(names(.x))))
  
  pb <- progress::progress_bar$new(total = length(processed.query.data)*length(pathway.path))
  
  processed.query.data %>% 
    map(~sort(stringr::str_to_upper(names(.x)))) %>% 
    map(function(names){
      pathway.path %>% 
        keep(~{
          pb$tick()
          all(data.table::`%chin%`(.x, names))
        })
    })
}



#' @title get_pathway_activity_score
#'
#' @description
#' This function generates pathway activity score of the active pathways for each cell/tissue type. It uses active pathway path and processed query data with a high expression threshold to generate the ranking metric.
#'
#' @rdname get_pathway_activity_score
#' @name get_pathway_activity_score
#'
#' @details
#' This function generates pathway activity score of the active pathways for each cell/tissue type. It uses active pathway path and processed query data with a high expression threshold to generate the ranking metric.
#'
#' @return This function returns a list of pathway activity score for each cell/tissue type.
#'
#' @param active.pathway.path A list of active pathway path data for each cell/tissue type as returned by the function 'identify_active_pathway_path'.
#' @param processed.query.data A list with expressed query data where each sublist corresponds for each cell/tissue type as returned by the function 'preprocess_querydata'.
#' @param high.exp.th A high expression threshold value for the processed query data. It is used to get active (i.e., highly expressed) molecule proportion for each pathway path.
#'
#' @export
#'
#' @examples
#' ## Here we will use "pathway.path" as background data from the SPAGI repository.
#' ## Also we will use "ROR1.data" as query RNA-seq gene expression data. This data is for ocular lens epithelial cell differentiated from human pluripotent stem cells.
#' ## These data sets are loaded automatically with the package.
#' 
#' ## Pre-process the query data (ROR1.data), the data has already been made in CPM and log2 normalized format. Also we have already made the replicate names same for the data.
#' ROR1.processed.data <- preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
#' ## Identify active pathway paths of the processed query data
#' ROR1.active.pathway <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)
#' ## Get activity score of the active pathways
#' ROR1.pathway.activity.score <- get_pathway_activity_score(active.pathway.path = ROR1.active.pathway, processed.query.data = ROR1.processed.data, high.exp.th = 7)
#' head(ROR1.pathway.activity.score$ROR1_LEC)
get_pathway_activity_score <- function(active.pathway.path, processed.query.data, high.exp.th) {
  ## process separately each cell/tissue to get active pathway ranking metric
  pathway.activity.score <- list()
  for (i in 1:length(active.pathway.path)) {
    # get each cell/tissue active pathway paths
    each.cell.active.pathway.path <- active.pathway.path[[i]]
    
    # take the cell/tissue name from the pathway to get that cell/tissue processed.query.data and high.exp.th
    tmp.cell.name <- names(active.pathway.path[i])
    
    # take the respective cell/tissue processed data
    tmp.cell.processed.data <- processed.query.data[[tmp.cell.name]]
    
    ## get the pathway paths average active gene count proportion
    pathway.path.active.gene.count.proportion <- lapply(each.cell.active.pathway.path, function(y) {
      individual.pathway.active.gene.count <- lapply(y, function(x) {
        # get active (i.e., highly expressed) gene count proportion for each path and return
        tmp.path.gene.exp <- tmp.cell.processed.data[unlist(x)]
        tmp.path.active.genes <- names(tmp.path.gene.exp[tmp.path.gene.exp >= high.exp.th])
        return(length(tmp.path.active.genes) / length(unlist(x)))
      })
      # calcaulate average active gene count proportion for each pathway and return
      return(sum(unlist(individual.pathway.active.gene.count)) / length(y))
    })
    ##
    
    # assign average active gene count proportion for each cell/tissue
    pathway.activity.score[[tmp.cell.name]] <- unlist(pathway.path.active.gene.count.proportion)
  }
  # return average active gene count proportion for all cell/tissue
  return(pathway.activity.score)
}

#' @title get_pathway_downstream_tf_number
#'
#' @description
#' This function calculates the number of downstrean transcription factors of each active pathway for each cell/tissue type. It uses the active pathway path object to calculate the number of downstrean transcription factor.
#'
#' @rdname get_pathway_downstream_tf_number
#' @name get_pathway_downstream_tf_number
#'
#' @details
#' This function calculates the number of downstrean transcription factors of each active pathway for each cell/tissue type. It uses the active pathway path object to calculate the number of downstrean transcription factor.
#'
#' @return This function returns a list of number of downstrean transcription factors for each cell/tissue type.
#'
#' @param active.pathway.path A list of active pathway path data for each cell/tissue type as returned by the function 'identify_active_pathway_path'.
#'
#' @export
#'
#' @examples
#' ## Here we will use "pathway.path" as background data from the SPAGI repository.
#' ## Also we will use "ROR1.data" as query RNA-seq gene expression data. This data is for ocular lens epithelial cell differentiated from human pluripotent stem cells.
#' ## These data sets are loaded automatically with the package.
#' 
#' ## Pre-process the query data (ROR1.data), the data has already been made in CPM and log2 normalized format. Also we have already made the replicate names same for the data.
#' ROR1.processed.data <- preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
#' ## Identify active pathway paths of the processed query data
#' ROR1.active.pathway <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)
#' ## Get the number of downstream transcription factors of the active pathways
#' ROR1.pathway.downstream.tf.count <- get_pathway_downstream_tf_number(active.pathway.path = ROR1.active.pathway)
#' head(ROR1.pathway.downstream.tf.count$ROR1_LEC)
get_pathway_downstream_tf_number <- function(active.pathway.path) {
  ## get each pathway number of downstream transcription factors for each cell type
  pathway.downstream.tf.count <- lapply(active.pathway.path, function(y) {
    # here, y is each cell type pathway path
    tmp.downstream.tf.count <- lapply(y, function(x) {
      # here x is each pathway for that cell type
      return(length(x))
    })
    return(unlist(tmp.downstream.tf.count))
  })
  return(pathway.downstream.tf.count)
}


#' @title get_pathway_ranking_metric
#'
#' @description This function generates pathway activity score and number of
#' downstream transcription factors of the active pathways for each cell/tissue
#' type. It uses active pathway path and processed query data with a high
#' expression threshold to generate the activity score. However, it uses only
#' the active pathway path data to calculate the number of downstream
#' transcription factors.
#'
#' @rdname get_pathway_ranking_metric
#' @name get_pathway_ranking_metric
#'
#' @details This function generates pathway activity score and number of
#' downstream transcription factors of the active pathways for each cell/tissue
#' type. It uses active pathway path and processed query data with a high
#' expression threshold to generate the activity score. However, it uses only
#' the active pathway path data to calculate the number of downstream
#' transcription factors.
#'
#' @return This function returns a list of sublist with pathway activity score
#'   and the number of downstream transcription factors for each cell/tissue
#'   type.
#'
#' @param active.pathway.path A list of active pathway path data for each
#'   cell/tissue type as returned by the function
#'   'identify_active_pathway_path'.
#' @param processed.query.data A list with expressed query data where each
#'   sublist corresponds for each cell/tissue type as returned by the function
#'   'preprocess_querydata'.
#' @param high.exp.th A high expression threshold value for the processed query
#'   data. It is used to get active (i.e., highly expressed) molecule proportion
#'   for each pathway path.
#'
#' @export
#'
#' @examples
#' ## Here we will use "pathway.path" as background data from the SPAGI repository.
#' ## Also we will use "ROR1.data" as query RNA-seq gene expression data. This data is for ocular lens epithelial cell differentiated from human pluripotent stem cells.
#' ## These data sets are loaded automatically with the package.
#'
#' ## Pre-process the query data (ROR1.data), the data has already been made in CPM and log2 normalized format. Also we have already made the replicate names same for the data.
#' ROR1.processed.data <- preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
#' ## Identify active pathway paths of the processed query data
#' ROR1.active.pathway <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)
#' ## Get active pathway ranking metric (i.e., activity score and number of downstream transcription factors)
#' ROR1.active.pathway.ranking.metric <- get_pathway_ranking_metric(active.pathway.path = ROR1.active.pathway, processed.query.data = ROR1.processed.data, high.exp.th = 7)
#' head(ROR1.active.pathway.ranking.metric$activity.score$ROR1_LEC)
#' head(ROR1.active.pathway.ranking.metric$downstream.tf.count$ROR1_LEC)
get_pathway_ranking_metric <- function(active.pathway.path, processed.query.data, high.exp.th) {
  ## First get the activity score of each pathway for each cell/tissue type
  activity.score <- get_pathway_activity_score(active.pathway.path = active.pathway.path, processed.query.data = processed.query.data, high.exp.th = high.exp.th)
  
  ## Next get the number of downstream transcription factors of each pathway for each cell/tissue type
  downstream.tf.count <- get_pathway_downstream_tf_number(active.pathway.path = active.pathway.path)
  
  ## Finally combine the above results in a list and return
  pathway.ranking.metric <- list()
  pathway.ranking.metric[["activity.score"]] <- activity.score
  pathway.ranking.metric[["downstream.tf.count"]] <- downstream.tf.count
  
  return(pathway.ranking.metric)
}


#' @title display_pathway_ranking_metric
#'
#' @description
#' This function plots the pathway ranking metric for each cell type or tissue in a 2D plane.
#'
#' @rdname display_pathway_ranking_metric
#' @name display_pathway_ranking_metric
#'
#' @details
#' This function plots the pathway ranking metric for each cell type or tissue in a 2D plane.
#'
#' @return NULL
#'
#' @param pathway.ranking.metric The ranking metric result returned by 'get_pathway_ranking_metric' function.
#'
#' @export
#'
#' @examples
#' ## Here we will use "pathway.path" as background data from the SPAGI repository.
#' ## Also we will use "ROR1.data" as query RNA-seq gene expression data. This data is for ocular lens epithelial cell differentiated from human pluripotent stem cells.
#' ## These data sets are loaded automatically with the package.
#' 
#' ## Pre-process the query data (ROR1.data), the data has already been made in CPM and log2 normalized format.
#' ROR1.processed.data <- preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
#' ## Identify active pathway paths of the processed query data
#' ROR1.active.pathway <- identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)
#' ## Get active pathway ranking metric (i.e., activity score and number of downstream transcription factors)
#' ROR1.active.pathway.ranking.metric <- get_pathway_ranking_metric(active.pathway.path = ROR1.active.pathway, processed.query.data = ROR1.processed.data, high.exp.th = 7)
#' ## Plot the ranking metric result (i.e., activity score and number of downstream transcription factors) in a 2D plane
#' display_pathway_ranking_metric(pathway.ranking.metric = ROR1.active.pathway.ranking.metric)
#' ## To separate the top ranked pathways we can do this
#' abline(v = 45, h = 0.2, lty = 2, col = "black")
display_pathway_ranking_metric <- function(pathway.ranking.metric) {
  # in each loop, the result of one query cell/tissue type is processed and plotted in a 2D plane
  for (i in 1:length(pathway.ranking.metric$activity.score)) {
    # get the name of each cell type
    cell.tissue.names <- names(pathway.ranking.metric$activity.score[i])
    
    ## for setting the title
    if (nchar(cell.tissue.names) > 40) {
      title <- paste("The result of:\n", cell.tissue.names, sep = " ")
    }
    else {
      title <- paste("The result of", cell.tissue.names, sep = " ")
    }
    
    ## plot the result in a 2D plane - number of downstream TF in x axis and activity score in y axis
    plot(pathway.ranking.metric$downstream.tf.count[[i]], pathway.ranking.metric$activity.score[[i]],
         xlim = c(0, 200), ylim = c(0, 1.0), type = "n", bty = "n", main = title,
         xlab = "Number of downstream transcription factor", ylab = "Pathway activity score"
    )
    text(pathway.ranking.metric$downstream.tf.count[[i]], pathway.ranking.metric$activity.score[[i]],
         names(pathway.ranking.metric$downstream.tf.count[[i]]),
         cex = 0.5, col = "black"
    )
    
    # for console message for each cell/tissue type ploting
    consol.msg <- paste(cell.tissue.names, "-- result plotting done!!", sep = " ")
    print(consol.msg)
  }
}

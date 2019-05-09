#' @title preprocess_querydata
#'
#' @description
#' This function preprocesses the query data to calculate average value of the replicates and then get expressed genes based on expression cut-off.
#'
#' @rdname preprocess_querydata
#' @name preprocess_querydata
#'
#' @details
#' This function preprocesses the query data to calculate average value of the replicates and then get expressed genes based on expression cut-off.
#'
#' @return This function returns a list with specifically expressed genes for each cell type / tissue
#'
#' @param cell.tissue.data An expression matrix with replicated column headers per-replicate. It is assumed that all query data are in RPKM/FPKM/CPM and log normalized form. Also assume that gene ids are official gene symbols. For the matrix, rows denote the genes and the columns denote the cell types or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue.
#' @param exp.cutoff.th An expression cut-off threshold value for the query data.
#' @param species The species abbreviation of the query data (cell.tissue.data). Default is "hsapiens".
#' @param data.format Format of cell.tissue.data. Default is "matrix".
#' @param experiment.descriptor A vector corresponding to the matrix column names of cell.tissue.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to NULL.
#'
#' @export
#'
#' @examples
#' query.data <- matrix(sample(1:10, 100, replace = TRUE), 10, 10)
#' rownames(query.data) <- c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' colnames(query.data) <- c("cell1", "cell1", "cell1", "cell2", "cell2", "cell2", "cell3", "cell3", "cell3", "cell3")
#' preprocess_querydata(cell.tissue.data = query.data, exp.cutoff.th = 5)
preprocess_querydata <- function(cell.tissue.data, exp.cutoff.th, species = "hsapiens", data.format = "matrix", experiment.descriptor = NULL) {
  # browser()
  if (data.format == "matrix") {
    ## to make all genes uppercase for mmusculus as PPI data are all uppercase
    if (species == "mmusculus") {
      rownames(cell.tissue.data) <- toupper(rownames(cell.tissue.data))
    }
    else if (species == "hsapiens") {
      cell.tissue.data <- cell.tissue.data
    }
    else {
      print("ERROR: could not support other species at this moment!!")
      return(NULL)
    }
    ##
    
    
    ## process the query data
    # calculate average value for each gene based on replicates for each cell type
    if (is.null(experiment.descriptor)) experiment.descriptor <- colnames(cell.tissue.data)
    cell.data.formatted <- format_matrix_data(cell.tissue.data, experiment.descriptor)
    
    # get expressed genes for each cell based on cut-off threshold
    expressed.cell.data <- list()
    cell.names <- colnames(cell.data.formatted)
    for (i in 1:ncol(cell.data.formatted)) {
      each.cell.data <- cell.data.formatted[, i]
      names(each.cell.data) <- rownames(cell.data.formatted)
      # now get expressed genes based on cut-off value
      each.cell.exp.data <- each.cell.data[each.cell.data >= exp.cutoff.th]
      expressed.cell.data[[cell.names[i]]] <- each.cell.exp.data
    }
    
    # return the list with expressed genes data for each cell
    return(expressed.cell.data)
    ##
  }
  else {
    print("ERROR: could not support other data format at this moment!!")
    return(NULL)
  }
}

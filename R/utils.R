#' @title compute_mean
#'
#' @description
#' This function calculates average expression values of multiple replicate samples
#'
#' @rdname compute_mean
#' @name compute_mean
#'
#' @details This function calculates average expression values of multiple
#' replicate samples using rowMeans. If the data has no replicates (a single
#' vector or a matrix with only one column) the original values are returned.
#'
#' @return This function returns a vector with the mean expression values.
#'
#' @param data.matrix A matrix with gene expression values, rows denote genes
#'   and columns denote replicate samples. Alternatively a vector of gene
#'   expression values.
#'   
#' @export
#'
#' @examples
#' dd <- matrix(sample(1:10, 30, replace = TRUE), 10, 3)
#' rownames(dd) <- c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' compute_mean(dd)
compute_mean <- function(data.matrix) {
  if ((is.vector(data.matrix) == TRUE) || (ncol(data.matrix) < 2)) {
    data.matrix <- cbind(val1 = data.matrix, vall2 = data.matrix)
  }
  return(rowMeans(data.matrix))
}

#' @title format_matrix_data
#'
#' @description
#' This function re-formats gene expression matrix data by averaging replicate samples. The cell type identifiers should be the column names of the expression matrix, otherwise the experiment.descriptor parameter is provided to specify cell-types of the samples (i.e. from a separate metadata file).
#'
#' @rdname format_matrix_data
#' @name format_matrix_data
#'
#' @details
#' This function re-formats gene expression matrix data by averaging replicate samples. The cell type identifiers should be the column names of the expression matrix, otherwise the experiment.descriptor parameter is provided to specify cell-types of the samples (i.e. from a separate metadata file).
#'
#' @return This function returns a matrix with mean value of the replicates for each cell type.
#'
#' @param matrix.data A matrix with gene expression profiles where rows denote the genes and the columns denote the cells and/or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue.
#' @param experiment.descriptor A vector corresponding to the columns of matrix.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to the column names of the matrix.data.
#'
#' @export
#'
#' @examples
#' mm <- matrix(sample(1:10, 100, replace = TRUE), 10, 10)
#' rownames(mm) <- c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' colnames(mm) <- c("cell2", "cell3", "cell1", "cell2", "cell1", "cell3", "cell3", "cell2", "cell3", "cell1")
#' format_matrix_data(mm)
format_matrix_data <- function(matrix.data, experiment.descriptor = colnames(matrix.data)) {
  colnames(matrix.data) <- experiment.descriptor
  
  # sorting the matrix data sets according to the column names
  ordered.column.names <- colnames(matrix.data)[order(colnames(matrix.data))]
  matrix.data.ordered <- matrix.data[, order(colnames(matrix.data))]
  colnames(matrix.data.ordered) <- ordered.column.names
  
  repCellTissue <- table(colnames(matrix.data.ordered))
  nCellTissue <- length(repCellTissue)
  nameCellTissue <- names(repCellTissue)
  
  
  # processing the matrix data to make a new matrix with number of rows same as before, but
  # number of columns equal to the number of cells/tissues by calculating the average value of the replicates
  if (ncol(matrix.data.ordered) == nCellTissue) {
    f.mat <- matrix.data.ordered # here each cell/tissue has only 1 replicate, with existing cell names
  }
  else { # here each cell/tissue has one or more replicates
    f.mat <- matrix(0, nrow = nrow(matrix.data.ordered), ncol = nCellTissue)
    r <- repCellTissue
    s <- 1
    e <- r[1]
    for (i in 1:nCellTissue) {
      mean.data <- compute_mean(matrix.data.ordered[, s:e])
      f.mat[, i] <- mean.data
      if (i != nCellTissue) {
        s <- s + r[i]
        e <- e + r[i + 1]
      }
    }
    rownames(f.mat) <- rownames(matrix.data.ordered)
    colnames(f.mat) <- nameCellTissue
  }
  
  return(f.mat)
}

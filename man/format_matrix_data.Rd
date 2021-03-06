% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{format_matrix_data}
\alias{format_matrix_data}
\title{format_matrix_data}
\usage{
format_matrix_data(matrix.data,
  experiment.descriptor = colnames(matrix.data))
}
\arguments{
\item{matrix.data}{A matrix with gene expression profiles where rows denote the genes and the columns denote the cells and/or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue.}

\item{experiment.descriptor}{A vector corresponding to the columns of matrix.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to the column names of the matrix.data.}
}
\value{
This function returns a matrix with mean value of the replicates for each cell type.
}
\description{
This function re-formats gene expression matrix data by averaging replicate samples. The cell type identifiers should be the column names of the expression matrix, otherwise the experiment.descriptor parameter is provided to specify cell-types of the samples (i.e. from a separate metadata file).
}
\details{
This function re-formats gene expression matrix data by averaging replicate samples. The cell type identifiers should be the column names of the expression matrix, otherwise the experiment.descriptor parameter is provided to specify cell-types of the samples (i.e. from a separate metadata file).
}
\examples{
mm <- matrix(sample(1:10, 100, replace = TRUE), 10, 10)
rownames(mm) <- c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
colnames(mm) <- c("cell2", "cell3", "cell1", "cell2", "cell1", "cell3", "cell3", "cell2", "cell3", "cell1")
format_matrix_data(mm)
}

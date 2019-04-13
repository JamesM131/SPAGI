#' @title SPAGI: Identification of active signalling pathways by integrating gene expression and protein interaction network
#'
#' @description SPAGI is an R package for active pathway identification for RNA-seq gene expression profiles. This package contains the neccessary R code to run SPAGI as described in "SPAGI: Identification of active signalling pathways using RNA-seq gene expression profiles". SPAGI is implemented within a helpful framework to identify active pathway paths from RNA-seq gene expression profiles.
#'
#'
#' @author Md Humayun Kabir <humayun.mkabir@gmail.com>
#'
#' @docType package
#' @name spagi-package
#' @aliases spagi SPAGI
#'
#' @examples
#' ## Do a sample analysis using human ocular lens epithelial cell (LEC) RNA-seq gene expression data.
#' 
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
#' ## Plot the ranking metric result (i.e., activity score and number of downstream transcription factors) in a 2D plane
#' display_pathway_ranking_metric(pathway.ranking.metric = ROR1.active.pathway.ranking.metric)
#' ## To separate the top ranked pathways we can do this
#' abline(v = 45, h = 0.2, lty = 2, col = "black")
NULL






##################################################################################.
# Copyright Victor Chang Cardiac Research Institute 2018
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################.






######################################################################################.
# Users should also check the names of functions loaded by this script
#
#
###############################################################.
# Please note that it is necessary to pre-process the RNA-seq query data to calculate RPKM/FPKM/CPM of raw read count data and make log2 normalization before utilizing these data with the SPAGI package.
# Also note that the background pathway path and housekeeping gene data are in official gene symbol format. So please make your gene ids as official gene symbols before using with SPAGI package.
# The SPAGI package does not perform any normalizations. It assumes that all query data are in RPKM/FPKM/CPM and log2 normalized format.
# To utilize the SPAGI package, you have to provide an expression cut-off threshold and a high expression threshold (i.e., an expression value that is high enough) for your query data.
###############################################################.

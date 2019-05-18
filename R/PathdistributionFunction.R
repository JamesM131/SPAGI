#function that loops through a nested list and gets the length of each set of elements in the list
#Pathway is a list from the spagi library where each 
distpathlength <- function(pathway)
{
  myvector <- c()
  for(i in 1:length(pathway))
  {
    for(j in 1:length(pathway[[i]]))
    {
      myvector <- c(myvector,(length(unlist(pathway[[i]][j]))))
    }
  }
  return(myvector)
}


library(spagi)
#Generate pathway with greater than 7 paths
mm.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="mmusculus")
#' ## Get PPI data for the protein molecules of species "hsapiens".
hs.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="hsapiens")
#' ## Now combine and get the filtered PPI and the RP and TF proteins of the combined filtered PPI
comb.ppi.result<-combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
#' ##Generate the pathway path data using the comb.ppi.result and housekeeping.gene data sets
#Original cut-off
pathway.path7 <- generate_pathway_path(ppi.result = comb.ppi.result, housekeeping.gene)

#effectively no cut-off
pathway.path100<-generate_pathway_path(ppi.result=comb.ppi.result, housekeeping.gene, max.path.length = 100)

#dataframe of the distribution of the path lengths
dist7 <- as.data.frame(distpathlength(pathway.path7))
dist100 <- as.data.frame(distpathlength(pathway.path100))


#Count of each number of paths 
count7 <- as.data.frame(count(dist7, `distpathlength(pathway.path7)`))
count100 <- as.data.frame(count(dist100,`distpathlength(pathway.path100)`))

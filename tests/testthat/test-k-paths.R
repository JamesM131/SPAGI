context("test-k-paths")



test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})



# test_that("shortest 1 path is a list of length 4", {
#   path <- shortest_paths(g1, from = RPs[i], to = TFs, mode = "out")
#   expect_is(path, "list")
#   expect_length(path, 4)
# })
comb.ppi.result <- combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
## Generate the pathway path data using the comb.ppi.result and housekeeping.gene data sets
pathway.path <- generate_pathway_path(ppi.result = comb.ppi.result, housekeeping.gene)

test_that("There are 89161 backgronud paths", {
  # mm.ppi <- get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species = "mmusculus")
  # ## Get PPI data for the protein molecules of species "hsapiens".
  # hs.ppi <- get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species = "hsapiens")
  # ## Now combine and get the filtered PPI and the RP and TF proteins of the combined filtered PPI
  expect_length(pathway.path, 89161)
})

test_that("Active Path Ouput Remains the same", {
  #Here I need to use spagi::pathway.path  and  the ROR1 data to get the identification pathway result. This can  be used to access the
  ROR1.processed.data <- preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
  ROR1.active.pathway <- identify_active_pathway_path(pathway.path, processed.query.data = ROR1.processed.data)
  
  expect_length(pathway.path, 23665)
})



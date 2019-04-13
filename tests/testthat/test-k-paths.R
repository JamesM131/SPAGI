context("test-k-paths")



test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})



test_that("shortest 1 path is a list of length 4", {
  path <- shortest_paths(g1, from = RPs[i], to = TFs, mode = "out")
  expect_is(path, "list")
  expect_length(path, 4)
})
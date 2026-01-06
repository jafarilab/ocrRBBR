# tests/testthat/test-data.R
test_that("multiome_human_mouse loads correctly", {
  data("multiome_human_mouse", package = "ocrRBBR")
  expect_true(exists("mouse_atacseq_data"))
  expect_true(exists("human_rnaseq_data"))
})
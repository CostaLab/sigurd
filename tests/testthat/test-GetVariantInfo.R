test_that("Testing GetVariantInfo.R", {
  # Loading the input object.
  input_object <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  test_consensus <- GetVariantInfo(SE = input_object, information = "consensus", variants = "chrM_1_G_A", cells = NULL)
  test_fraction  <- GetVariantInfo(SE = input_object, information = "fraction",  variants = "chrM_1_G_A", cells = "Test_Cell_1")
  test_coverage  <- GetVariantInfo(SE = input_object, information = "coverage",  variants = "chrM_1_G_A", cells = c("Test_Cell_1", "Test_Cell_3"))
  test_alts      <- GetVariantInfo(SE = input_object, information = "alts",      variants = "chrM_1_G_A", cells = NULL)
  test_refs      <- GetVariantInfo(SE = input_object, information = "refs",      variants = "chrM_1_G_A", cells = NULL)
  # We generate the expected output.
  expected_output_consensus <- Matrix::sparseMatrix(i = c(1,1,1,1), j = 1:4, x = rep(3,4), dims = c(1,4), dimnames = list("chrM_1_G_A", paste0("Test_Cell_", 1:4)))
  expected_output_fraction  <- Matrix::sparseMatrix(i = 1, j = 1, x = 0.07168459, dims = c(1,1), dimnames = list("chrM_1_G_A", "Test_Cell_1"))
  expected_output_coverage  <- Matrix::sparseMatrix(i = c(1,1), j = 1:2, x = c(279,230), dims = c(1,2), dimnames = list("chrM_1_G_A", paste0("Test_Cell_", c(1,3))))
  expected_output_alts      <- Matrix::sparseMatrix(i = c(1,1,1,1), j = 1:4, x = c(20,50,50,10), dims = c(1,4), dimnames = list("chrM_1_G_A", paste0("Test_Cell_", 1:4)))
  expected_output_refs      <- Matrix::sparseMatrix(i = c(1,1,1,1), j = 1:4, x = c(157,25,75,116), dims = c(1,4), dimnames = list("chrM_1_G_A", paste0("Test_Cell_", 1:4)))
  # We perform the tests.
  expect_equal(test_consensus, expected_output_consensus, tolerance = 1e-6)
  expect_equal(test_fraction, expected_output_fraction, tolerance = 1e-6)
  expect_equal(test_coverage, expected_output_coverage, tolerance = 1e-6)
  expect_equal(test_alts, expected_output_alts, tolerance = 1e-6)
  expect_equal(test_refs, expected_output_refs, tolerance = 1e-6)
})

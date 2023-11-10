test_that("Testing getAltMatrix.R", {
  # Loading the test SE object.
  input_object <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_Data.rds"))
  test_A <- getAltMatrix(SE_object = input_object, letter = "A", chromosome_prefix = "chrM")
  test_T <- getAltMatrix(SE_object = input_object, letter = "T", chromosome_prefix = "chrM")
  test_C <- getAltMatrix(SE_object = input_object, letter = "C", chromosome_prefix = "chrM")
  test_G <- getAltMatrix(SE_object = input_object, letter = "G", chromosome_prefix = "chrM")
  # We generate the expected results.
  expected_result_A <- Matrix::sparseMatrix(i = c(1,2,3,1,1:3,1:3), j = c(rep(1,3),2,rep(3,3),rep(4,3)), x = c(20,60,80,50,50,50,50,10,20,60), dims = c(3,4),
                                            dimnames = list(c("chrM_1_G>A", "chrM_3_T>A", "chrM_4_C>A"), paste0("Cell_", 1:4)))
  expected_result_T <- Matrix::sparseMatrix(i = rep(1:3, 4), j = rep(1:4,each=3), x = c(69,145,93,112,84,98,41,166,85,6,86,98), dims = c(3,4),
                                            dimnames = list(c("chrM_1_G>T", "chrM_2_A>T", "chrM_4_C>T"), paste0("Cell_", 1:4)))
  expected_result_C <- Matrix::sparseMatrix(i = rep(1:3, 4), j = rep(1:4,each=3), x = c(33,76,36,141,141,63,64,35,89,168,48,91), dims = c(3,4),
                                            dimnames = list(c("chrM_1_G>C", "chrM_2_A>C", "chrM_3_T>C"), paste0("Cell_", 1:4)))
  expected_result_G <- Matrix::sparseMatrix(i = rep(1:3, 4), j = rep(1:4,each=3), x = c(115,153,36,76,148,168,87,50,123,110,104,110), dims = c(3,4),
                                            dimnames = list(c("chrM_2_A>G", "chrM_3_T>G", "chrM_4_C>G"), paste0("Cell_", 1:4)))
  expect_equal(test_A, expected_result_A, tolerance = 1e-6)
  expect_equal(test_T, expected_result_T, tolerance = 1e-6)
  expect_equal(test_C, expected_result_C, tolerance = 1e-6)
  expect_equal(test_G, expected_result_G, tolerance = 1e-6)
})

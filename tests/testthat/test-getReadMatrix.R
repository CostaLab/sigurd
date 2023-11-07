test_that("Testing getReadMatrix.R", {
  # We load the input object.
  input_object <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_Data.rds"))
  test_A <- sigurd::getReadMatrix(SE = input_object, letter = "A", chromosome_prefix = "chrM")
  test_C <- sigurd::getReadMatrix(SE = input_object, letter = "C", chromosome_prefix = "chrM")
  test_G <- sigurd::getReadMatrix(SE = input_object, letter = "G", chromosome_prefix = "chrM")
  test_T <- sigurd::getReadMatrix(SE = input_object, letter = "T", chromosome_prefix = "chrM")
  # We generate the expected results.
  expected_result_A <- Matrix::sparseMatrix(i = c(1:4,1:2,rep(1:4,2)), j = c(rep(1,4),rep(2,2),rep(3,4),rep(4,4)), x = c(20,40,60,80,50,50,50,50,50,50,10,30,20,60),                 dims = c(4,4), dimnames = list(c("chrM_1_G_A", "chrM_2_A_A", "chrM_3_T_A", "chrM_4_C_A"), paste0("Cell_", 1:4)))
  expected_result_C <- Matrix::sparseMatrix(i = rep(1:4,4),            j = rep(1:4,each=4),                        x = c(33,76,36,114,141,141,63,168,64,35,89,89,168,48,91,114),     dims = c(4,4), dimnames = list(c("chrM_1_G_C", "chrM_2_A_C", "chrM_3_T_C", "chrM_4_C_C"), paste0("Cell_", 1:4)))
  expected_result_G <- Matrix::sparseMatrix(i = rep(1:4,4),            j = rep(1:4,each=4),                        x = c(157,115,153,36,25,76,148,168,75,87,50,123,116,110,104,110), dims = c(4,4), dimnames = list(c("chrM_1_G_G", "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_G"), paste0("Cell_", 1:4)))
  expected_result_T <- Matrix::sparseMatrix(i = rep(1:4,4),            j = rep(1:4,each=4),                        x = c(69,145,86,93,112,84,53,98,41,166,122,85,6,86,116,98),       dims = c(4,4), dimnames = list(c("chrM_1_G_T", "chrM_2_A_T", "chrM_3_T_T", "chrM_4_C_T"), paste0("Cell_", 1:4)))
  # We perform the tests.
  expect_equal(test_A, expected_result_A, tolerance = 1e-6)
  expect_equal(test_C, expected_result_C, tolerance = 1e-6)
  expect_equal(test_G, expected_result_G, tolerance = 1e-6)
  expect_equal(test_T, expected_result_T, tolerance = 1e-6)
})









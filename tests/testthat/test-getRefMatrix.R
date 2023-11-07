test_that("Testing getRefMatrix.R", {
  # We load the input object.
  input_object <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_Data.rds"))
  test_A <- sigurd::getRefMatrix(SE = input_object, letter = "A", chromosome_prefix = "chrM")
  test_C <- sigurd::getRefMatrix(SE = input_object, letter = "C", chromosome_prefix = "chrM")
  test_G <- sigurd::getRefMatrix(SE = input_object, letter = "G", chromosome_prefix = "chrM")
  test_T <- sigurd::getRefMatrix(SE = input_object, letter = "T", chromosome_prefix = "chrM")
  # We generate the expected results.
  expected_result_A <- c(Cell_1 =  40, Cell_2 = 50,  Cell_3 =  50, Cell_4 =  30)
  expected_result_C <- c(Cell_1 = 114, Cell_2 = 168, Cell_3 =  89, Cell_4 = 114)
  expected_result_G <- c(Cell_1 = 157, Cell_2 = 25,  Cell_3 =  75, Cell_4 = 116)
  expected_result_T <- c(Cell_1 =  86, Cell_2 = 53,  Cell_3 = 122, Cell_4 = 116)
  # We perform the tests.
  expect_equal(test_A, expected_result_A, tolerance = 1e-6)
  expect_equal(test_C, expected_result_C, tolerance = 1e-6)
  expect_equal(test_G, expected_result_G, tolerance = 1e-6)
  expect_equal(test_T, expected_result_T, tolerance = 1e-6)
})

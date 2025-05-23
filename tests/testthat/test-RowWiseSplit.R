test_that("Testing RowWiseSplit.R", {
  inputobject <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  test <- sigurd::RowWiseSplit(se = inputobject, n_cores = 1, remove_nocalls = FALSE)
  test_parallel <- sigurd::RowWiseSplit(se = inputobject, n_cores = 2, remove_nocalls = FALSE)
  test_removed_nocalls <- sigurd::RowWiseSplit(se = inputobject, n_cores = 1, remove_nocalls = TRUE)
  test_removed_nocalls_parallel <- sigurd::RowWiseSplit(se = inputobject, n_cores = 2, remove_nocalls = TRUE)
  # We generate the expected results.
  expected_result <- list("chrM_1_G_A" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_3_T_A" = c("Test_Cell_1" = 1, "Test_Cell_2" = 0, "Test_Cell_3" = 1, "Test_Cell_4" = 1),
                          "chrM_4_C_A" = c("Test_Cell_1" = 1, "Test_Cell_2" = 0, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_1_G_C" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_2_A_C" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_3_T_C" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_2_A_G" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_3_T_G" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_4_C_G" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_1_G_T" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_2_A_T" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_4_C_T" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1))
  expect_equal(test, expected_result, tolerance = 1e-6)
  expect_equal(test_parallel, expected_result, tolerance = 1e-6)
  expected_result <- list("chrM_1_G_A" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_3_T_A" = c("Test_Cell_1" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1),
                          "chrM_4_C_A" = c("Test_Cell_1" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_1_G_C" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_2_A_C" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_3_T_C" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_2_A_G" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_3_T_G" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_4_C_G" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_1_G_T" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_2_A_T" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1), 
                          "chrM_4_C_T" = c("Test_Cell_1" = 1, "Test_Cell_2" = 1, "Test_Cell_3" = 1, "Test_Cell_4" = 1))
  expect_equal(test_removed_nocalls, expected_result, tolerance = 1e-6)
  expect_equal(test_removed_nocalls_parallel, expected_result, tolerance = 1e-6)
})


















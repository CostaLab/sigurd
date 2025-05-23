test_that("Testing VariantBurden.R", {
  inputobject <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  test <- sigurd::VariantBurden(inputobject)
  # We generate the expected result.
  expected_result <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  expected_result$Burden <- c(Test_Cell_1 = 2.72123540439267, Test_Cell_2 = 3.19347598570439, Test_Cell_3 = 2.87721693179241, Test_Cell_4 = 2.85496178156784)
  expect_equal(test, expected_result, tolerance = 1e-6)
})

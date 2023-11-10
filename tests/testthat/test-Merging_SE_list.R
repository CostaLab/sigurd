test_that("Testing Merging_SE_list.R", {
  inputobject_1 <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  inputobject_2 <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  combined_object <- sigurd::Merging_SE_list(se = list(inputobject_1, inputobject_2))
  # saveRDS(combined_object, paste0(getwd(), "/test_data/Merging_SE_list_ExpectedResults.rds"))
  expected_result <- readRDS(paste0(getwd(), "/test_data/Merging_SE_list_ExpectedResults.rds"))
  expect_equal(combined_object, expected_result, tolerance = 1e-6)
})

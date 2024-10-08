test_that("Testing GetCellInfoPerVariant.R", {
  # Loading the test object.
  se <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  vois <- c("chrM_1_G_A", "chrM_3_T_C")
  test <- GetCellInfoPerVariant(se = se, voi_ch = vois, verbose = FALSE)
  # We generate the expected result.
  expected_result <- tibble::tibble(cell = paste0("Test_Cell_", 1:4), cov_chrM_1_G_A = c(279,328,230,300), af_chrM_1_G_A = c(0.07168459, 0.15243902, 0.21739130, 1/30),
                                    cov_chrM_3_T_C = c(335,264,311,331), af_chrM_3_T_C = c(0.1074627, 0.2386364, 0.2861736, 0.2749245))
  expect_equal(test, expected_result, tolerance = 1e-6)
})

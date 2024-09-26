test_that("Testing Filtering.R", {
  # Testing if blacklisting works.
  test_blacklist <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  blacklist <- read.table(paste0(getwd(), "/test_data/Filtering_Blacklisted_Barcodes.tsv"))[,1]
  test_blacklist <- sigurd::Filtering(SE = test_blacklist, cells_exclude = blacklist, min_cells_per_variant = 2, min_variants_per_cell = 1, verbose = FALSE)
  # saveRDS(test_blacklist, paste0(getwd(), "/test_data/Filtering_Blacklist_ExpectedResults.rds"))
  expected_result_blacklist <- readRDS(paste0(getwd(), "/test_data/Filtering_Blacklist_ExpectedResults.rds"))
  expect_equal(test_blacklist, expected_result_blacklist, tolerance = 1e-6)

  # Testing if fraction thresholding works.
  test_fraction_threshold <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  test_fraction_threshold <- sigurd::Filtering(SE = test_fraction_threshold, fraction_threshold = 0.25, min_cells_per_variant = 2, min_variants_per_cell = 1, verbose = FALSE)
  # saveRDS(test_fraction_threshold, paste0(getwd(), "/test_data/Filtering_Fraction_Threshold_ExpectedResults.rds"))
  expected_result_fraction_threshold <- readRDS(paste0(getwd(), "/test_data/Filtering_Fraction_Threshold_ExpectedResults.rds"))
  expect_equal(test_fraction_threshold, expected_result_fraction_threshold, tolerance = 1e-6)

  # Testing if alt reads thresholding works.
  test_alts_threshold <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  test_alts_threshold <- sigurd::Filtering(SE = test_alts_threshold, alts_threshold = 50, min_cells_per_variant = 2, min_variants_per_cell = 1, verbose = FALSE)
  # saveRDS(test_alts_threshold, paste0(getwd(), "/test_data/Filtering_Alts_Threshold_ExpectedResults.rds"))
  expected_result_alts_threshold <- readRDS(paste0(getwd(), "/test_data/Filtering_Alts_Threshold_ExpectedResults.rds"))
  expect_equal(test_alts_threshold, expected_result_alts_threshold, tolerance = 1e-6)

  # Testing if min cells per variant thresholding works.
  test_min_cells_per_variant <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  test_min_cells_per_variant <- sigurd::Filtering(SE = test_min_cells_per_variant, min_cells_per_variant = 4, min_variants_per_cell = 1, verbose = FALSE)
  # saveRDS(test_min_cells_per_variant, paste0(getwd(), "/test_data/Filtering_CellThreshold_ExpectedResults.rds"))
  expected_result_min_cells_per_variant <- readRDS(paste0(getwd(), "/test_data/Filtering_CellThreshold_ExpectedResults.rds"))
  expect_equal(test_min_cells_per_variant, expected_result_min_cells_per_variant, tolerance = 1e-6)

  # Testing if min variants per cell thresholding works.
  test_min_variants_per_cell <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  test_min_variants_per_cell <- sigurd::Filtering(SE = test_min_variants_per_cell, min_cells_per_variant = 2, min_variants_per_cell = 4, verbose = FALSE)
  # saveRDS(test_min_variants_per_cell, paste0(getwd(), "/test_data/Filtering_VariantThreshold_ExpectedResults.rds"))
  expected_result_min_variants_per_cell <- readRDS(paste0(getwd(), "/test_data/Filtering_VariantThreshold_ExpectedResults.rds"))
  expect_equal(test_min_variants_per_cell, expected_result_min_variants_per_cell, tolerance = 1e-6)
})

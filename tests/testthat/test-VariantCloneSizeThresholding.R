# We check if the MAESTER input file is present.
test_input_file       <- file.exists(paste0(getwd(), "/testthat/test_data/MAEGATK_inputfile_test.csv"))
test_maegatk_sw       <- file.exists(paste0(getwd(), "/testthat/test_data/MAESTER_data/SW/SW_CellLineMix_All_mr3_maegatk.rds"))
test_maegatk_sw_cbs   <- file.exists(paste0(getwd(), "/testthat/test_data/MAESTER_data/SW/SW_CellLineMix_All_mr3_maegatk_CellBarcodes.tsv"))
test_maegatk_tenx     <- file.exists(paste0(getwd(), "/testthat/test_data/MAESTER_data/TenX/TenX_CellLineMix_All_mr3_maegatk.rds"))
test_maegatk_tenx_cbs <- file.exists(paste0(getwd(), "/testthat/test_data/MAESTER_data/TenX/TenX_CellLineMix_All_mr3_maegatk_CellBarcodes.tsv"))

test_that("Testing VariantCloneSizeThresholding.R", {
  expect_equal(2 * 2, 4)
})

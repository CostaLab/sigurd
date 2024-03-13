test_that("Testing LoadingVCF_typewise.R", {
  # We test for a sample input file.
  test <- LoadingVCF_typewise(samples_file = paste0(getwd(), "/test_data/VCF_Test/VCF_inputfile_test.csv"), vcf_path = paste0(getwd(), "/test_data/VCF_Test/test.vcf"), patient = "Test",
                              samples_path = NULL, type_use = "scRNAseq_Somatic", min_reads = NULL, min_cells = 2, verbose = FALSE)
  # saveRDS(test, paste0(getwd(), "/test_data/VCF_Test/VCF_ExpectedResults.rds"))
  expected_result <- readRDS(paste0(getwd(), "/test_data/VCF_Test/VCF_ExpectedResults.rds"))
  expect_equal(test, expected_result, tolerance = 1e-6)

  # We test for a single sample.
  test <- LoadingVCF_typewise(samples_file = NULL, vcf_path = paste0(getwd(), "/test_data/VCF_Test/test.vcf"), patient = "Test", samples_path = paste0(getwd(), "/test_data/VCF_Test/cellSNP.cells.vcf"),
                              type_use = "scRNAseq_Somatic", min_reads = NULL, min_cells = 2, verbose = FALSE)
  expected_result <- readRDS(paste0(getwd(), "/test_data/VCF_Test/VCF_ExpectedResults.rds"))
  expect_equal(test, expected_result, tolerance = 1e-6)
})

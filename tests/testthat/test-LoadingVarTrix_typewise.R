test_that("Testing LoadingVarTrix_typewise.R", {
  # We perform the test for the a single sample.
  input_vatrix <- sigurd::LoadingVarTrix_typewise(samples_file = NULL, samples_path = paste0(getwd(), "/test_data/VarTrix_Test/"), barcodes_path = paste0(getwd(), "/test_data/VarTrix_Test/barcodes.tsv"),
                                                  snp_path = paste0(getwd(), "/test_data/VarTrix_Test/SNV.loci.txt"), vcf_path = paste0(getwd(), "/test_data/VarTrix_Test/test.vcf"), 
                                                  patient = "Test", type_use = "scRNAseq_Somatic", min_reads = NULL, min_cells = 2, verbose = FALSE)
  # saveRDS(input_vatrix, paste0(getwd(), "/test_data/VarTrix_Test/test.rds"))
  expected_results <- readRDS(paste0(getwd(), "/test_data/VarTrix_Test/test.rds"))
  expect_equal(input_vatrix, expected_results, tolerance = 1e-6)
  # We perform the test for the a sample input file.
  input_vatrix <- sigurd::LoadingVarTrix_typewise(samples_file = paste0(getwd(), "/test_data/VarTrix_Test/inputfile_test.csv"), samples_path = NULL, barcodes_path = NULL,
                                                  snp_path = NULL, vcf_path = paste0(getwd(), "/test_data/VarTrix_Test/test.vcf"), 
                                                  patient = "Test", type_use = "scRNAseq_Somatic", min_reads = NULL, min_cells = 2, verbose = FALSE)
  expected_results <- readRDS(paste0(getwd(), "/test_data/VarTrix_Test/test.rds"))
  expect_equal(input_vatrix, expected_results, tolerance = 1e-6)
})

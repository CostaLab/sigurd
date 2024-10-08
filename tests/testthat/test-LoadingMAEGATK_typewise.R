test_that("Testing LoadingMAEGATK_typewise.R", {
  # These are the reference alleles for the first 4 positions.
  ref_allele <- c("G", "A", "T", "C")
  # All possible combinations. The first allele means the reference allele at this position.
  # The second allele means the reads we are observing.
  variants <- c("chrM_1_G_A", "chrM_2_A_A", "chrM_3_T_A", "chrM_4_C_A",
                "chrM_1_G_G", "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_G",
                "chrM_1_G_T", "chrM_2_A_T", "chrM_3_T_T", "chrM_4_C_T",
                "chrM_1_G_C", "chrM_2_A_C", "chrM_3_T_C", "chrM_4_C_C")
  # We generate a matrix to have an overview over the reads we observe per cell and variant.
  reads_per_variant <- matrix(0, nrow = 16, ncol = 4, dimnames = list(variants, paste0("Cell_", 1:4)))
  reads_per_variant[,1] <- c( 0,20,0,0,  23, 0,0,0,  0,0,45, 0, 0,0,0,25)
  reads_per_variant[,2] <- c( 0, 0,0,0,  20,20,0,0, 40,0, 0, 0, 0,0,0, 0)
  reads_per_variant[,3] <- c( 0, 0,0,0,   0, 0,0,0,  0,0, 0,30, 0,0,0, 0)
  reads_per_variant[,4] <- c(20, 0,0,0,   0, 0,0,0,  0,0, 0,30, 0,0,0, 0)
  # We generate sparse matrices for the SummarizedExperimentObject.
  A_counts_fw  <- Matrix::sparseMatrix(i = c(1,2,3,4,1,1,2,3,4,1,3,4), j = c(1,1,1,1,2,3,3,3,3,4,4,4), x = c(10,20,30,40,50,10,20,30,40,10,20,30), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  A_qual_fw    <- Matrix::sparseMatrix(i = c(1,2,3,4,1,1,2,3,4,1,3,4), j = c(1,1,1,1,2,3,3,3,3,4,4,4), x = c(23,37,20,31,17,15,22,40,30,10,16,14), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  A_counts_rev <- Matrix::sparseMatrix(i = c(1,2,3,4,2,1,2,3,4,2,2,4), j = c(1,1,1,1,2,3,3,3,3,4,4,4), x = c(10,20,30,40,50,40,30,20,10,10,20,30), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  A_qual_rev   <- Matrix::sparseMatrix(i = c(1,2,3,4,1,1,2,3,4,1,3,4), j = c(1,1,1,1,2,3,3,3,3,4,4,4), x = c(17,31,23,27,37,5,24,8,32,26,33,28), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  C_counts_fw  <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(26,58,59,81,18,87,27,39,34,62,86,47,85,78,25,43), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  C_qual_fw    <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(17,36,37,38,40,26,29,35,11,32,22,19,20,21,25,10), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  C_counts_rev <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(7,83,5,87,58,54,8,9,2,1,3,44,29,90,64,71), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  C_qual_rev   <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(27,25,12,22,23,14,17,5,15,26,30,33,11,16,9,40), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  G_counts_fw  <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(90,6,17,74,35,76,61,98,81,70,33,59,12,77,27,48), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  G_qual_fw    <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(7,35,36,40,13,5,37,12,15,25,26,9,17,18,10,31), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  G_counts_rev <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(67,19,58,42,80,0,26,12,72,78,17,45,24,91,96,62), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  G_qual_rev   <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(40,39,15,5,30,21,14,33,36,27,18,28,10,11,6,35), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  T_counts_fw  <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(31, 70, 38, 4, 50, 77, 95, 41, 28, 20, 92, 65, 2, 90, 6, 17), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  T_qual_fw    <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(32,35,19,5,23,12,25,6,9,18,24,34,28,21,14,37), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  T_counts_rev <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(38, 42, 3, 2, 95, 7, 71, 45, 58, 33, 30, 51, 91, 8, 79, 81), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  T_qual_rev   <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(22,8,13,7,36,28,31,18,30,19,17,27,35,16,23,40), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))

  # We add the sparse matrices together to get the coverage.
  As <- as.matrix(A_counts_fw) + as.matrix(A_counts_rev)
  Cs <- as.matrix(C_counts_fw) + as.matrix(C_counts_rev)
  Gs <- as.matrix(G_counts_fw) + as.matrix(G_counts_rev)
  Ts <- as.matrix(T_counts_fw) + as.matrix(T_counts_rev)
  coverage <- As + Cs
  coverage <- coverage + Gs
  coverage <- coverage + Ts
  coverage <- as(coverage, "CsparseMatrix")

  # We generate a GRanges object for the SummarizedExperimentObject.
  rowRanges <- GenomicRanges::GRanges(seqnames = "chrM", ranges = IRanges::IRanges(start = 1:4, end = 1:4, width = 1), strand = "*", refAllele = ref_allele)
  # We generate the actual SummarizedExperimentObject.
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(A_counts_fw = A_counts_fw, A_counts_rev = A_counts_rev, A_qual_fw = A_qual_fw, A_qual_rev = A_qual_rev,
                                                                 C_counts_fw = C_counts_fw, C_counts_rev = C_counts_rev, C_qual_fw = C_qual_fw, C_qual_rev = C_qual_rev,
                                                                 G_counts_fw = G_counts_fw, G_counts_rev = G_counts_rev, G_qual_fw = G_qual_fw, G_qual_rev = G_qual_rev,
                                                                 T_counts_fw = T_counts_fw, T_counts_rev = T_counts_rev, T_qual_fw = T_qual_fw, T_qual_rev = T_qual_rev,
                                                                 coverage = coverage),
                                                   rowRanges = rowRanges)
  # saveRDS(se, paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_Data.rds"))
  barcodes <- data.frame(paste0("Cell_", 1:4))
  # write.table(barcodes, paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_Barcodes.tsv"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  test <- sigurd::LoadingMAEGATK_typewise(samples_file = NULL, samples_path = paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_Data.rds"), patient = "Test", type_use = "scRNAseq_MT",
                                          chromosome_prefix = "chrM", min_cells = 2, barcodes_path = paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_Barcodes.tsv"), verbose = FALSE)
  # saveRDS(test, paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  expected_result <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  # We perform the test for a single sample.
  expect_equal(test, expected_result, tolerance = 1e-6)
  # We perform the test for a sample input file.
  test <- sigurd::LoadingMAEGATK_typewise(samples_file = paste0(getwd(), "/test_data/MAEGATK_inputfile_test.csv"), samples_path = NULL, patient = "Test", type_use = "scRNAseq_MT",
                                          chromosome_prefix = "chrM", min_cells = 2, barcodes_path = paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_Barcodes.tsv"), verbose = FALSE)
  expect_equal(test, expected_result, tolerance = 1e-6)
})

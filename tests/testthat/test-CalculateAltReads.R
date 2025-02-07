test_that("Testing CalculateAltReads.R", {
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
  reads_per_variant[,1] <- c(0,20,0,0, 23,0,0,0, 0,0,45,0, 0,0,0,25)
  reads_per_variant[,2] <- c(0,0,0,0, 20,20,0,0, 40,0,0,0, 0,0,0,0)
  reads_per_variant[,3] <- c(0,0,0,0, 0,0,0,0, 0,0,0,30, 0,0,0,0)
  reads_per_variant[,4] <- c(20,0,0,0, 0,0,0,0, 0,0,0,30, 0,0,0,0)
  # We generate sparse matrices for the SummarizedExperimentObject.
  A_counts_fw  <- Matrix::sparseMatrix(i = c(2,1), j = c(1,4), x = c(20,20), dims = c(4,4),
                                       dimnames = list(NULL, paste0("Cell_", 1:4)))
  A_counts_rev <- Matrix::sparseMatrix(i = c(2,1), j = c(1,4), x = c(20,20), dims = c(4,4),
                                       dimnames = list(NULL, paste0("Cell_", 1:4)))
  C_counts_fw  <- Matrix::sparseMatrix(i = 4, j = 1, x = 25, dims = c(4,4),
                                       dimnames = list(NULL, paste0("Cell_", 1:4)))
  C_counts_rev <- Matrix::sparseMatrix(i = 4, j = 1, x = 25, dims = c(4,4),
                                       dimnames = list(NULL, paste0("Cell_", 1:4)))
  G_counts_fw  <- Matrix::sparseMatrix(i = c(1,1,2), j = c(1,2,2), x = c(23,20,20), dims = c(4,4),
                                       dimnames = list(NULL, paste0("Cell_", 1:4)))
  G_counts_rev <- Matrix::sparseMatrix(i = c(1,1,2), j = c(1,2,2), x = c(23,20,20), dims = c(4,4),
                                       dimnames = list(NULL, paste0("Cell_", 1:4)))
  T_counts_fw  <- Matrix::sparseMatrix(i = c(3,1,4,4), j = c(1,2,3,4), x = c(45,40,30,30), dims = c(4,4),
                                       dimnames = list(NULL, paste0("Cell_", 1:4)))
  T_counts_rev <- Matrix::sparseMatrix(i = c(3,1,4,4), j = c(1,2,3,4), x = c(45,40,30,30), dims = c(4,4),
                                       dimnames = list(NULL, paste0("Cell_", 1:4)))
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
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(A_counts_fw = A_counts_fw, A_counts_rev = A_counts_rev,
                                                                 C_counts_fw = C_counts_fw, C_counts_rev = C_counts_rev,
                                                                 G_counts_fw = G_counts_fw, G_counts_rev = G_counts_rev,
                                                                 T_counts_fw = T_counts_fw, T_counts_rev = T_counts_rev,
                                                                 coverage = coverage),
                                                   rowRanges = rowRanges)

  # We generate the expected results.
  expected_result <- Matrix::sparseMatrix(i = c(7,10,12,1,12), j = c(2,2,3,4,4), x = c(40,80,60,40,60), dims = c(12,4),
                                          dimnames = list(c("chrM_1_G_A", "chrM_3_T_A", "chrM_4_C_A", "chrM_1_G_C", "chrM_2_A_C", "chrM_3_T_C", "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_G", "chrM_1_G_T", "chrM_2_A_T", "chrM_4_C_T"), paste0("Cell_", 1:4)))
  expect_equal(sigurd::CalculateAltReads(SE = se, chromosome_prefix = "chrM"),  expected_result, tolerance = 1e-6)
})













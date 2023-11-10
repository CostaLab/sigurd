test_that("Testing CalculateStrandCorrelation.R", {
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
  A_counts_fw  <- Matrix::sparseMatrix(i = c(1,2,3,4,1,1,2,3,4,1,3,4), j = c(1,1,1,1,2,3,3,3,3,4,4,4), x = c(10,20,30,40,50,10,20,30,40,10,20,30), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  A_counts_rev <- Matrix::sparseMatrix(i = c(1,2,3,4,2,1,2,3,4,2,2,4), j = c(1,1,1,1,2,3,3,3,3,4,4,4), x = c(10,20,30,40,50,40,30,20,10,10,20,30), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))

  C_counts_fw  <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(26, 58, 59, 81, 18, 87, 27, 39, 34, 62, 86, 47, 85, 78, 25, 43), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  C_counts_rev <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(7, 83, 5, 87, 58, 54, 8, 9, 2, 1, 3, 44, 29, 90, 64, 71), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))

  G_counts_fw  <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(90, 6, 17, 74, 35, 76, 61, 98, 81, 70, 33, 59, 12, 77, 27, 48), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  G_counts_rev <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(67, 19, 58, 42, 80, 0, 26, 12, 72, 78, 17, 45, 24, 91, 96, 62), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))

  T_counts_fw  <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(31, 70, 38, 4, 50, 77, 95, 41, 28, 20, 92, 65, 2, 90, 6, 17), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))
  T_counts_rev <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4,4), x = c(38, 42, 3, 2, 95, 7, 71, 45, 58, 33, 30, 51, 91, 8, 79, 81), dims = c(4,4), dimnames = list(NULL, paste0("Cell_", 1:4)))

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
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(A_counts_fw = A_counts_fw, A_counts_rev = A_counts_rev,# A_qual_fw = A_qual_fw, A_qual_rev = A_qual_rev,
                                                                 C_counts_fw = C_counts_fw, C_counts_rev = C_counts_rev,# C_qual_fw = C_qual_fw, C_qual_rev = C_qual_rev,
                                                                 G_counts_fw = G_counts_fw, G_counts_rev = G_counts_rev,# G_qual_fw = G_qual_fw, G_qual_rev = G_qual_rev,
                                                                 T_counts_fw = T_counts_fw, T_counts_rev = T_counts_rev,# T_qual_fw = T_qual_fw, T_qual_rev = T_qual_rev,
                                                                 coverage = coverage),
                                                   rowRanges = rowRanges)
  # We generate the expected results.
  expected_result <- c(NA, NA, -0.1889822, 0.6921127, 0.3135761, -0.2888473, -0.8573596, 0.9442110, 0.5970276, 0.6795482, -0.1622359, -0.9911656)
  names(expected_result) <- c("chrM_1_G_A", "chrM_3_T_A", "chrM_4_C_A", "chrM_1_G_C", "chrM_2_A_C", "chrM_3_T_C", "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_G", "chrM_1_G_T",
                              "chrM_2_A_T", "chrM_4_C_T")
  expect_equal(sigurd::CalculateStrandCorrelation(SE = se, chromosome_prefix = "chrM"), expected_result, tolerance = 1e-6)
})

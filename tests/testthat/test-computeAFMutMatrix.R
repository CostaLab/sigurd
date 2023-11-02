test_that("Testing if computeAFMutMatrix works.", {
  # This is the test data.
  # A G T C
  ref_allele <- c("G", "A", "T", "C")
  variants <- c("chrM_1_G_A", "chrM_3_T_A", "chrM_4_C_A",
                "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_G",
                "chrM_1_G_T", "chrM_2_A_T", "chrM_4_C_T",
                "chrM_1_G_C", "chrM_2_A_C", "chrM_3_T_C")
  test <- matrix(0, nrow = 12, ncol = 4, dimnames = list(variants, paste0("Cell_", 1:4)))
  A_counts_fw  <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(3,3,1,1,3,2,2,1,3,1,2,2),
                                       dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  A_counts_rev <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(3,3,1,1,3,2,2,1,3,1,2,2),
                                       dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  C_counts_fw  <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(3,3,1,1,3,2,2,1,3,1,2,2),
                                       dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  C_counts_rev <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(3,3,1,1,3,2,2,1,3,1,2,2),
                                       dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  G_counts_fw  <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(3,3,1,1,3,2,2,1,3,1,2,2),
                                       dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  G_counts_rev <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(3,3,1,1,3,2,2,1,3,1,2,2),
                                       dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  T_counts_fw  <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(3,3,1,1,3,2,2,1,3,1,2,2),
                                       dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  T_counts_rev <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(3,3,1,1,3,2,2,1,3,1,2,2),
                                       dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))

  scrna_data <- rbind(consensus    = c(   0,    3,    2,   1,    3,   1,   2, 0, 0,   1,   1,    2, 0,    3,    3,    2), 
                      fraction     = c(   0, 0.54,    1,   0, 0.61,   0,   1, 0, 0,   0,   0,    1, 0, 0.95, 0.95,    1),
                      coverage     = c(   0,  100,   16,  32, 1000,  64, 128, 0, 0, 256, 512, 1024, 0, 2000, 4000, 8192),
                      alts         = c(   0,   54,   16,   0,  610,   0, 128, 0, 0,   0,   0, 1024, 0, 1900, 3800, 8192),
                      refs         = c(   0,   46,    0,  32,  390,  64,   0, 0, 0, 256, 512,    0, 0,  100,  200,    0))
  scrna_coldata <- S4Vectors::DataFrame(Cell = paste0("Cell_", 1:4), Type = "scRNAseq_Somatic", AverageCoverage = Matrix::rowMeans(scrna_coverage))
  names(scrna_coldata$AverageCoverage) <- NULL
  scrna_rowdata <- S4Vectors::DataFrame(VariantName = paste0("Variant_", 1:4), Concordance = 0:3, Depth = Matrix::colMeans(scrna_coverage))
  names(scrna_rowdata$Depth) <- NULL
  scRNAseq <- SummarizedExperiment::SummarizedExperiment(assays = list(consensus = scrna_consensus, fraction = scrna_fraction, coverage = scrna_coverage, alts = scrna_alts, refs = scrna_refs),
                                                         colData = scrna_coldata, rowData = scrna_rowdata)
  fraction <- computeAFMutMatrix(SE = se_merged, chromosome_prefix = chromosome_prefix)
})

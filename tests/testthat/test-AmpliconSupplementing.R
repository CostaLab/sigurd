test_that("Testing AmpliconSupplementing.", {
  # This is the test data.
  scrna_data <- rbind(consensus    = c(   0,    3,    2,   1,    3,   1,   2, 0, 0,   1,   1,    2, 0,    3,    3,    2), 
                      fraction     = c(   0, 0.54,    1,   0, 0.61,   0,   1, 0, 0,   0,   0,    1, 0, 0.95, 0.95,    1),
                      coverage     = c(   0,  100,   16,  32, 1000,  64, 128, 0, 0, 256, 512, 1024, 0, 2000, 4000, 8192),
                      alts         = c(   0,   54,   16,   0,  610,   0, 128, 0, 0,   0,   0, 1024, 0, 1900, 3800, 8192),
                      refs         = c(   0,   46,    0,  32,  390,  64,   0, 0, 0, 256, 512,    0, 0,  100,  200,    0))
  amplicon_data <- rbind(consensus = c(   2,    3,    3,   1,    3,   1,   2, 0, 0,   1,   1,    2, 0,    3,    3,    2), 
                         fraction  = c(   1, 0.54,  0.8,   0, 0.91,   0,   1, 0, 0,   0,   0,    1, 0, 0.95, 0.95,    1),
                         coverage  = c(1000,  100,  500,  32, 1000, 500, 128, 0, 0, 256, 512, 1024, 0, 3000, 8000, 8192),
                         alts      = c(1000,   54,  400,   0,  910,   0, 128, 0, 0,   0,   0, 1024, 0, 2850, 7600, 8192),
                         refs      = c(   0,   46,  100,  32,   90, 500,   0, 0, 0, 256, 512,    0, 0,  150,  400,    0))
  # Generating sparse matrices to generate a SummarizedExperiment object.
  scrna_consensus <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(3,3,1,1,3,2,2,1,3,1,2,2),
                                          dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  scrna_fraction  <- Matrix::sparseMatrix(i = c(2,1,4,1,2,4,3,4), j = c(1,2,2,3,3,3,4,4), x = c(0.61, 0.54, 0.95, 1, 1, 0.95, 1, 1),
                                          dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  scrna_coverage  <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,2,2,2,2,3,3,3,3,4,4,4), x = c(1000,100,64,256,2000,16,128,512,4000,32,1024,8192),
                                          dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  scrna_alts      <- Matrix::sparseMatrix(i = c(2,1,4,1,2,4,3,4), j = c(1,2,2,3,3,3,4,4), x = c(610,54,1900,16,128,3800,1024,8192),
                                          dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  scrna_refs      <- Matrix::sparseMatrix(i = c(2,1,2,3,4,3,4,1), j = c(1,2,2,2,2,3,3,4), x = c(390,46,64,256,100,512,200,32),
                                          dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  scrna_coldata <- S4Vectors::DataFrame(Cell = paste0("Cell_", 1:4), Patient = "Patient1", Sample = "Sample1", Type = "scRNAseq_Somatic", AverageCoverage = Matrix::rowMeans(scrna_coverage))
  names(scrna_coldata$AverageCoverage) <- NULL
  scrna_rowdata <- S4Vectors::DataFrame(VariantName = paste0("Variant_", 1:4), Concordance = 0:3, Depth = Matrix::colMeans(scrna_coverage))
  names(scrna_rowdata$Depth) <- NULL
  scRNAseq <- SummarizedExperiment::SummarizedExperiment(assays = list(consensus = scrna_consensus, fraction = scrna_fraction, coverage = scrna_coverage, alts = scrna_alts, refs = scrna_refs),
                                                         colData = scrna_coldata, rowData = scrna_rowdata)

  amplicon_consensus <- Matrix::sparseMatrix(i = c(1,2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,1,2,2,2,2,3,3,3,3,4,4,4), x = c(2,3,3,1,1,3,3,2,1,3,1,2,2),
                                             dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  amplicon_fraction  <- Matrix::sparseMatrix(i = c(1,2,1,4,1,2,4,3,4), j = c(1,1,2,2,3,3,3,4,4), x = c(1,0.91,0.54,0.95,0.8,1,0.95,1,1),
                                             dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  amplicon_coverage  <- Matrix::sparseMatrix(i = rep(1:4, each = 4), j = rep(1:4, 4), x = amplicon_data["coverage",],
                                             dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  amplicon_coverage  <- Matrix::sparseMatrix(i = c(1,2,1,2,3,4,1,2,3,4,1,3,4), j = c(1,1,2,2,2,2,3,3,3,3,4,4,4), x = c(1000,1000,100,500,256,3000,500,128,512,8000,32,1024,8192),
                                             dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  amplicon_alts      <- Matrix::sparseMatrix(i = c(1,2,1,4,1,2,4,3,4), j = c(1,1,2,2,3,3,3,4,4), x = c(1000,910,54,2850,400,128,7600,1024,8192),
                                             dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))
  amplicon_refs      <- Matrix::sparseMatrix(i = c(2,1,2,3,4,1,3,4,1), j = c(1,2,2,2,2,3,3,3,4), x = c(90,46,500,256,150,100,512,400,32),
                                             dimnames = list(paste0("Variant_", 1:4), paste0("Cell_", 1:4)))

  amplicon_coldata <- S4Vectors::DataFrame(Cell = paste0("Cell_", 1:4), Patient = "Patient1", Sample = "Sample1", Type = "scRNAseq_Amplicon", AverageCoverage = Matrix::rowMeans(amplicon_coverage))
  names(amplicon_coldata$AverageCoverage) <- NULL
  amplicon_rowdata <- S4Vectors::DataFrame(VariantName = paste0("Variant_", 1:4), Concordance = 3:0, VariantQuality = 1:4, Depth = Matrix::colMeans(amplicon_coverage))
  names(amplicon_rowdata$Depth) <- NULL
  amplicon <- SummarizedExperiment::SummarizedExperiment(assays = list(consensus = amplicon_consensus, fraction = amplicon_fraction, coverage = amplicon_coverage, alts = amplicon_alts, refs = amplicon_refs),
                                                         colData = amplicon_coldata, rowData = amplicon_rowdata)
  # The objects for the test.
  test_result <- sigurd::AmpliconSupplementing(scRNAseq, amplicon, verbose = FALSE)
  test_coldata <- S4Vectors::DataFrame(Cell = scrna_coldata$Cell, Patient = "Patient1", Sample = "Sample1", SamplescRNAseq = "Sample1", TypescRNAseq = scrna_coldata$Type, AverageCoveragescRNAseq = scrna_coldata$AverageCoverage,
                                       TypeAmplicon = amplicon_coldata$Type, AverageCoverageAmplicon = amplicon_coldata$AverageCoverage, AverageCoverage = amplicon_coldata$AverageCoverage, 
                                       row.names = scrna_coldata$Cell)
  test_rowdata <- S4Vectors::DataFrame(VariantName = paste0("Variant_", 1:4), ConcordancescRNAseq = 0:3, DepthscRNAseq = c(250, 605, 1164, 2312), ConcordanceAmplicon = 3:0, VariantQuality = 1:4, DepthAmplicon = c(500, 964, 2285, 2312), Concordance = 3:0, Depth = c(500, 964, 2285, 2312),
                                       row.names = paste0("Variant_", 1:4))

  testthat::expect_equal(names(assays(test_result)),  c("consensus", "fraction", "coverage", "alts", "refs"))
  testthat::expect_equal(colData(test_result), test_coldata)
  testthat::expect_equal(rowData(test_result), test_rowdata)
  testthat::expect_equal(assays(test_result)[["consensus"]], amplicon_consensus)
  testthat::expect_equal(assays(test_result)[["fraction"]], amplicon_fraction)
  testthat::expect_equal(assays(test_result)[["coverage"]], amplicon_coverage)
  testthat::expect_equal(assays(test_result)[["alts"]], amplicon_alts)
  testthat::expect_equal(assays(test_result)[["refs"]], amplicon_refs)
})

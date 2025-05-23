test_that("Testing SetVariantInfo.R", {
  inputobject <- readRDS(paste0(getwd(), "/test_data/LoadingMAEGATK_typewise_Test_ExpectedResults.rds"))
  counts_matrix <- Matrix::sparseMatrix(i = rep(1:4,6), j = rep(1:6,each=4), x = c(118,112,85,85,98,96,106,86,103,80,87,91,75,120,78,110,98,87,109,82,122,114,91,111),
                                        dims = c(4,6), dimnames = list(paste0("Gene_", 1:4), paste0("Test_Cell_", 1:6)))
  seurat_object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts_matrix, assay = "RNA"))
  test <- SetVariantInfo(SE = inputobject, seurat_object = seurat_object, information = "consensus", variants = "chrM_3_T_A")
  # We generate the expected result.
  expected_result <- seurat_object
  expected_result <- Seurat::AddMetaData(object = seurat_object, metadata = c("Alt", "NoCall", "Alt", "Alt", NA, NA), col.name = "chrM_3_T_A_consensus")
  expect_equal(test, expected_result, tolerance = 1e-6)
})

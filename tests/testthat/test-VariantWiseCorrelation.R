test_that("Testing VariantWiseCorrelation.R", {
  # We generate the input lists.
  variants_list <- list(JAK2_V617F         = c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1), 
                        chr11_61796992_G_C = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0),
                        chrM_1_G_A         = c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1),
                        chrM_2_A_G         = c(0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0),
                        chrM_3_T_G         = c(0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0),
                        chrM_4_C_T         = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1))
  for(i in 1:length(variants_list)) names(variants_list[[i]]) <- paste0("Cell_", 1:20)
  # We generate the expected output.
  expected_result1 <- data.frame(Variant1 = c(rep("JAK2_V617F", 3), rep("chr11_61796992_G_C", 3), rep(c("chrM_1_G_A", "chrM_4_C_T"), each = 3)),
                                 Variant2 = c("chr11_61796992_G_C", "chrM_1_G_A", "chrM_4_C_T", "JAK2_V617F", "chrM_1_G_A", "chrM_4_C_T", "JAK2_V617F", "chr11_61796992_G_C", "chrM_4_C_T", "JAK2_V617F", "chr11_61796992_G_C", "chrM_1_G_A"),
                                 P = c(0.4175949, 0, 0.4629459, 0.4175949, 0.4175949, 0.4629459, 0, 0.4175949, 0.4629459, 0.4629459, 0.4629459, 0.4629459),
                                 Corr = c(0.1919192, 1, 0.1740777, 0.1919192, 0.1919192, 0.1740777, 1, 0.1919192, 0.1740777, 0.1740777, 0.1740777, 0.1740777),
                                 Cells_1_Alt = c(rep(11,9), rep(15,3)), Cells_1_Ref = c(rep(9,9), rep(5,3)),
                                 Cells_2_Alt = c(rep(c(11,11,15), 3), rep(11, 3)), Cells_2_Ref = c(rep(c(9,9,5),3), rep(9,3)),
                                 P_adj = c(0.4629459, 0, 0.4629459, 0.4629459, 0.4629459, 0.4629459, 0, 0.4629459, 0.4629459, 0.4629459, 0.4629459, 0.4629459))
  expected_result2 <- data.frame(Variant1 = c(rep("JAK2_V617F", 3), rep("chr11_61796992_G_C", 3), rep(c("chrM_1_G_A", "chrM_4_C_T"), each = 3)),
                                 Variant2 = c("chr11_61796992_G_C", "chrM_1_G_A", "chrM_4_C_T", "JAK2_V617F", "chrM_1_G_A", "chrM_4_C_T", "JAK2_V617F", "chr11_61796992_G_C", "chrM_4_C_T", "JAK2_V617F", "chr11_61796992_G_C", "chrM_1_G_A"),
                                 P = c(0.4175949, 0, 0.4629459, 0.4175949, 0.4175949, 0.4629459, 0, 0.4175949, 0.4629459, 0.4629459, 0.4629459, 0.4629459),
                                 Corr = c(0.1919192, 1, 0.1740777, 0.1919192, 0.1919192, 0.1740777, 1, 0.1919192, 0.1740777, 0.1740777, 0.1740777, 0.1740777),
                                 Cells_1_Alt = c(rep(11,9), rep(15,3)), Cells_1_Ref = c(rep(9,9), rep(5,3)),
                                 Cells_2_Alt = c(rep(c(11,11,15), 3), rep(11, 3)), Cells_2_Ref = c(rep(c(9,9,5),3), rep(9,3)),
                                 P_adj = c(1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1))
  expect_equal(sigurd::VariantWiseCorrelation(variants_list, n_cores = 1, p_value_adjustment = "fdr", verbose = FALSE), expected_result1, tolerance = 1e-6)
  expect_equal(sigurd::VariantWiseCorrelation(variants_list, n_cores = 2, p_value_adjustment = "fdr", verbose = FALSE), expected_result1, tolerance = 1e-6)
  expect_equal(sigurd::VariantWiseCorrelation(variants_list, n_cores = 1, p_value_adjustment = "bonferroni", verbose = FALSE), expected_result2, tolerance = 1e-6)
})

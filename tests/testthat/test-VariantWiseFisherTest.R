test_that("Testing VariantWiseFisherTest.R", {
  # We generate the input lists.
  variants_list <- list(JAK2_V617F         = c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1), 
                        chr11_61796992_G_C = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0),
                        chrM_1_G_A         = c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1),
                        chrM_2_A_G         = c(0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0),
                        chrM_3_T_G         = c(0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0),
                        chrM_4_C_T         = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1))
  for(i in 1:length(variants_list)) names(variants_list[[i]]) <- paste0("Cell_", 1:20)
  # We generate the expected output.
  test1 <- sigurd::VariantWiseFisherTest(variants_list, n_cores = 1, p_value_adjustment = "fdr", verbose = FALSE)
  test2 <- sigurd::VariantWiseFisherTest(variants_list, n_cores = 2, p_value_adjustment = "fdr", verbose = FALSE)
  test3 <- sigurd::VariantWiseFisherTest(variants_list, n_cores = 1, p_value_adjustment = "bonferroni", verbose = FALSE)
  expected_result1 <- data.frame(Variant1 = c(rep("JAK2_V617F", 5), rep("chr11_61796992_G_C", 5)),
                                 Variant2 = c("chr11_61796992_G_C", "chrM_1_G_A", "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_T", "JAK2_V617F", "chrM_1_G_A", "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_T"),
                                 P = c(6.534175e-01, 5.953799e-06, 4.058585e-01, 1, 6.168731e-01, 6.534175e-01, 6.534175e-01, 1, 4.058585e-01, 6.168731e-01),
                                 OddsRatio = c(2.1011776, Inf, 0.4357796, 0.9619440, 2.1582974, 2.1011776, 2.1011776, 0.9619440, 0.4357796, 2.1582974),
                                 Cells_Alt_1_2 = c(7,11,5,6,9,7,7,6,5,9), Cells_Alt_1_Ref_2 = c(4,0,6,5,2,4,4,5,6,2), Cells_Alt_2_Ref_1 = c(4,0,6,5,6,4,4,5,6,6),
                                 Cells_Ref_1_2 = c(5,9,3,4,3,5,5,4,3,3), P_adj = c(8.167719e-01, 5.953799e-05, 8.167719e-01, 1, 8.167719e-01, 8.167719e-01, 8.167719e-01, 1, 8.167719e-01, 8.167719e-01))
  expected_result2 <- data.frame(Variant1 = c(rep("JAK2_V617F", 5), rep("chr11_61796992_G_C", 5)),
                                 Variant2 = c("chr11_61796992_G_C", "chrM_1_G_A", "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_T", "JAK2_V617F", "chrM_1_G_A", "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_T"),
                                 P = c(6.534175e-01, 5.953799e-06, 4.058585e-01, 1, 6.168731e-01, 6.534175e-01, 6.534175e-01, 1, 4.058585e-01, 6.168731e-01),
                                 OddsRatio = c(2.1011776, Inf, 0.4357796, 0.9619440, 2.1582974, 2.1011776, 2.1011776, 0.9619440, 0.4357796, 2.1582974),
                                 Cells_Alt_1_2 = c(7,11,5,6,9,7,7,6,5,9), Cells_Alt_1_Ref_2 = c(4,0,6,5,2,4,4,5,6,2), Cells_Alt_2_Ref_1 = c(4,0,6,5,6,4,4,5,6,6),
                                 Cells_Ref_1_2 = c(5,9,3,4,3,5,5,4,3,3), P_adj = c(1, 5.953799e-05, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(sigurd::VariantWiseFisherTest(variants_list, n_cores = 1, p_value_adjustment = "fdr", verbose = FALSE),  expected_result1, tolerance = 1e-6)
  expect_equal(sigurd::VariantWiseFisherTest(variants_list, n_cores = 2, p_value_adjustment = "fdr", verbose = FALSE),  expected_result1, tolerance = 1e-6)
  expect_equal(sigurd::VariantWiseFisherTest(variants_list, n_cores = 1, p_value_adjustment = "bonferroni", verbose = FALSE),  expected_result2, tolerance = 1e-6)
})

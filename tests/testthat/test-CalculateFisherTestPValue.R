test_that("Test CalculateFisherTestPValue.R", {
  # We generate the input lists.
  variants_list <- list(JAK2_V617F         = c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1), 
                        chr11_61796992_G_C = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0),
                        chrM_1_G_A         = c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1),
                        chrM_2_A_G         = c(0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0),
                        chrM_3_T_G         = c(0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0),
                        chrM_4_C_T         = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1))
  for(i in 1:length(variants_list)) names(variants_list[[i]]) <- paste0("Cell_", 1:20)
  
  # We generate the expected output.
  expected_result1 <- c(0.6534175, 2.1011776, 7, 4, 4, 5)
  names(expected_result1) <- c("", "odds ratio", "", "", "", "")
  expect_equal(sigurd::CalculateFisherTestPValue(variant_values = variants_list[["JAK2_V617F"]], other_mutation = "chr11_61796992_G_C", all_variants_list = variants_list[names(variants_list) != "JAK2_V617F"]),  expected_result1, tolerance = 1e-6)
})

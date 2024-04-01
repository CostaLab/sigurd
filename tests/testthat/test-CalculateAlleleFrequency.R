test_that("Testing CalculateAlleleFrequency.R", {
  # We generate the input matrix.
  variants <- c("chrM_1_G_A", "chrM_3_T_A", "chrM_4_C_A",
                "chrM_2_A_G", "chrM_3_T_G", "chrM_4_C_G",
                "chrM_1_G_T", "chrM_2_A_T", "chrM_4_C_T",
                "chrM_1_G_C", "chrM_2_A_C", "chrM_3_T_C")
  sample(seq(from = 0, to = 100, by = 10), size = 48, replace = TRUE)
  ref_reads <- matrix(c( 0, 90, 80, 100,  90,  80,   0, 100, 80,   0, 100, 90,
                         0, 50, 40,  50,  50,  40,   0,  50, 40,   0,  50, 50,
                        80, 50, 90,  40,  50,  90,  80,  40, 90,  80,  40, 50,
                        70, 20, 80,   0,  20,  80,  70,   0, 80,  70,   0, 20), 
                      ncol = 4, nrow = 12, dimnames = list(variants, paste0("Cell_", 1:4)))
  alt_reads <- matrix(c(90, 10, 10,  50,   0,   0,   0,   0,  0,   0,   0,  0,
                        90, 50,  0,   0, 100,   0,   0,   0,  0,   0,  40, 60,
                         0,  0,  0,   0,   0,   0,   0,   0, 90,  10,   0, 80,
                        60, 70, 30,  10,  80, 100, 100,  20,  0, 100,   0,  0),
                      ncol = 4, nrow = 12, dimnames = list(variants, paste0("Cell_", 1:4)))
  
  expected_result1 <- matrix(c(1, 0.1, 1/9, 1/3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.5, 0, 0, 2/3, 0, 0, 0, 0, 0, 0.4444444, 0.5454545, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1/9, 0, 0.6153846, 0.4615385, 0.7777778, 0.2727273, 1, 0.8, 0.5555556,
                               0.5882353, 1, 0, 0.5882353, 0, 0),
                             ncol = 4, nrow = 12, dimnames = list(variants, paste0("Cell_", 1:4)))
  expected_result2 <- matrix(c(1, 0.1, 1/9, 1/3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.5, 0, 0, 2/3, 0, 0, 0, 0, 0, 0.4444444, 0.5454545, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1/9, 0, 0.6153846, 0.4615385, 0.7777778, 0.2727273, 0.9999999, 0.8,
                               0.5555556, 0.5882353, 1, 0, 0.5882353, 0, 0),
                             ncol = 4, nrow = 12, dimnames = list(variants, paste0("Cell_", 1:4)))
  expected_result3 <- matrix(c(0.989011, 0.0990099, 0.1098901, 0.3311258, 0, 0, 0, 0, 0, 0, 0, 0, 0.989011, 0.4950495, 0, 0, 0.6622517, 0, 0, 0, 0, 0, 0.4395604, 0.5405405, 0, 0, 0, 0, 0, 0, 0, 0, 0.4972376, 0.1098901,
                               0, 0.6106870, 0.4580153, 0.7692308, 0.2702703, 0.9090909, 0.7920792, 0.5524862, 0.5847953, 0.952381, 0, 0.5847953, 0, 0),
                             ncol = 4, nrow = 12, dimnames = list(variants, paste0("Cell_", 1:4)))
  testthat::expect_equal(sigurd::CalculateAlleleFrequency(reference_reads = ref_reads, alternative_reads = alt_reads, pseudo_count = 0),  expected_result1, tolerance = 1e-6)
  testthat::expect_equal(sigurd::CalculateAlleleFrequency(reference_reads = ref_reads, alternative_reads = alt_reads, pseudo_count = 0.000001),  expected_result2, tolerance = 1e-6)
  testthat::expect_equal(sigurd::CalculateAlleleFrequency(reference_reads = ref_reads, alternative_reads = alt_reads, pseudo_count = 1),  expected_result3, tolerance = 1e-6)
})

test_that("multiplication works", {
  x <- c("JAK2_p.V617I_c.1849G>A", "ABL1_p.E274K_c.820G>A")
  y <- c("chrM_1_G_A", "chrM_3_T_A")
  expected_result <- c(x, y)
  expect_equal(sigurd::combine_NAMES(x = x, y = y), expected_result)
})

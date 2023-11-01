testthat::test_that("ggsci_pal returns a valid color palette for aaas", {
  colors <- sigurd::ggsci_pal("aaas")(10)
  colors_expected <- c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", "#808180FF", "#1B1919FF")
  testthat::expect_true(all(colors == colors_expected), "Testing if the colors of the aaas palette are as expected.")
})

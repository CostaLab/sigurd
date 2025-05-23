test_that("Checking if genotyping conversion to numeric works.", {
  expect_equal(sigurd::char_to_numeric("1/1"),  2)
  expect_equal(sigurd::char_to_numeric("1/0"),  2)
  expect_equal(sigurd::char_to_numeric("0/1"),  2)
  expect_equal(sigurd::char_to_numeric("0/0"),  1)
  expect_equal(sigurd::char_to_numeric("asdf"), 0)
})

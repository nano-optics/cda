context("Checking the behavior of small utility functions")

test_that("Euler rotation returns a matrix", {
  m <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  m2 <- rbind(c(0, 1, 0), c(-1, 0, 0), c(0, 0, 1))
  expect_equal(cda$euler(0, 0, 0), m)
  expect_equal(cda$euler(pi/2, 0, 0), m2)
})

context("Alpha diversity")
library(mbm)

x <- matrix(c(
	0,1,0,3,10,
	0,0,0,15,0,
	2,7,3,3,1), byrow=TRUE, nrow=3)

test_that("Species richness", {
	expect_error(richness(-1 * x))
	expect_identical(richness(x), rowSums(x > 0))
})

test_that("Simpson", {
	expect_error(simpson(x))
	expect_equal(simpson(x, proportion = FALSE), c(0.56, 1, 0.28), tol=0.005)
	expect_equal(simpson(x, proportion = FALSE), simpson(x/rowSums(x), proportion = TRUE))
})

test_that("Shannon", {
	expect_warning(shannon(x))
	expect_equal(shannon(x/rowSums(x)), c(0.76, 0, 1.42), tol=0.005)
})

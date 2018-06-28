context("MBM Models")
library(mbm)

# make symmetric y matrix and some random covariates
y <- matrix(rbeta(100, 1, 1), ncol=10); diag(y) <- 0; y[upper.tri(y)] <- t(y)[upper.tri(y)]
x <- matrix(rnorm(30), ncol=3)
test_that("basic model produces no errors", {
	expect_error(mbm(y,x), NA)
	})
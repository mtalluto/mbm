context("MBM Models")
library(mbm)

# make symmetric y matrix and some random covariates
y <- matrix(rbeta(100, 1, 1), ncol=10); diag(y) <- 0; y[upper.tri(y)] <- t(y)[upper.tri(y)]
x <- matrix(rnorm(30), ncol=3)

test_that("basic model produces no errors", {
	expect_error(mbm(y,x), NA)
	})

# test link functions
test_that("adding link functions succeed", {
	expect_match(mbm(y,x)$pyclasses$link, "Identity", all=FALSE)
	expect_match(mbm(y,x,link='probit')$pyclasses$link, "Probit", all=FALSE)
	})
# test lengthscale
# test scaling
# test samples
# test (or drop) response curves
# test setting lengthscale
# test mean function
# test svgp

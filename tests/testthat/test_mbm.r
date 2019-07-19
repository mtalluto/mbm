context("MBM Models")
library(mbm)

# make symmetric y matrix and some random covariates
y <- matrix(rbeta(100, 1, 1), ncol=10); diag(y) <- 0; y[upper.tri(y)] <- t(y)[upper.tri(y)]
x <- matrix(rnorm(30), ncol=3)

test_that("basic model options produce no errors", {
	expect_error(mod <- mbm(y,x), NA)
	# expect_error(mod <- mbm(y,x, sparse = TRUE), NA)
	expect_error(mod <- mbm(y,x, force_increasing = TRUE), NA)
	expect_error(mod <- mbm(y,x, sparse = TRUE), NA)
	})

# test link functions
test_that("adding link functions succeeds", {
	link = 'identity'
	mod <- mbm(y,x, link = link)
	expect_match(mod$link, link)
	expect_match(tolower(class(mod$pyobj$likelihood$gp_link)), link, all=FALSE)
	
	link = 'probit'
	mod <- mbm(y,x, link = link)
	expect_match(mod$link, link)
	expect_match(tolower(class(mod$pyobj$likelihood$gp_link)), link, all=FALSE)
})

# test prediction
test_that("mbm prediction", {
	newx <- cbind(seq(-2, 2, length.out=20), seq(-1, 1, length.out=20), 
		seq(-0.5, 0.5, length.out=20))
	mod <- mbm(y,x)
	expect_error(pr <- predict(mod), regex=NA)
	expect_equal(nrow(pr), nrow(mod$covariates))
	expect_error(pr <- predict(mod, newdata = newx), regex=NA)
	expect_equal(nrow(pr), length(dist(newx)))
})
# test lengthscale
# test scaling
# test samples
# test (or drop) response curves
# test setting lengthscale
# test mean function
# test svgp
# test methods (printing etc)
# test spatial predict
## ESSENTIAL - test that y transformations work as expected
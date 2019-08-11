context("MBM Models")
library(mbm)

# make symmetric y matrix and some random covariates
set.seed(100)
y <- matrix(rbeta(100, 1, 1), ncol=10); diag(y) <- 0; y[upper.tri(y)] <- t(y)[upper.tri(y)]
x <- matrix(rnorm(30), ncol=3)

test_that("Python loads", {
	expect_error(check_python(), regex = NA)
})


test_that("Common errors", {
	expect_error(mod <- mbm(y,x, link = "log"))
	expect_error(mod <- mbm(y,x, likelihood = "poisson"))

	yy <- matrix(rbeta(200^2, 1, 1), ncol=200); diag(y) <- 0; y[upper.tri(y)] <- t(y)[upper.tri(y)]
	xx <- matrix(rnorm(200*3), ncol=3)
	expect_error(mbm(yy,xx), regex = 'memory')
	mod <- mbm(y,x)
	expect_error(predict(mod, newdata=x[,-1]), "variables")
	expect_warning(confint(mod, method="sample"), "not implemented")
	expect_error(confint(mod, parm="newdata"), "supported")
})



test_that("basic model options and standard methods produce no errors", {
	expect_error(mod <- mbm(y,x), NA)
	expect_error(gp_params(mod), regex = 'sparse')
	expect_error(mod <- mbm(y,x, force_increasing = TRUE), NA)

	## svgps with small numbers of iterations always throw a warning
	expect_warning(mod <- mbm(y,x, sparse = TRUE, sparse_iter = 10), regex = "iterations")
	expect_warning(mod <- mbm(y,x, sparse = TRUE, sparse_iter = 10, link = 'probit'),
		regex = "iterations")
	expect_warning(mod <- mbm(y,x, sparse = TRUE, sparse_iter = 10, force_increasing = TRUE), 
		regex = "iterations")
	# expect_error({sink("/dev/null"); summary(mod); sink()}, regex=NA)
	expect_equal(dim(inducing(mod)), dim(mod$inducing_inputs))
	chol_pars <- gp_params(mod)
	expect_equal(length(chol_pars$mean), attr(mod, "inducing_inputs"))
	expect_equal(dim(chol_pars$cholesky), rep(attr(mod, "inducing_inputs"), 2))
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
	mod <- mbm(y,x, link = 'probit')
	prResp <- predict(mod, type='response')[,1]
	expect_true(all(prResp >= 0 & prResp <= 1))
})

test_that("methods", {
	mod <- mbm(y,x)
	expect_true(is.mbm(mod))
	vdiffr::expect_doppelganger('print-method', mod)
	vdiffr::expect_doppelganger('plot-noline-noserr', plot(mod, line = FALSE, sterr=FALSE))
	vdiffr::expect_doppelganger('plot-noline-serr', plot(mod, line = FALSE, sterr=TRUE))
	vdiffr::expect_doppelganger('plot-line-noserr', plot(mod, line = TRUE, sterr=FALSE))
	vdiffr::expect_doppelganger('plot-line-serr', plot(mod, line = TRUE, sterr=TRUE))
	expect_error(ci <- confint(mod), NA)
	expect_true(all(is.finite(ci)))

})
# test samples
# test (or drop) response curves
# test setting lengthscale
# test methods (printing etc)
# test spatial predict

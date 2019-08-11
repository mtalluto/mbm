context("Beta diversity")
library(mbm)

nsites <- 10
nspecies <- 5
pr <- 0.3
yy <- matrix(sample(c(0,1), nspecies*nsites, replace = TRUE, prob = c(1-pr,pr)), 
	nrow = nsites, ncol=nspecies)
y <- matrix(1, nrow = nsites, ncol = nspecies)

test_that("Beta diversity works as expected", {
	expect_error(jc <- jaccard(y), NA)
	expect_true(all(jc == 0))
	jc <- jaccard(yy)
	expect_true(all(jc >= 0 & jc <= 1 | is.nan(jc)))

	expect_error(bray <- bc(y), NA)
	expect_true(all(bray == 0))
	bray <- bc(yy)
	expect_true(all(bray >= 0 & bray <= 1 | is.nan(bray)))

	expect_error(sord <- sorensen(y), NA)
	expect_true(all(sord == 0))
	sord <- sorensen(yy)
	expect_true(all(sord >= 0 & sord <= 1 | is.nan(sord)))
	expect_lt(sqrt(mean((sord - mbm:::sor(yy))^2, na.rm=TRUE)), 1e-4)
	expect_error(sord <- sorensen(yy, TRUE), NA)
	expect_true(all(sapply(sord, function(x) x >= 0 & x <= 1 | is.nan(x))))
})
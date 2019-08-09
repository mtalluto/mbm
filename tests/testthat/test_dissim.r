context("Beta diversity")
library(mbm)

nsites <- 10
nspecies <- 5
pr <- 0.3
yy <- matrix(sample(c(0,1), nspecies*nsites, replace = TRUE, prob = c(1-pr,pr)), 
	nrow = nsites, ncol=nspecies)

test_that("Beta diversity works as expected", {
	y <- matrix(1, nrow = nsites, ncol = nspecies)
	expect_error(jc <- jaccard(y), NA)
	expect_true(all(jc == 0))
	jc <- jaccard(yy)
	expect_true(all(jc >= 0 & jc <= 1 | is.nan(jc)))
})
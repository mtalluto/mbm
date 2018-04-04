set.seed(419124096)

# params have been chosen to fill a landscape going from -5 to 5 if 300 spp are used
# for covariance, best to use large params for rho; 100,100 makes for just a bit
# if the rho params are equal, mean correlation will be 0; the larger the numbers the
# smaller the variance
make_niche <- function(center_means = c(0,0), center_sig = 12*diag(2), 
	width_shapes = c(4,7), width_rates = c(1.2, 1.2), rho_ab = c(0,0),
	scale_ab = c(15,5))
# for width parameters, first entry in each describes width for env1, second
# for env2
# for rho params, rho ~ 2 * Beta(rho_ab[1], rho_ab[2]) - 1
# for the special case of rho = 0,0, covariance will be zero
{
	require(mvtnorm)
	# vector of two, the niche in each dimension
	cntrs <- rmvnorm(1, center_means, center_sig)
	# niche widths
	widths <- rgamma(2, width_shapes, width_rates)
	# niche axis correlations
	if(all(rho_ab == 0)) {
		rho <- 0
	} else {
		rho <- 2 * rbeta(1, rho_ab[1], rho_ab[2]) - 1
	}
	# how high is the niche; controls abundance and detectability of the species
	scale <- rbeta(1, scale_ab[1], scale_ab[2])

	data.frame(center1 = cntrs[1], center2 = cntrs[2], width1 = widths[1], 
		width2=widths[2], rho=rho, scale=scale)
}

simulate_landscape <- function(nsp, xrange=c(-5,5), yrange=c(-5,5), dims=c(50,50), ...)
{
	niches <- do.call(rbind, lapply(1:nsp, function(x) make_niche(...)))
	lscape <- expand.grid(env1=seq(xrange[1], xrange[2], length.out=dims[1]), 
		env2=seq(yrange[1], yrange[2], length.out=dims[2]))
	sitenames <- paste0('site', 1:nrow(lscape))
	spnames <- paste0('species', 1:nsp)
	rownames(lscape) <- sitenames
	probs <- apply(niches, 1, function(ni) niche_height(lscape, ni))
	colnames(probs) <- spnames
	site_sp <- matrix(rbinom(length(probs), 1, probs), ncol=ncol(probs), 
		dimnames = dimnames(probs))
	res <- list(site_sp = site_sp, niches=niches, landscape = lscape, probs = probs)
	class(res) <- c(class(res), "mbmSim")
	res
}

niche_height <- function(x, niche)
{
	require(mvtnorm)
	# computes the niche height as the multivariate normal density times the scale
	mu <- niche[1:2]
	sig <- matrix(NA, nrow=2, ncol=2)
	sig[1,1] <- niche[3]
	sig[2,2] <- niche[4]
	sig[1,2] <- sig[2,1] <- niche[3]*niche[4]*niche[5]
	sc <- niche[6]
	# scale height to one
	ht <- dmvnorm(x, mu, sig) / dmvnorm(mu, mu, sig)
	sc * ht
}



mbm_sims <- simulate_landscape(300, dims=c(100,100))
use_data(mbm_sims, overwrite=TRUE)

#' Construct (but do not fit) an mbm model
#' 
#' @param x Matrix giving a series of covariates (in columns) for all sites (in rows). 
#' 			Row names are required. All variables will be included in the model.
#' @param y Square dissimilarity or distance matrix, can be complete or lower triangular 
#' 			only. Row and column names are required and must match the site names in the 
#' 			rows of \code{x}.
#' @param y_name A name to give to the y variable
#' @param link Link function to use
#' @param likelihood Likelihood function to use
#' @param lengthscale Either missing (in which case all lengthscales will be optimized) or
#' 			a numeric vector of length \code{ncol(x)+1}. If a vector, the first entry 
#' 			corresponds to environmental distance, and entries \code{i = 1 + (1:n)} 
#' 			correspond to the variable in x[,i]. Values must be \code{NA} or positive 
#' 			numbers; if NA, the corresponding lengthscale will be set via optimization, 
#' 			otherwise it will be fixed to the value given.
#' @param sparse Should we use the stochastic variational GP.
#' @param force_increasing Boolean; if true, beta diversity will be constrained to 
#' 			increase with environmental distance
#' @param sparse_inducing Number of inducing inputs to use for the svgp
#' @param sparse_batch Batchsize to use for the svgp
#' @param sparse_iter Maximum number of optimizer iterations for svgp
#' @return An mbm object
#' @keywords internal
make_mbm <- function(y, x, y_name, link, likelihood, lengthscale, sparse, force_increasing, 
	sparse_inducing = 10, sparse_batch = 10, sparse_iter = 10000) {
	if(any(rownames(x) != rownames(y)))
		stop("rownames(x) must equal rownames(y)")

	GPy <- reticulate::import("GPy")

	model <- list()
	class(model) <- c('mbm', class(model))
	attr(model, 'y_name') <- y_name
	
	## process covariates
	model$x <- x
	x <- scale(x)
	model$x_scaling = function(xx) scale(xx, center = attr(x, "scaled:center"), 
		scale = attr(x, "scaled:scale"))
	model$x_unscaling = function(xx) xx*attr(x, "scaled:scale") + 
 		attr(x, "scaled:center")
	xDF <- env_dissim(x)
	
	## process response transformation & mean function
	## if a link function is desired with anything but a vanilla GP, we have to cheat
	## by pre-transforming the y-variable
	model$y <- y
	yDF <- reshape2::melt(y,varnames=c('site1', 'site2'), value.name = y_name)
	dat <- merge(xDF, yDF, all.x = TRUE, by=c('site1', 'site2'))
	if((sparse | force_increasing) & link != 'identity') {
		model <- setup_y_transform(model, link)
		link <- "identity"
	} else {
		model$y_transform <- model$y_rev_transform <- function(y) y
	}

	## add data to obj
	x_cols <- which(!colnames(dat) %in% c(y_name))
	model$response <- model$y_transform(as.matrix(dat[,y_name]))
	covars <- dat[,x_cols]
	names <- grep('site', colnames(covars))
	model$covariates <- as.matrix(covars[,-names])
	model$covar_sites <- covars[,names]
	
	##
	## SET UP SVGP
	##
	if(sparse) {
		attr(model, "batchsize") <- sparse_batch
		attr(model, "inducing_inputs") <- sparse_inducing
		attr(model, "svgp_maxiter") <- sparse_iter
		likelihood <- "gaussian"
		model$inducing_inputs <- apply(model$covariates, 2, 
			function(xx) runif(attr(model, "inducing_inputs"), min(xx), max(xx)))
	}

	##
	## SET UP PYTHON OBJECTS
	##
	model$pyobj <- list()
	model <- set_mbm_link(model, link)
	model$likelihood <- likelihood
	model$lengthscale <- lengthscale


	if(model$likelihood == 'gaussian') {
		model$pyobj$likelihood <- GPy$likelihoods$Gaussian(gp_link = model$pyobj$linkFun)
	} else {
		stop("Non-gaussian likelihoods are not supported")
	}

	if(model$likelihood == 'gaussian' & model$link == 'identity') {
		model$pyobj$inference <- GPy$inference$latent_function_inference$ExactGaussianInference()
		if(sparse) {
			attr(model, 'inference') <- 'svgp'
		} else 
			attr(model, 'inference') <- 'exact'
	} else {
		model$pyobj$inference <- GPy$inference$latent_function_inference$Laplace()
		attr(model, 'inference') <- 'laplace'
	}

	model$pyobj$kernel <- setup_mbm_kernel(dim = ncol(model$covariates), 
		lengthscale = model$lengthscale, sparse = sparse)

 	model$pyobj$mf <- set_mean_function(ncol(model$covariates), force_increasing)
 	attr(model, "mean_function") <- if(force_increasing) "increasing" else "0"

	return(model)	
}


#' Set up link functions for mbm objects
#' 
#' This is the only supported way for changing link functions in mbm objects; do not try
#' to do it by hand
#' @param x mbm model object
#' @param link Link function to use
#' @return A copy of the mbm object with the link function set
#' @keywords internal
set_mbm_link <- function(x, link) {
	GPy <- reticulate::import("GPy")
	x$link <- link
	if(link == 'identity') {
		x$pyobj$linkFun <- GPy$likelihoods$link_functions$Identity()
		x$inv_link <- function(y) y
	} else if(link == 'probit') {
		x$pyobj$linkFun <- GPy$likelihoods$link_functions$Probit()
		x$inv_link <- function(y) pnorm(y)
	} else {
		stop(link, 'is an unsupported link function')
	}
	return(x)
}

#' Set up MBM kernel
#' @param dim Kernel dimension (number of variables)
#' @param lengthscale Lengthscale parameter as from [mbm()]. 
#' 		Either NULL (in which case all lengthscales will be optimized) or 
#'		a numeric vector of length \code{ncol(x)+1}. If a vector, the first entry 
#' 		corresponds to environmental distance, and entries \code{i = 1 + (1:n)} 
#' 		correspond to the variable in x[,i]. Values must be \code{NULL} or positive 
#' 		numbers; if NULL, the corresponding lengthscale will be set via optimization, 
#' 		otherwise it will be fixed to the value given.
#' @param prior Prior distribution to use; currently ignored
#' @param sparse Logical, should a sparse GP be used?
#' @param which Which parameters to set up, either all, variance params, or lengthscale
#' @keywords internal
setup_mbm_kernel <- function(dim, lengthscale = NULL, prior, sparse,
		which = c('all', 'lengthscale', 'variance')) {

	reticulate::source_python(system.file("python/kernel.py", package="mbm"))
	k <- make_kernel(dim, sparse)
	k <- set_kernel_constraints_py(k, lengthscale, which)
	return(k)
}

#' Set up MBM mean function
#' @param dim Kernel dimension (number of variables)
#' @param useMeanFunction Logical, should the mean function be used?
#' @keywords internal
set_mean_function <- function(dim, useMeanFunction) {
	reticulate::source_python(system.file("python/mf.py", package="mbm"))
	mf <- set_mean_function_py(dim, useMeanFunction)
	return(mf)
}



#' Set y transformations for an MBM object
#' @details Sets up y transformations to use instead of a link function for probit links. When
#' 		the y-data include zeros and ones, an additional step is necessary to squeeze the data
#' 		from [0,1], (0,1], or [0,1) to the open interval (0,1). In the [0,1] case, the smithson
#' 		lemon-squeezer is used: y􏰂' = [y􏰀(N – 1) + 1/2]/N. Otherwise, eps is first added or
#' 		subtracted from all values.
#' @param x an MBM object
#' @param link character, link function to use
#' @param eps Constant to add to avoid infinite values for probit
#' @references Smithson, M. and Verkuilen, J. 2006. A Better Lemon Squeezer? Maximum-Likelihood
#'		Regression With Beta-Distributed Dependent Variables. Psychological Methods 11(1): 54-71.
#' @keywords internal
setup_y_transform <- function(x, link, eps = 0.001) {
	if(link != 'probit')
		stop("Currently only identity or probit links are supported with these options")

	if(min(x$y) == 0 & max(x$y) == 1) {
		# use the smithson transform when we have both 0s and 1s
		n <- length(model$y)
		x$y_transform <- function(p) qnorm((p * (n - 1) + 0.5) / n)
		x$y_rev_transform <- function(q) (pnorm(q) * n - 0.5) / (n-1)
	} else {
		eps <- if(min(x$y) == 0) eps else if(max(x$y == 1)) -eps else 0
		x$y_transform <- function(p) qnorm(p + eps)
		x$y_rev_transform <- function(q) pnorm(q) - eps
	}
	return(x)
}


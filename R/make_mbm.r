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
#' @param force_increasing Boolean; if true, beta diversity will be constrained to 
#' 			increase with environmental distance
#' 
#' 
#' 
#' 
#' 
# @param predictX List of prediction datasets; each list element is a matrix with same 
#' 			number of columns as \code{x}
# @param scale Boolean, if true, x values will be centered and scaled before fitting the 
#' 			model.
# @param response_curve The type of response curve to generate. The default 
#' 			(\code{distance}) will predict over a range of distances assuming pairs of 
#' 			sites equally spaced across the midpoint of environmental space. \code{none} 
#' 			Produces no response curve, while \code{all} creates response curves for all 
#' 			variables.
# @param svgp Should we use the stochastic variational GP.
# @param svgp_inducing Number of inducing inputs to use for the svgp
# @param svgp_batch Batchsize to use for the svgp
# @param svgp_iter Maximum number of optimizer iterations for svgp
#' @return An mbm object
#' @keywords internal
make_mbm <- function(y, x, y_name, link, likelihood, lengthscale, sparse, force_increasing)
# make_mbm <- function(x, y, y_name, predictX, link, scale, lengthscale, force_increasing, 
				# response_curve, svgp, svgp_inducing=10, svgp_batch=10, svgp_iter=10000)
{
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
	yDF <- reshape2::melt(y,varnames=c('site1', 'site2'), value.name = y_name)
	dat <- merge(xDF, yDF, all.x = TRUE, by=c('site1', 'site2'))

	# model$y_transform <- model$y_rev_transform <- function(y) y


	# if(svgp)
	# {
	# 	if(link != 'identity')
	# 	{
	# 		model <- set_ytrans(model, dat[,y_name], link)
	# 		link <- 'identity'
	# 		dat[,y_name] <- model$y_transform(dat[,y_name])
	# 	}
	# 	attr(model, "inference") <- "svgp"
	# 	attr(model, "batchsize") <- svgp_batch
	# 	attr(model, "inducing_inputs") <- svgp_inducing
	# 	attr(model, "svgp_maxiter") <- svgp_iter
	# } else {
	# 	attr(model, "inference") <- "exact"
	# }

	## add data to obj
	x_cols <- which(!colnames(dat) %in% c(y_name))
	model$response <- as.matrix(dat[,y_name])
	covars <- dat[,x_cols]
	names <- grep('site', colnames(covars))
	model$covariates <- as.matrix(covars[,-names])
	model$covar_sites <- covars[,names]
	
	##
	## SET UP PYTHON OBJECTS
	##
	model$pyobj <- list()

	model$link <- link
	model$likelihood <- likelihood
	model$lengthscale <- lengthscale

	if(model$link == 'identity') {
		linkFun <- GPy$likelihoods$link_functions$Identity()
	} else if(model$link == 'probit') {
		linkFun <- GPy$likelihoods$link_functions$Probit() 
	}

	if(model$likelihood == 'gaussian') {
		model$pyobj$likelihood <- GPy$likelihoods$Gaussian(gp_link = linkFun)
	} else {
		stop("Non-gaussian likelihoods are not supported")
	}

	if(model$likelihood == 'gaussian' & model$link == 'identity') {
		model$pyobj$inference <- GPy$inference$latent_function_inference$ExactGaussianInference()
	} else {
		model$pyobj$inference <- GPy$inference$latent_function_inference$Laplace()
	}

	model$pyobj$kernel <- setup_mbm_kernel(dim = ncol(model$covariates), 
		lengthscale = model$lengthscale, sparse = sparse)

 	model$pyobj$mf <- set_mean_function(ncol(model$covariates), force_increasing)

	return(model)	
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



# #' Set y transformations for an MBM object assuming a linear increasing mean function
# #' @param x an MBM object
# #' @param ydat vector of y data points
# #' @param link character, link function to use
# #' @param eps Constant to add to avoid infinite values for probit
# #' @return An mbm object
# #' @keywords internal
# set_ytrans <- function(x, ydat, link, eps = 0.001)
# {
# 	if(link != 'probit') {
# 		stop('mean function is only supported for probit or identity links')
# 	} else if(min(ydat) < 0 | max(ydat) > 1)
# 		stop('probit data must be on [0,1]')
		
# 	if(min(ydat) == 0 & max(ydat) == 1)
# 	{
# 		# use the smithson transform
# 		forward <- function(p) qnorm( (p * (length(p) - 1) + 0.5)/length(p))
# 		back <- function(q) (pnorm(q) * length(q) - 0.5)/(length(q) - 1)
# 	} else
# 	{
# 		# decide on constant to add depending on whether there are 0s or 1s
# 		eps <- if(min(ydat) == 0) eps else if(max(ydat == 1)) -eps else 0
# 		forward <- function(p) qnorm(p + eps)
# 		back <- function(q) pnorm(q) - eps
# 	}
# 	x$y_transform <- forward
# 	x$y_rev_transform <- back
# 	return(x)
# }

# #' Add \code{link} and \code{rev_link} closures to mbm model
# #' 
# #' If a mean function is used, the link will be set to identity, and y_transform and 
# #' y_untransform methods will also be added.
# #' @param x an mbm object
# #' @param link character; which link function to use
# #' @return \code{x} x with link and rev_link functions added
# #' @keywords internal
# set_link <- function(x, link = c('identity', 'probit', 'log'))
# {
# 	link <- match.arg(link)
# 	if(link == 'identity'){
# 		fun <- unfun <- function(x) x
# 	} else if(link == 'probit') {
# 		fun <- qnorm
# 		unfun <- pnorm
# 	} else if(link == 'log') {
# 		fun <- log
# 		unfun <- exp
# 	}
# 	x$link <- fun
# 	x$rev_link <- unfun
# 	attr(x, "link_name") <- link
# 	return(x)
# }


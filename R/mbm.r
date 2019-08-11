#' Create an mbm model
#' 
#' @param y Square dissimilarity or distance matrix, can be complete or lower triangular 
#' 			only. Row and column names are required and must match the site names in the 
#' 			rows of \code{x}.
#' @param x Matrix giving a series of covariates (in columns) for all sites (in rows). Row 
#' 			names are required. All variables will be included in the model.
#' @param y_name A name to give to the y variable
#' @param link Link function to use
#' @param likelihood Likelihood function to use
# @param n_samples NA or integer; if NA, analytical predictions with standard deviation 
# 			are returned, otherwise posterior samples are returned
# @param response_curve The type of response curve to generate. The default 
# 			(\code{distance}) will predict over a range of distances assuming pairs of 
# 			sites equally spaced across the midpoint of environmental space. \code{none} 
# 			Produces no response curve, while \code{all} creates response curves for all
# 			variables.
# @param pyMsg boolean, should we print messages from python? Useful for debugging
#' @param lengthscale Either NULL (in which case all lengthscales will be optimized) or 
#'		a numeric vector of length \code{ncol(x)+1}. If a vector, the first entry 
#' 		corresponds to environmental distance, and entries \code{i = 1 + (1:n)} 
#' 		correspond to the variable in x[,i]. Values must be \code{NULL} or positive 
#' 		numbers; if NULL, the corresponding lengthscale will be set via optimization, 
#' 		otherwise it will be fixed to the value given.
#' @param sparse Should we use the stochastic variational GP (see 'details').
#' @param force_increasing Boolean; if true, beta diversity will be constrained to 
#' 			increase with environmental distance
#' @param sparse_inducing Number of inducing inputs to use if `sparse = TRUE`
#' @param sparse_batch Batch size to use if `sparse = TRUE`
#' @param sparse_iter Maximum number of optimizer iterations if `sparse = TRUE`
#' @param exact_thresh integer; threshold at which mbm will refuse to run an exact gp.
#' @param verbose Should messages during model fitting be printed?
# @param ... Additional arguments to pass to MBM
#' 
#' @details For larger datasets (more than ~100 sites), it is recommended to use
#' 			`sparse=TRUE`. This will use a sparse approximation to the default
#' 			method, following the stochastical variational GP (Hensman et al 2013).
#' 			Note that if a link function is selected, it will be applied as a 
#' 			transformation of the y data--i.e., for link function L() we fit a SVGP to
#' 			describe the expectation E(L(y))--rather than as a true link function--
#' 			fitting L(E(y))--as is done when `svgp=FALSE`. This is due to a 
#' 			limitation in the underlying GP library.
#' 
#' @return An S3 object of class mbm, containing the following components:
#' 	* `x`: the original (untransformed) site by covariate matrix
#'	* `y`: the original (untransformed) site by site diversity data
#'  * `covariates`: Transformed x-variables supplied to mbm
#'  * `response`: Transformed response variable; this is the data supplied to mbm
#'  * `covar_sites`: Site names to match the covariate matrix
#'  * `y_transform`: transformation applied to y-data before modelling
#'  * `y_rev_transform`: reverse transformation to get y-data back on the original scale
#'  * `link`: a character string identifying the link function
#'  * `inv_link`: inverse of the link function
#'  * `pyobj`: A list of python objects used by the model; this is not meant for user interaction
#' @references Hensman J, Fusi N, and Lawrence ND. 2013. Gaussian Processes for Big Data.
#' 			In: In Proceedings of the 29th Conference on Uncertainty in Artificial 
#' 			Intelligence.
#' @export
# mbm <- function(y, x, predictX, link = c('identity', 'probit', 'log'), scale = TRUE, 
# 				n_samples = NA, response_curve = c('distance', 'none', 'all'),
# 				lengthscale, y_name = 'beta', force_increasing = FALSE, svgp = FALSE,
# 				pyMsg = FALSE, exact_thresh = 100, ...)
mbm <- function(y, x, y_name = 'beta', link = c('identity', 'probit'), likelihood = c('gaussian'),
	lengthscale = NULL, sparse = FALSE, force_increasing = FALSE, sparse_inducing = 10,
	sparse_batch = 10, sparse_iter = 10000, exact_thresh = 100, verbose = FALSE)
{
	link <- match.arg(link)
	likelihood <- match.arg(likelihood)
	# sanity check on sample sizes
	if(!sparse & nrow(y) > exact_thresh) {
		msg <- paste0("A model with n = ", nrow(y), " sites is too large and may run out", 
			" of memory. Try sparse = TRUE. If you really want an exact GP, you can ",
			"increase the exact_thresh parameter.")
		stop(msg) 
	}

	# load python MBM class and dependencies
	GPy <- reticulate::import("GPy")
	pyWarnings <- reticulate::import("warnings")
	# check for warnings
	if(verbose) {
		pyWarnings$simplefilter('default')		
	} else {
		pyWarnings$filterwarnings('ignore')
	}

	model <- make_mbm(y, x, y_name, link, likelihood, lengthscale, sparse, force_increasing, 
		sparse_inducing, sparse_batch, sparse_iter)

	# set up mbm object
	initialize <- TRUE ## no support for non-initialized models yet
	if(attr(model, 'inference') == 'svgp') {
		# climin <- reticulate::import("climin")
		# stop(climin[["__version__"]])
		reticulate::source_python(system.file("python/svgp_optim.py", package="mbm"))
		mod <- run_svgp(X=model$covariates, Y=model$response, 
				Z = model$inducing_inputs, kernel = model$pyobj$kernel, 
				likelihood = model$pyobj$likelihood,  mf = model$pyobj$mf,
				bs = as.integer(attr(model, "batchsize")), init = initialize, 
				maxiter = as.integer(attr(model, "svgp_maxiter")), verbose = verbose)
		model$n_iter <- mod[[2]]
		model$pyobj$gp <- mod[[1]]

		if(model$n_iter >= attr(model, "svgp_maxiter"))
			warning("Max nubmer of iterations reached indicating the model parameters may not", 
				"have converged; try increasing sparse_batch or sparse_iter")
	} else {
		model$pyobj$gp <- GPy$core$GP(X=model$covariates, Y=model$response, 
			kernel = model$pyobj$kernel, likelihood = model$pyobj$likelihood, 
			inference_method = model$pyobj$inference, initialize = initialize, 
			mean_function = model$pyobj$mf)
		model$pyobj$gp$optimize()
	}
	
	# set up parameters and parameter names
	model$params <- model$pyobj$gp$param_array
	names(model$params) <- parse_param_names(model)
	model
}



#' Determine parameter names to a fit mbm model
#' @param x A fit mbm model
#' @return A vector of parameter names that can be assigned to names(x$params)
#' @keywords internal
parse_param_names <- function(x) {
	pnames <- NULL
	## order matters
	## first inducing inputs if svgp
	if(attr(x, 'inference') == 'svgp') {
		pnames <- c(pnames, apply(expand.grid(1:ncol(x$inducing_inputs), 
			1:nrow(x$inducing_inputs)), 1, 
			function(xx) paste0('inducing_input.[', xx[2], ',', xx[1], ']')))
	}
	## 2 is mean function parameters if present
	if(attr(x, "mean_function") == "increasing") {
		pnames <- c(pnames, "prior_intercept", paste0("prior_slope_", colnames(x$covariates)))
	}
	## next is variance and lengthscales
	pnames <- c(pnames, 'rbf_variance', paste0('ls_', colnames(x$covariates)))
	## after lengthscale is white noise variance for svgp
	if(attr(x, 'inference') == 'svgp')
		pnames <- c(pnames, 'White_noise.variance')
	## then gaussiam noise variance
	pnames <- c(pnames, 'gaussian_noise_variance')
	## then cholesky decomposition parameters
	if(attr(x, 'inference') == 'svgp') {
		nz <- nrow(x$inducing_inputs)
		nch <- (nz * (nz+1))/2
		pnames <- c(pnames, paste0("u_cholesky.", 1:nch), paste0("u_mean.", 1:nz))
	}
	pnames
}






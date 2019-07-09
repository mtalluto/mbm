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
#' 			are returned, otherwise posterior samples are returned
# @param response_curve The type of response curve to generate. The default 
#' 			(\code{distance}) will predict over a range of distances assuming pairs of 
#' 			sites equally spaced across the midpoint of environmental space. \code{none} 
#' 			Produces no response curve, while \code{all} creates response curves for all
#' 			variables.
#' @param lengthscale Either NULL (in which case all lengthscales will be optimized) or 
#'		a numeric vector of length \code{ncol(x)+1}. If a vector, the first entry 
#' 		corresponds to environmental distance, and entries \code{i = 1 + (1:n)} 
#' 		correspond to the variable in x[,i]. Values must be \code{NULL} or positive 
#' 		numbers; if NULL, the corresponding lengthscale will be set via optimization, 
#' 		otherwise it will be fixed to the value given.
#' @param pyMsg boolean, should we print messages from python? Useful for debugging
#' @param sparse Should we use the stochastic variational GP (see 'details').
#' @param force_increasing Boolean; if true, beta diversity will be constrained to 
#' 			increase with environmental distance

# @param exact_thresh integer; threshold at which mbm will refuse to run an exact gp.
# @param ... Additional arguments to pass to MBM
#' 
#' @details For larger datasets (more than ~100 sites), it is recommended to use
#' 			\code{sparse=TRUE}. This will use a sparse approximation to the default
#' 			method, following the stochastical variational GP (Hensman et al 2013).
#' 			Note that if a link function is selected, it will be applied as a 
#' 			transformation of the y data--i.e., for link function L() we fit a SVGP to
#' 			describe the expectation E(L(y))--rather than as a true link function--
#' 			fitting L(E(y))--as is done when \code{svgp=FALSE}. This is due to a 
#' 			limitation in the underlying GP library.
#' @return An S3 object of class mbm. 
#' @references Hensman J, Fusi N, and Lawrence ND. 2013. Gaussian Processes for Big Data.
#' 			In: In Proceedings of the 29th Conference on Uncertainty in Artificial 
#' 			Intelligence.
#' @export
# mbm <- function(y, x, predictX, link = c('identity', 'probit', 'log'), scale = TRUE, 
# 				n_samples = NA, response_curve = c('distance', 'none', 'all'),
# 				lengthscale, y_name = 'beta', force_increasing = FALSE, svgp = FALSE,
# 				pyMsg = FALSE, exact_thresh = 100, ...)
mbm <- function(y, x, y_name = 'beta', link = c('identity', 'probit'), likelihood = c('gaussian'),
	lengthscale = NULL, pyMsg = FALSE, sparse = FALSE, force_increasing = FALSE)
{
	link <- match.arg(link)
	likelihood <- match.arg(likelihood)
	# load python MBM class and dependencies
	# reticulate::source_python(system.file("mbm.py", package="mbm"))
	GPy <- reticulate::import("GPy")

	# pyWarnings <- reticulate::import("warnings")
	# # check for warnings
	# if(pyMsg) {
	# 	pyWarnings$simplefilter('default')
		
	# } else {
	# 	pyWarnings$filterwarnings('ignore')
	# }



	# set up mbm object
	model <- make_mbm(y, x, y_name, link, likelihood, lengthscale, sparse, force_increasing)


	model$pyobj$gp <- GPy$core$GP(X=model$covariates, Y=model$response, 
		kernel = model$pyobj$kernel, likelihood = model$pyobj$likelihood, 
		inference_method = model$pyobj$inference, initialize = TRUE, 
		mean_function = model$pyobj$mf)
	model$pyobj$gp$optimize()
	

	# get results into a more r-like format
	model$params <- model$pyobj$gp$param_array
	pnames <- NULL
	if(force_increasing) {
		pnames <- c(pnames, "prior_intercept", paste0("prior_slope_", colnames(model$covariates)))
	}
	pnames <- c(pnames, 'rbf_variance', paste0('ls_', colnames(model$covariates)), 
		'gaussian_noise_variance')
	names(model$params) <- pnames




	# sanity check on sample sizes
	# if(!svgp & nrow(y) > exact_thresh) {
	# 	msg <- paste0("A model with n = ", nrow(y), " sites is too large and may run out", 
	# 		" of memory. Try svgp = TRUE. If you really want an exact GP, you can ",
	# 		"increase the exact_thresh parameter.")
	# 	stop(msg) 
	# }


	model
}







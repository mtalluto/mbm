#' Create an mbm model
#' 
#' @param y Square dissimilarity or distance matrix, can be complete or lower triangular 
#' 			only. Row and column names are required and must match the site names in the 
#' 			rows of \code{x}.
#' @param x Matrix giving a series of covariates (in columns) for all sites (in rows). Row 
#' 			names are required. All variables will be included in the model.
#' @param link Link function to use
#' @param predictX List of prediction datasets; each list element is a matrix with same 
#' 			number of columns as \code{x}. See details.
#' @param scale Boolean, if true, x values will be centered and scaled before fitting the
#' 			model.
#' @param n_samples NA or integer; if NA, analytical predictions with standard deviation 
#' 			are returned, otherwise posterior samples are returned
#' @param response_curve The type of response curve to generate. The default 
#' 			(\code{distance}) will predict over a range of distances assuming pairs of 
#' 			sites equally spaced across the midpoint of environmental space. \code{none} 
#' 			Produces no response curve, while \code{all} creates response curves for all
#' 			variables.
#' @param lengthscale Either missing (in which case all lengthscales will be optimized) or 
#'			a numeric vector of length \code{ncol(x)+1}. If a vector, the first entry 
#' 			corresponds to environmental distance, and entries \code{i = 1 + (1:n)} 
#' 			correspond to the variable in x[,i]. Values must be \code{NA} or positive 
#' 			numbers; if NA, the corresponding lengthscale will be set via optimization, 
#' 			otherwise it will be fixed to the value given.
#' @param y_name A name to give to the y variable
#' @param force_increasing Boolean; if true, beta diversity will be constrained to 
#' 			increase with environmental distance
#' @param svgp Should we use the stochastic variational GP (see 'details').
#' @param GPy_location Optional string giving the location of the user's GPy installaion
#' @param pyCmd Where to look for python; the version in use must have GPy installed
#' @param pyMsg boolean, should we print messages from python? Useful for debugging
#' @param exact_thresh integer; threshold at which mbm will refuse to run an exact gp.
#' @param ... Additional arguments to pass to MBM
#' 
#' @details Prediction datasets can either be supplied when the model is called, or by 
#' 			using the \code{predict} method on the \code{mbm} object. The former will 
#' 			generally be much faster to run; see \code{\link{predict.mbm}}. Note that 
#' 			predictions are always generated for the input dataset.
#' 
#' 			For larger datasets (more than ~100 sites), it is recommended to use
#' 			\code{svgp=TRUE}. This will use a sparse approximation to the default
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
mbm <- function(y, x, predictX, link = c('identity', 'probit', 'log'), scale = TRUE, 
				n_samples = NA, response_curve = c('distance', 'none', 'all'),
				lengthscale, y_name = 'beta', force_increasing = FALSE, svgp = FALSE,
				GPy_location = NA, pyCmd = 'python', pyMsg = FALSE, exact_thresh = 100, ...)
{

	# sanity check on sample sizes
	if(!svgp & nrow(y) > exact_thresh) {
		msg <- paste0("A model with n = ", nrow(y), " sites is too large and may run out", 
			" of memory. Try svgp = TRUE. If you really want an exact GP, you can ",
			"increase the exact_thresh parameter.")
		stop(msg) 
	}
	link <- match.arg(link)
	response_curve <- match.arg(response_curve)
	if(response_curve == 'all')
	{
		warning("response_curve = 'all' is not implemented; switching to 'distance'")
		response_curve <- 'distance'
	}

	model <- make_mbm(x, y, y_name, predictX, link, scale, lengthscale, force_increasing, 
		response_curve, svgp, ...)
	
	# generate temporary files
	files <- write_mbm_dat(model)

	# set up arguments to the python call
	mbmArgs <- make_args(model, files, GPy_location=GPy_location, n_samples = n_samples)

	# write out prediction datasets
	if('predictX' %in% names(model))
	{
		predictResults <- mapply(write_mbm_predict, model$predictX, names(model$predictX))
		predictFiles <- predictResults[1,]

		# concatinate the --predict args
		mbmArgs <- c(mbmArgs, predictResults[2,])
	}

	# run the model
	result <- system2('python', args=mbmArgs, stdout = TRUE)
	if(pyMsg) print(result)
	if("status" %in% names(attributes(result))) 
		stop("MBM returned an error: ", attr(result, "status"))
	# cat(result, '\n')

	# collect results
	params <- data.table::fread(files['params'], sep=',', data.table=FALSE)
	mfslopes <- grepl("mf.slope", params[,1])
	lscales <- grepl("rbf.lengthscale", params[,1])
	repl <- paste0("\\1", colnames(model$covariates))
	if(any(mfslopes)) {
		params[mfslopes,1] <- mapply(function(x,y) sub("(mf\\.slope\\.)([0-9]+)", x, y), 
			repl, params[mfslopes,1])
	}
	if(any(lscales)) {
		params[lscales,1] <- mapply(function(x,y) sub("(rbf\\.lengthscale\\.)([0-9]+)", 
			x, y), repl, params[lscales,1])
	}
	model$params <- params[,2]
	names(model$params) <- params[,1]
	# model$params <- unlist(data.table::fread(files['params'], sep=',', data.table=FALSE))
	# parnames <- c('rbf.variance', paste('lengthscale', colnames(model$covariates), 
	# 	sep='.'), 'noise.variance')
	# if(attr(model, "mean_function") == "linear_increasing") {
	# 	parnames <- c("mf.intercept", paste("mf.slope", colnames(model$covariates), 
	# 		sep='.'), parnames)
	# }
	# names(model$params) <- parnames
	model$linear.predictors <- read_mbm_predict(files['covariates'], nsamp=n_samples)
	model$fitted.values <- if('fit' %in% colnames(model$linear.predictors)) {
			model$y_rev_transform(model$rev_link(model$linear.predictors[,'fit']))
		} else model$y_rev_transform(model$rev_link(model$linear.predictors))

	if("predictX" %in% names(model))
		model$predictions <- lapply(predictFiles, read_mbm_predict, nsamp = n_samples)
	model
}


#' Produce arguments for the python call
#' 
#' @param x An MBM object
#' @param files Character vector of filenames generated from \code{\link{write_mbm_dat}}
#' @keywords internal
#' @return A character vector of arguments to a python call
make_args <- function(x, files, GPy_location = NA, n_samples = NA, tfOutput = '.out')
{
	args <- c(system.file('mbm.py', package='mbm', mustWork = TRUE),
		paste0('--y=', files['response']), paste0('--x=', files['covariates']), 
		paste0('--link=', attr(x, "link_name")),
		paste0('--par=', files['params']), paste0('--out=', tfOutput))
	if(!is.na(GPy_location))
		args <- c(args, paste0('--gpy=', GPy_location))
	if(!is.na(n_samples))
		args <- c(args, paste0('--sample=', n_samples))
	if('fixed_lengthscales' %in% names(x)) 
		args <- c(args, paste0('--ls=', prep_ls(x$fixed_lengthscale)))
	if(attr(x, "mean_function") == "linear_increasing")
		args <- c(args, '--mf')	
	if(attr(x, "inference") == "svgp")
	{
		args <- c(args, '--svgp')
		args <- c(args, paste0('--bs=', attr(x, "batchsize")))
		args <- c(args, paste0('--inducing=', attr(x, "inducing_inputs")))
		args <- c(args, paste0('--svgp_iter=', attr(x, "svgp_maxiter")))
	}
	return(args)
}


#' Convenience function to set up the lengthscale for passing to python
#' @param lengthscale A \code{lengthscale} from an \code{\link{mbm}} object
#' @keywords internal
#' @return A lengthscale string suitable for passing to python
prep_ls <- function(lengthscale) {
	lengthscale[is.na(lengthscale)] <- "nan"
	lengthscale <- paste(lengthscale, collapse=',')
	lengthscale
}

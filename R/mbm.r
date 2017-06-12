#' Create an mbm model
#' 
#' @param y Square dissimilarity or distance matrix, can be complete or lower triangular only. Row and column names are
#'          required and must match the site names in the rows of \code{x}.
#' @param x Matrix giving a series of covariates (in columns) for all sites (in rows). Row names are required. All 
#'          variables will be included in the model.
#' @param link Link function to use
#' @param predictX List of prediction datasets; each list element is a matrix with same number of columns as \code{x}. See details.
#' @param scale Boolean, if true, x values will be centered and scaled before fitting the model.
#' @param n_samples NA or integer; if NA, analytical predictions with standard deviation are returned, otherwise posterior samples are returned
#' @param response_curve The type of response curve to generate. The default (\code{distance}) will predict over a range of distances
#'          assuming pairs of sites equally spaced across the midpoint of environmental space. \code{none} Produces no response curve,
#'          while \code{all} creates response curves for all variables.
#' @param lengthscale Either missing (in which case all lengthscales will be optimized) or a numeric vector of length \code{ncol(x)+1}.
#'         If a vector, the first entry corresponds to environmental distance, and entries \code{i = 1 + (1:n)} correspond to the variable in 
#'         x[,i]. Values must be \code{NA} or positive numbers; if NA, the corresponding lengthscale will be set via optimization, otherwise
#'         it will be fixed to the value given.
#' @param y_name A name to give to the y variable
#' @param force_increasing Boolean; if true, beta diversity will be constrained to increase with environmental distance
#' @param GPy_location Optional string giving the location of the user's GPy installaion
#' @param pyCmd Where to look for python; the version in use must have GPy installed
#' @param pyMsg boolean, should we print messages from python? Useful for debugging
#' 
#' @details Prediction datasets can either be supplied when the model is called, or by using the \code{predict} method on the \code{mbm}
#'          object. The former will generally be much faster to run; see \code{\link{predict.mbm}}. Note that predictions are always
#'          generated for the input dataset.
#' @return An S3 object of class mbm. 
#' @export
mbm <- function(y, x, predictX, link = c('identity', 'probit', 'log'), scale = TRUE, n_samples = NA, response_curve = c('distance', 'none', 'all'),
				lengthscale, y_name = 'beta', force_increasing = FALSE, GPy_location, pyCmd = 'python', pyMsg = FALSE)
{
	link <- match.arg(link)
	response_curve <- match.arg(response_curve)
	if(response_curve == 'all')
	{
		warning("response_curve = 'all' is not implemented; switching to 'distance'")
		response_curve <- 'distance'
	}

	model <- make_mbm(x, y, y_name, predictX, link, scale, lengthscale, force_increasing, response_curve)

	tfBase <- "mbm_"
	tfExt <- ".csv"
	tfOutput <- '_out.csv'
	
	# generate temporary files
	yFile <- tempfile(paste0(tfBase, 'y_'), fileext=tfExt)
	xFile <- tempfile(paste0(tfBase, 'x'), fileext=tfExt)
	data.table::fwrite(as.data.frame(model$response), yFile)
	data.table::fwrite(model$covariates, xFile)
	parFile <- tempfile(paste0(tfBase, 'par_'), fileext = tfExt)

	# set up arguments to the python call
	mbmArgs <- c(system.file('mbm.py', package='mbm', mustWork = TRUE), # the file name of the python script
				paste0('--y=', yFile), paste0('--x=', xFile), paste0('--link=', attr(model, "link_name")), paste0('--par=', parFile),
				paste0('--out=', tfOutput))  # additional arguments
	if(!missing(GPy_location))
		mbmArgs <- c(mbmArgs, paste0('--gpy=', GPy_location))
	if(!is.na(n_samples))
		mbmArgs <- c(mbmArgs, paste0('--sample=', n_samples))
	if('fixed_lengthscales' %in% names(model)) mbmArgs <- c(mbmArgs, paste0('--ls=', prep_ls(model$fixed_lengthscale)))
	if(attr(model, "mean_function") == "linear_increasing")
		mbmArgs <- c(mbmArgs, paste0('--mf'))
	
	# write out prediction datasets
	if('predictX' %in% names(model))
	{
		predictFiles <- sapply(names(model$predictX), function(nm) tempfile(paste0(tfBase, 'pr_', nm, '_'), fileext=tfExt))
		mapply(function(fname, dat) {
			data.table::fwrite(as.data.frame(dat), fname)
		}, predictFiles, model$predictX)
		# concatinate the --predict args
		mbmArgs <- c(mbmArgs, sapply(predictFiles, function(prf) paste0('--pr=', prf)))
	}

	# run the model
	result <- system2('python', args=mbmArgs, stdout = TRUE)
	if(pyMsg) print(result)
	if("status" %in% names(attributes(result))) 
		stop("MBM returned an error: ", attr(result, "status"))
	
	# collect results
	model$params <- unlist(data.table::fread(parFile, sep=',', data.table=FALSE))
	parnames <- c('rbf.variance', paste('lengthscale', colnames(model$covariates), sep='.'), 'noise.variance')
	if(attr(model, "mean_function") == "linear_increasing") 
		parnames <- c("mf.intercept", paste("mf.slope", colnames(model$covariates), sep='.'), parnames)
	names(model$params) <- parnames
	model$linear.predictors <- get_predicts(paste0(xFile, tfOutput), n_samples)
	model$fitted.values <- if('fit' %in% colnames(model$linear.predictors)) {
			model$rev_link(model$linear.predictors[,'fit']) 
		} else model$rev_link(model$linear.predictors)

	if("predictX" %in% names(model))
		model$predictions <- lapply(predictFiles, function(fname) get_predicts(paste0(fname, tfOutput), n_samples))

	model
}

# just a wrapper for data.table that handles column names
get_predicts <- function(fname, nsamp)
{
	if(is.na(nsamp)) {
		colnames = c('fit', 'stdev')
	} else {
		colnames = paste0('samp', 1:nsamp)
	}
	data.table::fread(fname, sep=',', data.table=FALSE, col.names = colnames)
}



# convenience function to set up the lengthscale for passing to python
prep_ls <- function(lengthscale) {
	lengthscale[is.na(lengthscale)] <- "nan"
	lengthscale <- paste(lengthscale, collapse=',')
	lengthscale
}
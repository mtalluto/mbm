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
#' @param GPy_location Optional character giving the location of the user's GPy installaion
#' @param pyCmd Where to look for python; the version in use must have GPy installed
#' @param pyMsg boolean, should we print messages from python? Useful for debugging
#' 
#' @details Prediction datasets can either be supplied when the model is called, or by using the \code{predict} method on the \code{mbm}
#'          object. The former will generally be much faster to run; see \code{\link{predict.mbm}}. Note that predictions are always
#'          generated for the input dataset.
#' @return An S3 object of class mbm. 
#' @export
mbm <- function(y, x, predictX, link = c('identity', 'probit', 'log'), scale = TRUE, n_samples = NA, response_curve = c('distance', 'none', 'all'),
				lengthscale, y_name = 'beta', GPy_location, pyCmd = 'python', pyMsg = FALSE)
{
	link <- match.arg(link)
	response_curve <- match.arg(response_curve)
	model <- list()
	model$link <- set_link(link)
	model$rev_link <- set_unlink(link)
	
	if(response_curve == 'all')
	{
		warning("response_curve = 'all' is not implemented; swithincg to 'distance'")
		response_curve <- 'distance'
	}
	
	class(model) <- c('mbm', class(model))
	attr(model, 'y_name') <- y_name

	if(scale) {
		x <- scale(x)
		model$x_scaling = function(xx) scale(xx, center = attr(x, "scaled:center"), scale = attr(x, "scaled:scale"))
		model$x_unscaling = function(xx) xx*attr(x, "scaled:scale") + attr(x, "scaled:center")
		if(!missing(predictX))
		{
			if(!is.list(predictX))
				predictX <- list(predictX)
			if(scale) {
				warning("Prediction datasets will be scaled to the same scale as x")
				predictX <- lapply(predictX, model$x_scaling)
			}
		}
	}
	if(!missing(predictX))
	{
		if(is.null(names(predictX)))
			names(predictX) <- 1:length(predictX)
		predictX <- lapply(predictX, env_dissim, sitenames = FALSE)
	}
	
	xDF <- env_dissim(x)

	yDF <- reshape2::melt(y,varnames=c('site1', 'site2'), value.name = y_name)
	dat <- merge(xDF, yDF, all.x = TRUE)
	x_cols <- which(!colnames(dat) %in% c('site1', 'site2', y_name))
	model$response <- dat[,y_name]
	model$covariates <- dat[,x_cols]
	
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
	mbmArgs <- c(system.file('mbm.py', package='mbmtools', mustWork = TRUE), # the file name of the python script
				paste0('--y=', yFile), paste0('--x=', xFile), paste0('--link=', link), paste0('--par=', parFile),
				paste0('--out=', tfOutput))  # additional arguments
	if(!missing(GPy_location))
		mbmArgs <- c(mbmArgs, paste0('--gpy=', GPy_location))
	if(!is.na(n_samples))
		mbmArgs <- c(mbmArgs, paste0('--sample=', n_samples))
	
	if(!missing(lengthscale))
	{
		if(length(lengthscale) != ncol(model$covariates) | !(all(lengthscale > 0 | is.na(lengthscale))))
			stop("Invalid lengthscale specified; see help file for details")
		lengthscale[is.na(lengthscale)] <- "nan"
		lengthscale <- paste(lengthscale, collapse=',')
		mbmArgs <- c(mbmArgs, paste0('--ls=', lengthscale))
	}
	
	# set up response curve
	if(response_curve == 'distance')
	{
		rcX <- rc_data(model, 'distance')
		if(missing(predictX)) {
			predictX <- rcX
		} else
			predictX <- c(rcX, predictX)
	}
	# write out prediction datasets
	if(!missing(predictX))
	{
		model$predictX <- predictX
		predictFiles <- sapply(names(predictX), function(nm) tempfile(paste0(tfBase, 'pr_', nm, '_'), fileext=tfExt))
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

# convenience functions for getting the link transformations
set_link <- function(link)
{
	if(link == 'identity'){
		fun <- function(x) x
	} else if(link == 'probit') {
		fun <- qnorm
	} else if(link == 'log') {
		fun <- log
	} else 
		stop("unknown link ", link)
	return(fun)
}

set_unlink <- function(link)
{
	if(link == 'identity'){
		fun <- function(x) x
	} else if(link == 'probit') {
		fun <- pnorm
	} else if(link == 'log') {
		fun <- exp
	} else 
		stop("unknown link ", link)
	return(fun)
}
#' Predict method for MBM objects
#' 
#' @param x A previously-fit MBM object
#' @param newdata Optional dataset for prediction. If present, it should be either a character vector giving
#'    the name of one of the datasets specified for \code{predictX} when the model was fit, or  a new dataset
#'    in the same format used to fit the model (i.e., a site by covariate matrix). If missing, predictions 
#'    will be for the original data.
#' @param n_samples NA or integer; if NA, analytical predictions with standard deviation are returned, otherwise posterior samples are returned
#' @param GPy_location Optional string giving the location of the user's GPy installaion
#' @param pyMsg boolean, should we print messages from python? Useful for debugging
#' @details Prediction to new data is possible after the fact for mbm models, however there are significant
#'     performance penalties for doing so. Thus, whenever possible, it is preferable to predict during
#'     model fitting via the \code{predictX} argument to the \link{\code{mbm}} function. All prediction
#'     is done on the response scale.
#'     This function caches to disk, thus it is important to ensure that adequate disk space is
#'     available when using large prediction datasets.
#' @return A site by site matrix of predictions 
#' @export
predict.mbm <- function(x, newdata, n_samples = NA, GPy_location = NA, pyMsg = FALSE)
{
	tfOutput <- '_out.csv'

	# steps
	# 0. check new data to see what we have to do
	if(missing(newdata))
	{
		preds <- x$fitted.values
	} else if(is.character(newdata))
	{
		preds <- x$predictions[[newdata]]
	} else {

		# 1. parse x to recreate an mbm call similar to the mbm function
		files <- write_mbm_dat(x)
		# 2. set up mbm call to be a re-launch, not a new model
		mbmArgs <- make_args(x, files, GPy_location=GPy_location, n_samples = n_samples)
		mbmArgs <- c(mbmArgs, '--resume')

		# 3. parse newdata into predictX format
		dat <- prep_predict(newdata, x)
		prFile <- write_mbm_predict(dat)
		mbmArgs <- c(mbmArgs, sapply(prFile, function(prf) paste0('--pr=', prf)))

		# 4. run model
		result <- system2('python', args=mbmArgs, stdout = TRUE)
		if(pyMsg) print(result)
		if("status" %in% names(attributes(result))) 
			stop("MBM returned an error: ", attr(result, "status"))

		# 5. post-process
		preds <- read_mbm_predict(prFile, nsamp = n_samples)
	}
	return(preds)
}

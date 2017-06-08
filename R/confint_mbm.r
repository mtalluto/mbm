#' Confidence intervals for MBM models
#' 
#' Produce a plot of a response curve from an MBM object, with some sensible defaults.
#' If a response curve was not generated when the model was fit, teh default is to try to make an empirical curve
#' from the fits.
#' 
#' @param object An MBM object
#' @param parm Specification of what parameters to consider. Default \code{'fits'} or a value of 0 gives confidence intervals for fitted data points.
#'       Other strings (or numerical index) will be looked up in the model predictions.
#' @param level Confidence level required
#' @param method How to compute confidence interval. \code{'auto'} (the default) uses parametric estimates of the confidence intervals;
#'        this is the only method currently implemented.
#' @param col_line Line color
#' @param col_polygon Polygon color
#' @param cex_pt Size of point plotting symbols
#' @param ... Additional arguments to be passed to base graphics plotting commands
#' @return A matrix with upper and lower confidence limits
#' @export
confint.mbm <- function(object, parm, level = 0.95, method = c('auto', 'parametric', 'sample'))
{
	method <- match.arg(method)
	if(method == 'sample')
	{
		warning("Sample not implemented; changing to method = 'parametric'")
		method <- 'parametric'
	}
	if(method == 'auto') method <- 'parametric'
	if(method == 'parametric')
		ci_method <- ci_parametric
	
	if(parm == 'fits' | parm ==0)
	{
		ci <- ci_method(object$linear.predictors, level)
	} else
	{
		ci <- ci_method(object$predictions[[parm]], level)
	}
	object$y_rev_transform(object$rev_link(ci))
}

#' Parametric confidence intervals for mbm objects
#' 
#' Computed via normal approximation. 
#' @param object A fit from an mbm object; first column should be the prediction and second the standard deviation
#' @param level Confidence level required
#' @return A matrix with upper and lower confidence limits in columns and observations in rows
ci_parametric <- function(x, level)
{
	quants <- qnorm(c((1 - level)/2, 1 - (1 - level)/2))
	cbind(lower = x[,1] + x[,2] * quants[1], upper = x[,1] + x[,2] * quants[2])
}
#' Confidence intervals for MBM models
#' 

#' 
#' @param object An MBM object
#' @param parm Specification of what parameters to consider. Currently, only `fits` is 
#' 		supported; this produces confidence intervals for the predictions of the original
#' 		data.
#' @param level Confidence level required
#' @param method How to compute confidence interval. The default uses a parametric estimate;
#'        this is the only method currently implemented.
#' @return A matrix with upper and lower confidence limits
#' @export
confint.mbm <- function(x, parm='fits', level = 0.95, method = c('parametric', 'sample'))
{
	method <- match.arg(method)
	if(method == 'sample') {
		warning("Sample not implemented; changing to method = 'parametric'")
		method <- 'parametric'
	}
	if(method == 'parametric')
		ci_method <- ci_parametric
	
	if(parm == 'fits') {
		ci <- ci_method(predict(x), level)
	} else {
		stop("Only parm = 'fits' is currently supported")
	}
	return(ci)
}

#' Parametric confidence intervals for mbm objects
#' 
#' Computed via normal approximation. 
#' @param object A fit from an mbm object; a 'fits' column and a 'stdev' column is required
#' @param level Confidence level required
#' @keywords internal
#' @return A matrix with upper and lower confidence limits in columns and observations in rows
ci_parametric <- function(x, level)
{
	quants <- qnorm(c((1 - level)/2, 1 - (1 - level)/2))
	cbind(lower = x[,'fit'] + x[,'stdev'] * quants[1], upper = x[,'fit'] + x[,'stdev'] * quants[2])
}

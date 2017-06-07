#' Plot predictions vs fits for an mbm object
#' 
#' @param x An \code{\link{mbm}} object
#' @param line Boolean; should we show the 1:1 line?
#' @param sterr Boolean; if true, standard errors will be shown with the points
#' @param ... Additional arguments to be passed to base graphics plotting commands
#' @export
plot.mbm <- function(x, line = TRUE, sterr = FALSE, ...)
{
	xx <- x$response
	yy <- x$fitted.values
	# set some defaults if not overridden
	args <- list(...)
	dlims <- range(xx, 
				   if(sterr) c(x$rev_link(x$linear.predictors[,1] + x$linear.predictors[,2]), 
				   			x$rev_link(x$linear.predictors[,1] - x$linear.predictors[,2])) else yy)
	args <- add_default(args, 'ylim', dlims)
	args <- add_default(args, 'xlim', dlims)
	args <- add_default(args, 'xlab', 'Response')
	args <- add_default(args, 'ylab', 'Fitted Values')
	args <- add_default(args, 'lty', 2)
	
	do.call(plot, c(list(x=xx, y=yy), args))
	if(line) do.call(abline, c(list(a=0, b=1), args))
	if(sterr)
	{
		lineArgs = list(x0=x$response, x1=x$response, y0 = x$rev_link(x$linear.predictors[,1] + x$linear.predictors[,2]),
						y1 = x$rev_link(x$linear.predictors[,1] - x$linear.predictors[,2]))
		args$lty <- 1
		do.call(segments, c(lineArgs, args))
	}
}
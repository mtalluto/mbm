#' Standard R methods for mbm objects
#' 
#' @name plot.mbm
#' @aliases print.mbm
#' @aliases summary.mbm
#' @aliases is.mbm
#' @aliases format.mbm
#' @aliases coef.mbm
#' @param x An \code{\link{mbm}} object
#' @param line Boolean; should we show the 1:1 line?
#' @param sterr Boolean; if true, standard errors will be shown with the points
#' @param ... Additional arguments to be passed to base graphics plotting commands
#' @rdname methods
#' @seealso For mbm objects fit with \code{svgp=TRUE}, \code{\link{inducing}} for the 
#' 		matrix of inducing inputs, \code{\link{gp_params}} for the mean vector
#' 		and the variance-covariance matrix of the GP
#' @export
plot.mbm <- function(x, line = TRUE, sterr = FALSE, ...)
{
	ytrans <- function(yy) x$y_rev_transform(x$inv_link(yy))
	xx <- x$y_rev_transform(x$inv_link(x$response))
	yy <- predict(x)
	ymu <- ytrans(yy[,1])

	dlims <- range(ymu)
	if(sterr) {
		lineArgs = list(x0=xx, x1=xx, 
			y0 = ytrans(yy[,1] + yy[,2]),
			y1 = ytrans(yy[,1] - yy[,2]))
		dlims <- range(c(lineArgs$y0, lineArgs$y1))
	}
	
	# set some defaults if not overridden
	args <- list(...)
	dlims <- range(xx, if(sterr) c(ytrans(x$linear.predictors[,1] + 
		x$linear.predictors[,2]), ytrans(x$linear.predictors[,1] - 
		x$linear.predictors[,2])) else ymu)
	args <- add_default(args, 'ylim', dlims)
	args <- add_default(args, 'xlim', dlims)
	args <- add_default(args, 'xlab', 'Response')
	args <- add_default(args, 'ylab', 'Fitted Values')
	args <- add_default(args, 'lty', 2)
	
	do.call(plot, c(list(x=xx, y=ymu), args))
	if(line) do.call(abline, c(list(a=0, b=1), args))
	if(sterr)
	{
		args$lty <- 1
		do.call(segments, c(lineArgs, args))
	}
}

#' @rdname methods
#' @export
coef.mbm <- function(x)
{
	x$params[!grepl("(inducing|u_mean|u_cholesky)", names(x$params))]
}

#' @rdname methods
#' @export
print.mbm <- function(x)
{
	cat(format(x), sep='\n')
}

#' @rdname methods
#' @export
summary.mbm <- function(x) print(x)

#' @rdname methods
#' @export
is.mbm <- function(x) inherits(x, 'mbm')

#' @rdname methods
#' @export
format.mbm <- function(x)
{
	pars <- coef(x)
	c(paste("MBM model on ", ncol(x$covariates) - 1, "variables"), 
		paste(format(names(pars)), format(pars, digits=2)))
}


#' Standard R methods for mbmSP objects
#' 
#' @name plot.mbmSP
#' @aliases print.mbmSP
#' @aliases summary.mbmSP
#' @aliases is.mbmSP
#' @param x An \code{\link{mbmSP}} object
#' @param sterr Boolean; if true, standard errors will be mapped instead of fits
#' @param ... Additional arguments to be passed to plotting commands
#' @rdname spmethods
#' @export
plot.mbmSP <- function(x, sterr = FALSE, ...)
{
	col_default <- if(requireNamespace('viridis', quietly = TRUE)) {
		viridis::magma(100)
	} else heat.colors(100)
	if(sterr)
	{
		args <- list(...)
		args <- add_default(args, 'col', col_default)
		ras <- raster::raster(x$stdev)
		do.call(raster::plot, c(x=ras, args))
	} else {
		ras <- raster::stack(x$fits)
		raster::plotRGB(ras, scale = 1, ...)
	}
}


#' @rdname spmethods
#' @export
print.mbmSP <- function(x) print(x$pcoa)

#' @rdname spmethods
#' @export
summary.mbmSP <- function(x) summary(x$pcoa)

#' @rdname spmethods
#' @export
is.mbmSP <- function(x) inherits(x, 'mbmSP')


#' Matrix of inducing inputs for an MBM object
#' @param An \code{\link{mbm}} object
#' @return A matrix of inducing inputs
#' @export
inducing <- function(x)
{
	mat <- matrix(x$params[grepl("inducing", names(x$params))], 
		ncol = ncol(x$covariates), byrow=TRUE)
	colnames(mat) <- colnames(x$covariates)
	mat
}

#' GP parameters for an MBM object
#' @param An \code{\link{mbm}} object
#' @return A list with 2 named items; \code{$mean} is the mean vector and \code{$cholesky}
#' 		is the Cholesky decomposition of the variance-covariance matrix of the SVGP.
#' @export
gp_params <- function(x)
{
	if(attr(x, 'inference') != 'svgp')
		stop("gp_params() is only defined for sparse gps")
	meanvec <- as.vector(x$params[grepl('u_mean', names(x$params))])
	cholesky <- matrix(NA, nrow = length(meanvec), ncol=length(meanvec))
	cholesky[upper.tri(cholesky, diag=TRUE)] <- x$params[grepl('u_cholesky', names(x$params))]
	cholesky <- t(cholesky)
	list(mean=meanvec, cholesky =cholesky)
}

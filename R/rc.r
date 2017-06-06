#' Generate x-data for response curves
#' 
#' Currently this is very simple, does not handle the 'all' case and it assumes covariates are 0-centered
#' 
#' @param x An \code{\link{mbm}} object; it is not necessary for it to be complete (i.e., fitted)
#' @param type The kind of response curve(s) to make
#' @param res How many points per response curve
#' 
#' @return A list of data frames with x-values to be fit to response curves
rc_data <- function(x, type=c('distance', 'all'), res=200)
{
	type <- match.arg(type)
	if(type == 'all')
		stop("type=='all' is not implemented")

	dat <- x$covariates
	rcDat <- matrix(0, ncol=ncol(dat), nrow=res, dimnames = list(1:res, colnames(dat)))
	rcDat[,'distance'] <- seq(0, max(dat[,'distance']), length.out = res)
	
	list(rc_distance = rcDat)
}


#' Plot response curves
#' 
#' Produce a plot of a response curve from an MBM object, with some sensible defaults.
#' If a response curve was not generated when the model was fit, teh default is to try to make an empirical curve
#' from the fits.
#' 
#' @param x An MBM object
#' @param missing_action What to do if no response curve was fit with the model. \code{empirical} will make an empirical
#'       curve from the model fits. \code{fit} will fit a response curve, which can be slow. \code{error} (the only implemented
#'       option) quits with an error message.
#' @param rc_name A character string, integer, or NA. If a string or integer, the name or index of the response curve to fit. If NA,
#'       (not yet implemented), will attempt to find all response curves and plot them all.
#' @param col_pt Point color
#' @param col_line Line color
#' @param col_polygon Polygon color
#' @param cex_pt Size of point plotting symbols
#' @param ... Additional arguments to be passed to base graphics plotting commands
#' @export
rc <- function(x, missing_action = c('error', 'empirical', 'fit'), rc_name = 'rc_distance', col_pt = '#666666', col_line = '#ce5336', 
			   col_polygon = paste0(col_line, '66'), cex_pt = 0.4, ...)
{
	if(is.na(rc_name))
		stop("Plotting all curves is not yet implemented")

	# see if the response curve is present and handle it if not
	if((!('predictions' %in% names(x)))|| (!(rc_name %in% names(x$predictions))))
	{
		# no predictions made
		stop("Response curve must be present to proceed")
	}

	varname <- gsub('rc_(.+)', '\\1', rc_name)
	
	xx <- x$predictX[[rc_name]][,varname]	
	yy <- x$rev_link(x$predictions[[rc_name]][,1])
	datX <- x$covariates[,varname]
	datY <- x$response
	
	interval <- confint(x, rc_name)
	
	# set some defaults if not overridden
	args <- list(...)
	args <- add_default(args, 'ylim', range(c(yy, datY)))
	args <- add_default(args, 'xlim', range(c(xx, datX)))
	args <- add_default(args, 'pch', 20)
	args <- add_default(args, 'xlab', rc_name)
	args <- add_default(args, 'ylab', attr(model, 'y_name'))
	
	do.call(plot, c(list(x=0, y=0, type='n'), args))
	do.call(points, c(list(x=datX, y=datY, col=col_pt, cex=cex_pt), args))
	polyList <- list(x = c(xx, rev(xx)), y = c(interval[,1], rev(interval[,2])), col=col_polygon, border=NA)
	do.call(polygon, c(polyList, args))
	do.call(lines, c(list(x=xx, y=yy, col=col_line),args))
}

add_default <- function(args, x, val)
{
	if(!(x %in% names(args)))
		args[[x]] <- val
	args
}

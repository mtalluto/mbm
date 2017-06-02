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
	rcDat[,'distance'] <- seq(0, max(dat[,'distance']))
	
	list(rc_distance = rcDat)
}
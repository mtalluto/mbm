#' Community weighted mean trait values
#' 
#' @param communities Site by species weights matrix; if presence-absence (1-0), the function returns the community mean
#' @param traits Matrix of traits; the rows MUST be in the same order as the columns of community
#' @param na.rm Logical, whether to remove NAs when computing the mean
#' @return site by trait CWM matrix
#' @export
cwm <- function(communities, traits, na.rm=FALSE)
{
	apply(traits, 2, function(x) colMeans(t(communities) * x, na.rm=na.rm))
}

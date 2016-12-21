#' Compute simpson diversity
#' 
#' @param x A matrix with rows as sites and species as columns. Values are abundances
#' @param proportion Logical. If true, abundances in x are considered proportional abundances. If false, x will be
#'                   transformed to proportional abundance
#' @value A vector of diversity values
#' @export
simpson <- function(x, proportion=TRUE) {
	if((proportion & any(x > 1) | any(x < 0))) stop("Abundances must be >= 0 and must be <= 1 if proportion=TRUE")
	if(! proportion) {
		N <- matrix(rowSums(x), nrow=nrow(x), ncol=ncol(x))
		x <- x/N
	}
	apply(x, 1, function(xx) sum(xx^2))
}


#' Compute species richness
#' 
#' @param x A matrix with rows as sites and species as columns. Values are abundances or presence/absence
#' @value A vector of diversity values
#' @export
richness <- function(x) {
	if(any(x < 0)) stop("Values in x must be positive")
	x[x > 0] <- 1
	x[x != 1] <- 0
	rowSums(x)
}
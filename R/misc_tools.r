#' @name smithson
#' @rdname smithson
#' @title Smithson transformation
#' 
#' @param x A vector
#' @return Transformation of x
#' @references Smithson M and Verkuilen J. 2006. A Better Lemon Squeezer? Maximum-Likelihood
#'               Regression With Beta-Distributed Dependent Variables. Psych Meth 11:54-71.
#' @export
smithson <- function(x) {
	N <- length(x)
	(x * (N-1) + 0.5)/N
}

#' @rdname smithson
#' @export
rev_smithson <- function(x) {
	N <- length(x)
	(N*x - 0.5)/(N-1)
}
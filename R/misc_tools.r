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


# order rows or columns of x based on the row/columns names in y
#' @param x,y matrices to match
#' @param by.x should we match rows, columns, or both?
match_mat_dims <- function(x, y, by.x = c('r', 'c', 'rc'), by.y = c('r', 'c'))
{
	by.x <- match.arg(by.x)
	by.y <- match.arg(by.y)
	
	if(by.y == 'c') y <- t(y)
	if(grepl('r', x)) 
	{
		ind <- match(rownames(y), rownames(x))
		x <- x[ind,]
	}
	if(grepl('c', x)) 
	{
		ind <- match(rownames(y), colnames(x))
		x <- x[,ind]
	}
	x
}

# make sure row and columns of a square matrix match
# by: determines whether row or column order is preserved; if neither, then first we sort by rows
order_symmetric_mat <- function(x, by=c('rows', 'columns', 'neither'))
{
	if(nrow(x) != ncol(x)) stop("matrix must be square")
	if((!all(rownames(x) %in% colnames(x))) | !(all(colnames(x) %in% rownames(x))))
		stop("all row and column names must have a match")
	
	by <- match.arg(by)
	
	if(by == 'cols')
		x <- t(x)
	if(by == 'neither')
	{
		ind <- order(rownames(x))
		x <- x[ind,]
	}
	
	ind <- match(rownames(x), colnames(x))
	x <- x[,ind]
	
	if(by == 'cols') x <- t(x)
	x
}


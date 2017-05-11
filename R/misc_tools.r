
# order rows or columns of x based on the row/columns names in y
#' @title Match matrix dimensions
#' @param x,y matrices to match
#' @param by.x should we match rows, columns, or both?
match_mat_dims <- function(x, y, by.x = c('r', 'c', 'rc'), by.y = c('r', 'c'))
{
	by.x <- match.arg(by.x)
	by.y <- match.arg(by.y)
	
	if(by.y == 'c') y <- t(y)
	if(grepl('r', by.x)) 
	{
		ind <- match(rownames(y), rownames(x))
		x <- x[ind,]
	}
	if(grepl('c', by.x)) 
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


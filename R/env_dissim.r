#' Compute a dataframe of environmental dissimilarity values
#' 
#' @param x A data frame of environmental covariates.
#' @param sites Either a vector of site names or a column index (or column name) giving the location of the site names
#'              in x. If 0, row names will be used
#' @param sitenames Boolean; should site names be returned with the data frame?
#' @return A data frame including both sites' names, environmental dissimilarity between the two sites, and midpoints
#'         between the two sites for each variable
#' @export
env_dissim <- function(x, sites = 0, sitenames = TRUE) {
	if(length(sites) == nrow(x)) {
		rownames(x) <- sites
	} else if(sites != 0) {
		rownames(x) <- x[,sites]
		x <- x[,-sites]
	}

	midpoints <- apply(x, 2, function(v) {
		mat <- sapply(v, function(x) (x+v)/2)
 		mat[lower.tri(mat)]
	})
	distance <- as.vector(dist(x))
	covars <- as.data.frame(cbind(distance, midpoints))
	colnames(covars) <- c('distance', colnames(x))
	
	if(sitenames)
	{
		nms <- t(do.call(cbind, sapply(1:(nrow(x) - 1), function(i) {
			sapply((i+1):nrow(x), function(j) {
				c(rownames(x)[i], rownames(x)[j])
			}) 
		})))
		covars$site1 <- nms[,1]
		covars$site2 <- nms[,2]
	}
	covars
}
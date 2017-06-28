#' @name bc
#' @aliases jaccard
#' @title Compute dissimilarity indices
#' @param comm Community matrix giving sites (rows) by species (columns); values are abundances. Column and row labels are highly recommended
#' @return Site by site dissimilarity matrix
#' @rdname bc
#' @export
bc <- function(comm)
{
	# get min abundance of shared species
	# pretty slow right now because it uses nested applies
	# a fast way to do this would be great
	shared_abundance <- apply(comm, 1, function(site1) {
		apply(comm, 1, function(site2) {
			mask <- site1 <= site2
			sum((site1 * mask) + ((!mask) * site2))
		})
	})
	
	# total abundance at each site for each pairwise combo
	site1 <- matrix(rowSums(comm), nrow=nrow(comm), ncol=nrow(comm))
	site2 <- t(site1)
	
	1.0 - ((2.0 * shared_abundance)/(site1 + site2))
}



#' @rdname bc
#' @export
jaccard <- function(comm)
{
	# get the intersection of all sites - number of sp in common
	site_similarity <- comm %*% t(comm)
	richness <- rowSums(comm)
	# make richness matrices for comparisons
	ra <- matrix(richness, nrow=length(richness), ncol=length(richness))
	1 - (site_similarity / (ra + t(ra) - site_similarity))
}


#' Compute a dataframe of environmental dissimilarity values
#' 
#' @param x A data frame of environmental covariates.
#' @param sites Either a vector of site names or a column index (or column name) giving the location of the site names
#'              in x. If 0, row names will be used; if rownames are missing rows will be numbered
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
	if(all(is.null(rownames(x))))
		rownames(x) <- 1:nrow(x)
	
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
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

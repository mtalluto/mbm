#' Sorensen phylogenetic distance
#' 
#' @param community Community matrix giving sites (rows) by species (columns). Column names should match the
#'                   names in the phylogeny
#' @param binomial Should we return the data in binomial format (matrix of X successes and matrix of N trials)
#' @details The index is defined as 1 - (2 * (A ∩ B) / (A + B)). For binomial data, we return (A + B) - (2 * (A ∩ B)) as the 
#' success matrix, so that a binomial regression will make predictions on the dissimilarity (rather than the similarity)
#' @return If binomial is FALSE, matrix of pairwise dissimilarity values between communities. If binomial is TRUE, 
#'             list of matrices, the first being the numerator of the sorensen index and the second being the denominator
#' @export
sorensen <- function(community, binomial=FALSE)
{
	if(!is.matrix(community)) community <- as.matrix(community)

	# convert to matrix of 0 and 1	
	if(any(! community %in% c(0,1) )) community <- 1 * (community > 0)
	# get the intersection of all sites - number of sp in common
	# this is the 0.5 * numerator of sorensen similarity
	site_similarity <- community %*% t(community)
	
	# compute species richness and repeat into a matrix
	richness <- rowSums(community)
	ra <- matrix(richness, nrow=length(richness), ncol=length(richness))
	
	# denominator - ra is the number of spp at site1, t(ra) is the number at site 2
	denom <- ra + t(ra)
	dimnames(denom) <- dimnames(site_similarity)
	
	if(binomial)
	{
		# return the numerator for *dissimilarity* (which is 1 - similarity, hence denom - sim in the top)
		list(x = denom - 2*site_similarity, n = denom)
	} else {
		1 - ((2 * site_similarity) / denom)
	}
}
#' Sorensen dissimilarity index
#' 
#' @param comm Community matrix giving sites (rows) by species (columns).
#' @param dist Optional species by species distance matrix; margin names must match colnames in comm.
#' @references Pavoine, S. and Ricotta, C. 2014. Functional and phylogenetic similarity among communities.
#' Methods in Ecology and Evolution 5(7): 666–675.
#' @return Matrix of pairwise dissimilarity values between communities.
#' @export
sorensen <- function(comm, dis)
{
	if(missing(dis))
		return(sor(comm))

	if((nrow(dis) != ncol(dis)) | any((colnames(dis) != rownames(dis))))
		stop("dis must be square with identical row and column names")
	
	taxa <- sort(intersect(colnames(comm), rownames(dis)))
	if(length(taxa) < ncol(comm)) 
		warning(ncol(comm) - length(taxa), " taxa were dropped from comm because they were not in dis")
	if(length(taxa) < nrow(dis))
		warning(nrow(dis) - length(taxa), " taxa were dropped from dis because they were not in comm")
	dis <- dis[taxa, taxa]
	comm <- comm[,taxa]

	sim <- 1 - dis
	
	num <- matrix(nrow = nrow(comm), ncol=nrow(comm))
	for(a in 1:nrow(comm)) {
		for(b in a:nrow(comm)) {
			num[a,b] <- sum(outer(comm[a,], comm[b,]) * sim)
		}
	}
	num[lower.tri(num)] <- t(num)[lower.tri(num)]
	denom <- outer(0.5*diag(num), 0.5*diag(num), '+')
	res <- 1 - num/denom
	rownames(res) <- colnames(res) <- rownames(comm)
	res

}


#' Sorensen dissimilarity (without distance)
#' 
#' Compute original sorensen dissimilarity index, with no option for a distance matrix among species
#' 
#' @param comm Community matrix giving sites (rows) by species (columns).
#' @details The index is defined as 1 - (2 * (A ∩ B) / (A + B)). 
#' @references Sørensen, T. (1948) A method of establishing groups of equal amplitude in plant sociology 
#' based on similarity of species content. Kongelige Danske Videnskabernes Selskabs Biologiske Skrifter, 5, 1–34.
#' @return Matrix of pairwise dissimilarity values between communities.
sor <- function(comm)
{
	if(!is.matrix(comm)) comm <- as.matrix(comm)

	# convert to matrix of 0 and 1	
	if(any(! comm %in% c(0,1) )) comm <- 1 * (comm > 0)
	# get the intersection of all sites - number of sp in common
	# this is the 0.5 * numerator of sorensen similarity
	site_similarity <- comm %*% t(comm)
	
	# compute species richness and repeat into a matrix
	richness <- rowSums(comm)
	ra <- matrix(richness, nrow=length(richness), ncol=length(richness))
	
	# denominator - ra is the number of spp at site1, t(ra) is the number at site 2
	denom <- ra + t(ra)
	dimnames(denom) <- dimnames(site_similarity)
	
	1 - ((2 * site_similarity) / denom)
}
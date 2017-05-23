#' Compute bray-curtis dissimilarity based on species abundances
#' 
#' @param comm Community matrix giving sites (rows) by species (columns); values are abundances. Column and row labels are highly recommended
#' @return Site by site dissimilarity matrix
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


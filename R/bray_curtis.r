#' Compute bray-curtis dissimilarity based on species abundances
#' 
#' @param x a data frame. Minimally containing a column for site names, a column for species names present at
#'           each site, and a column for abundances
#' @param A optional site by species abundance matrix. Site names, if desired, can be specified as the row names. 
#'             If A is specified, then x is optional
#' @param sites column name or index of x containing site names
#' @param species column name or index of x containing species names
#' @param abundance column name or index of abundance column
#' @details If both x and A are provided, x will be used
#' @value Site by site dissimilarity matrix
#' @export
bc <- function(x, A=NULL, sites=1, species=2, abundance=3)
{
	if(missing(x)) {
		if(is.null(A)) stop("Either x or A must be specified")
	} else {
		if(!is.null(A)) warning("Both x and A were specified; using x")
		# compute abundance matrix from tall data
		if(is.numeric(sites)) sites = colnames(x)[sites]
		if(is.numeric(species)) species = colnames(x)[species]
		if(is.numeric(abundance)) abundance = colnames(x)[abundance]
		A <- reshape2::acast(x, as.formula(paste(sites, '~', species)), fill=0, value.var=abundance)
	}

	# get min abundance of shared species
	# pretty slow right now because it uses nested applies
	# a fast way to do this would be great
	shared_abundance <- apply(A, 1, function(site1) {
		apply(A, 1, function(site2) {
			mask <- site1 <= site2
			sum((site1 * mask) + ((!mask) * site2))
		})
	})
	
	# total abundance at each site for each pairwise combo
	site1 <- matrix(rowSums(A), nrow=nrow(A), ncol=nrow(A))
	site2 <- t(site1)
	
	1.0 - ((2.0 * shared_abundance)/(site1 + site2))
}


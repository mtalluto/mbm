#' Compute jaccard dissimilarity
#' 
#' @param x a data frame. Minimally containing a column for site names, a column for species names present at
#'           each site, and a column for abundances
#' @param community optional site by species presence-absence matrix. Site names, if desired, can be specified as the row names. 
#'             If community is specified, then x is optional
#' @param sites column name or index of x containing site names
#' @param species column name or index of x containing species names
#' @param abundance column name or index of abundance column
#' @details If both x and community are provided, x will be used
#' @return Site by site dissimilarity matrix
#' @export
jaccard <- function(x, community, site_col=1, species_col=2, pres_abs=NULL)
{
	if(missing(x)) {
		if(missing(community)) stop("Either x or community must be specified")
	} else {
		if(!missing(community)) warning("Both x and community were specified; using x")
		# compute presence absence matrix from tall data
		if(is.numeric(site_col)) site_col = colnames(x)[site_col]
		if(is.numeric(species_col)) species_col = colnames(x)[species_col]
		community <- reshape2::acast(x, as.formula(paste(site_col, '~', species_col)), fun.aggregate= function(a)
			length(a) > 0, fill=0, value.var=species_col)
	}
	
	# get the intersection of all sites - number of sp in common
	site_similarity <- community %*% t(community)
	
	# compute species richness
	richness <- rowSums(community)
	
	# make richness matrices for comparisons
	ra <- matrix(richness, nrow=length(richness), ncol=length(richness))
	
	1 - (site_similarity / (ra + t(ra) - site_similarity))
}

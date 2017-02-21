#' Mean pairwise distance
#' 
#' @param community Community matrix (required), giving sites (rows) by species (columns). Column names should match the
#'                   names in the phylogeny, if supplied.
#' @param phylogeny Phylogeny (optional) in the format used by \code{\link{ape}}; if supplied, phylogenetic diversity
#'                  will be computed.
#' @param traits Matrix of species (rows) by trait values (columns); traits may be a mix of quantitative
#'                 and categorical.
#' @param type What kind of diversity to compute; either alpha (a), beta (b) or both (ab)
#' @param ... Additional arguments to pass to \code{picante} functions.
#' @details This is a generic wrapper around several functions from the \code{\link{picante}} package. It has facilities
#'           for computing mpd for both alpha and beta diversity, for functional and phylogenetic diversity.
#' @return If a single index is selected a vector(for alpha) or matrix (for beta) with the selected diversity metric.
#'          If multiple indices are selected, a list of diversity vectors/matrices
#' @export
mpd <- function(community, phylogeny, traits, type=c('a', 'b', 'ab'), ...) {
	type <- match.arg(type)
	result <- list()
	if(!missing(phylogeny)) {
		phydist <- cophenetic(phylogeny)
		if(type == 'a' | type == 'ab')
			result$pd_alpha <- picante::mpd(community, phydist, ...)
		if(type == 'b' | type == 'ab')
			result$pd.beta <- picante::comdist(community, phydist, ...)
	}
	if(!missing(traits)) {
		stop("FD is not yet implemented")
	}
	if(length(result) == 1) result <- result[[1]]
	return(result)
}
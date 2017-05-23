#' Community weighted mean trait values
#' 
#' @param comm Site by species weights matrix; if presence-absence (1-0), the function returns the community mean
#' @param traits Matrix of traits; rownames will be matched to colnames in \code{comm}
#' @param na.rm Logical, whether to remove NAs when computing the mean
#' @return site by trait CWM matrix
#' @export
cwm <- function(comm, traits, na.rm=FALSE)
{
	# make sure all species in comm and in traits match
	taxa <- sort(intersect(colnames(comm), rownames(traits)))
	if(length(taxa) < ncol(comm)) 
		warning(ncol(comm) - length(taxa), " taxa were dropped from comm because they were not in traits")
	if(length(taxa) < nrow(traits))
		warning(nrow(traits) - length(taxa), " taxa were dropped from traits because they were not in comm")
	traits <- traits[taxa, taxa]
	comm <- comm[,taxa]
	
	apply(traits, 2, function(x) colMeans(t(comm) * x, na.rm=na.rm))
}

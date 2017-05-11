#' Mean pairwise distance
#' 
#' @param comm Community matrix giving sites (rows) by species (columns). Column names should match the
#'                   names in the phylogeny
#' @param phylogeny Phylogeny in the format used by \code{\link{ape}}; if missing dis is required
#' @param dis Species distance matrix; if missing will be extracted from the phylogeny
#' @param dis.transform Function (e.g., sqrt) to transform distance before computing mpd
#' @details Computes mean pairwise distance (MPD) for all site pairs (including within and between site MPD)
#' @return Matrix of pairwise distances; diagonal is witin site (alpha) MPD, other values are pairwise distances (beta MPD)
#' @useDynLib mbmtools c_mpd
#' @export
mpd <- function(comm, phylogeny, dis, dis.transform)
{
	if(missing(dis))
	{
		if(missing(phylogeny)) stop("Either phylogeny or dis matrix must be specified")
		dis <- cophenetic(phylogeny)
	} else {
		# do some error checking
		if(! (nrow(dis) == ncol(dis) & nrow(dis) == ncol(comm)))
			stop("dis must be square, with dims == ncol(comm)")
	}
	if(! missing(dis.transform))
		dis <- dis.transform(dis)
	
	# make sure columns line up
	if(any(rownames(dis) != colnames(dis)))
		dis <- order_symmetric_mat(dis)

	# drop any taxa from phylogeny that aren't present in community
	if(!all(rownames(dis) %in% colnames(comm)))
	{
		keep <- which(rownames(dis) %in% colnames(comm))
		warning(nrow(dis) - length(keep), " species were dropped from the phylogeny because they were not present in the community matrix")
		dis <- dis[keep,keep]
	}
	
	# drop taxa from community that are missing in the phylogeny
	if(!all(colnames(comm) %in% rownames(dis)))
	{
		keep <- which(colnames(comm) %in% rownames(dis))
		warning(ncol(comm) - length(keep), " species were dropped from the community because they were not present in the phylogeny")
		comm <- comm[,keep]
	}
		
	# make sure community and distance matrix labels line up
	if(any(colnames(comm) != colnames(dis)))
		comm <- match_mat_dims(comm, dis, by.x='c')
	
	spMat <- which(comm > 0, arr.ind = TRUE)
	spMat <- spMat[order(spMat[,1]),]
	sites = spMat[,1] - 1	# subtract one from indices because C is 0 indexed
	species = spMat[,2] - 1
	cDist <- numeric(nrow(comm)^2)
	res <- .C("c_mpd", species=as.integer(species), sites=as.integer(sites), N=as.integer(length(sites)), 
		dist=as.double(dis), K=as.integer(nrow(dis)), dest=as.double(cDist), S=as.integer(nrow(comm)))
	matrix(res$dest, nrow=nrow(comm), ncol=nrow(comm), dimnames=list(rownames(comm), rownames(comm)))
}


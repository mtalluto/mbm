#' Sorensen dissimilarity index
#' 
#' @param comm Community matrix giving sites (rows) by species (columns).
#' @param partition Logical; if true, the index will be partitioned into turnover and
#' 			nestedness.
#' @references Baselga A. 2010. Partitioning the turnover and nestedness components of 
#' 		beta diversity. Global Ecol. Biogeogr. 19: 134–143
#' @references Sørensen, T. (1948) A method of establishing groups of equal amplitude in 
#' 		plant sociology based on similarity of species content. Kongelige Danske 
#' 		Videnskabernes Selskabs Biologiske Skrifter, 5, 1–34.
#' @return If \code{partition == TRUE}, list of three matrices giving the turnover,
#' 		nestedness, and total components of dissimilarity. Otherwise, a matrix of
#' 		pairwise dissimilarity values between communities.
#' @export
sorensen <- function(comm, partition = FALSE)
{
	if(!is.matrix(comm)) comm <- as.matrix(comm)

	a <- comm %*% t(comm > 0)
	b <- comm %*% t(!comm)
	cc <- t(b)
	bsor <- (b + cc) / (2*a + b + cc)

	if(partition) {
		mbc <- pmin(b, cc)
		bturn <- mbc / (a + mbc)
		bnes <- bsor - bturn
		return(list(turnover = bturn, nestedness = bnes, sorensen = bsor))
	} else
		return(bsor)
}

#' Phylogenetic sorensen dissimilarity index
#' 
#' @param comm Community matrix giving sites (rows) by species (columns).
#' @param tree Dendrogram of class \code{phylo}
#' @param partition Logical; if true, the index will be partitioned into turnover and
#' 			nestedness.
#' @references Pavoine, S. and Ricotta, C. 2014. Functional and phylogenetic similarity 
#' 		among communities. Methods in Ecology and Evolution 5(7): 666–675.
#' @return If \code{partition == TRUE}, list of three matrices giving the turnover,
#' 		nestedness, and total components of dissimilarity. Otherwise, a matrix of
#' 		pairwise dissimilarity values between communities.
#' @export
phylosor <- function(comm, tree, partition = FALSE)
{
	if(is.null(tree$edge.length))
		stop("The tree must have branch lengths")

	if(any(! comm %in% c(0,1) )) {
		warning("Community matrix can only contain 0 and 1; assigning 1 to all nonzero 
			entries")
		comm <- 1 * (comm > 0)
	}

	# node labels and community columns must match
	commDrop <- which(!colnames(comm) %in% tree$tip.label)
	if(length(commDrop) > 0)
	{
		if(length(commDrop) == ncol(comm))
			stop("No columns in comm match tip labels in tree")

		warning("Dropping ", length(commDrop), " columns from comm that were 
			not found in tree")

		comm <- comm[,-commDrop]
	}
	treeDrop <- which(!tree$tip.label %in% colnames(comm))
	if(length(treeDrop) > 0)
	{
		if(length(treeDrop) == length(tree$tip.labe))
			stop("No tip labels in tree were found in column labels in comm")

		warning("Dropping ", length(treeDrop), " modes from tree that were 
			not found in comm")

		tree <- ape::drop.tip(tree, treeDrop)
	}

	if(any(rowSums(comm) == 0))
	{
		drop <- which(rowSums(comm) == 0)
		warning(length(drop), " sites had 0 species and were dropped")
		comm <- comm[-drop,]
	}

	branch_sp <- branch_adjacency(tree, TRUE)
	# save weights
	wts <- attr(branch_sp, "branch.lengths")
	# match the margins
	branch_sp <- branch_sp[,match(colnames(comm), colnames(branch_sp))]

	site_branch <- comm %*% t(branch_sp)
	# convert site by branch to presence absence (instead of sp count)
	site_branch <- 1 * (site_branch > 0)
	# now weight by branch lengths
	site_branch <- sweep(site_branch, 2, wts, `*`)
	sorensen(site_branch, partition)
}


#' Sorensen dissimilarity index based on a dendrogram
#' 
#' Produces dissimilarity matrices using a species by species distance matrix; useful
#' for example for functional diversity.
#' 
#' @param comm Community matrix giving sites (rows) by species (columns).
#' @param distmat Distance matrix; margins must match exactly the columns in comm 
#' @param partition Logical; if true, the index will be partitioned into turnover and
#' 			nestedness.
#' @param ... Additional arguments to pass to \code{\link{hclust}}.
#' @references Pavoine, S. and Ricotta, C. 2014. Functional and phylogenetic similarity 
#' 		among communities. Methods in Ecology and Evolution 5(7): 666–675.
#' @return If \code{partition == TRUE}, list of three matrices giving the turnover,
#' 		nestedness, and total components of dissimilarity. Otherwise, a matrix of
#' 		pairwise dissimilarity values between communities.
#' @export
dendrosor <- function(comm, distmat, partition = TRUE, ...)
{
	if(!inherits(distmat, "dist"))
		distmat <- as.dist(distmat)

	if(!attr(distmat, "Size") == ncol(comm) & 
		any(attr(distmat, "Labels") != colnames(comm)))
		stop("ncol(comm) must equal the number of rows/columns in distmat and the
			margin labels must be identical")

	tree <- hclust(distmat, ...)
	phylosor(comm, ape::as.phylo(tree), partition)
}

#' Branch adjacency
#' 
#' Produce a branch adjacency matrix (i.e., branch by node)
#' 
#' @param x Object of class phylo
#' @param branch.lengths Logical; if true, the values of the adjacency matrix will be the
#' 		branch lengths, otherwise they will be 1 or 0
#' @return A branch by node adjacency matrix
#' @keywords internal
branch_adjacency <- function(x, branch.lengths = !is.null(x$edge.length)) {
    m <- matrix(0, ncol=length(x$tip), nrow=nrow(x$edge))
    from <- x$edge[,1]
    to <- x$edge[,2]
    g <- seq_along(x$tip.label)
    while (any(!is.na(g))) {
        i <- match(g, to)
        m[cbind(i, seq_along(i))] <- 1
        g <- from[i]
    }
    rownames(m) <- paste0("Branch", seq.int(nrow(m)))
    colnames(m) <- x$tip.label
    if(branch.lengths)
    {
    	attr(m, "branch.lengths") <- x$edge.length
    }
    return(m)
}


#' Sorensen dissimilarity (without partitioning)
#' 
#' Compute sorensen dissimilarity index; deprecated method (but faster than new one)
#' 
#' @param comm Community matrix giving sites (rows) by species (columns).
#' @details The index is defined as 1 - (2 * (A ∩ B) / (A + B)). 
#' @return Matrix of pairwise dissimilarity values between communities.
#' @keywords internal
sor <- function(comm)
{
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
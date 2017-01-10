#' Sorensen phylogenetic distance
#' 
#' @param community Community matrix giving sites (rows) by species (columns). Column names should match the
#'                   names in the phylogeny
#' @param phylogeny Phylogeny in the format used by \link{\code{ape}}
#' @param na.species.action Action to take if species are found in the community matrix that are missing in the
#'         phylogeny. See details.
#' @details Compute an analogue of the sorensen distance, using shared pd instead of taxonomic diversity
#' @details If na.species.action is set to 'omit,' missing species will be dropped from the community matrix
#' @return Matrix of pairwise dissimilarity values between communities
#' @export
sorensen_pd <- function(communities, phylogeny, na.species.action = c('error', 'omit')) {
	na.species.action <- match.arg(na.species.action)
	
	# remove unneeded tips from the phylogeny
	drop.phylogeny <- phylogeny$tip.label[ which(! phylogeny$tip.label %in% colnames(communities))] 
	if(length(drop.phylogeny) > 0) phylogeny <- ape::drop.tip(phylogeny, drop.phylogeny)

	# check for missing species in the phylogeny and drop if necessary
	drop.communities <- which(! colnames(communities) %in% phylogeny$tip.label)
	if(na.species.action == 'error' & length(drop.communities) > 0) {
		stop("Missing species found in community matrix; set na.species.action to 'omit' to ignore")
	 } else if(length(drop.communities) > 0) communities <- communities[,-drop.communities]

	# set up data frames of branches
	branches <- cbind(phylogeny$edge, data.frame(phylogeny$edge.length))
	ext.branches <- branches[ which(branches[,2] <= length(phylogeny$tip.label) ),]
	labs <- data.frame(phylogeny$tip.label)
	ext.branches.n <- merge(ext.branches, labs, by.x='2', by.y=0)
	int.branches <- branches[which (!branches[ ,2] %in% 1:length(phylogeny$tip.label)), ]

	#####Make species by node matrix
	dimnames <- list(as.character(ext.branches.n$"phylogeny.tip.label"), as.character(int.branches$"2")) #changed from '2'
	tip.branch.mat <- matrix(0,nrow(ext.branches.n), nrow(int.branches), dimnames=dimnames) 

	## make site by branch matrix
	idx.first.branch <- min(int.branches$"2")
	for (brch in int.branches$"2") {
	  tips.nums <- picante::internal2tips(phylogeny, brch)
	  tip.branch.mat[tips.nums, brch-idx.first.branch+1] <- 1
	  ###reorder tip.branch to match the site.species  
	  tip.branch.mat.order <- tip.branch.mat[match(colnames(communities),rownames(tip.branch.mat)),]
	  site.branch <- communities %*% tip.branch.mat.order
	  site.branch.ones <- (site.branch&TRUE)*1
	}
	
	#make a three column matrix of branches and branch lengths
	weights.mat=c(ext.branches.n[,3],int.branches[,3])
  
	#bin tips and internal
	communities[communities > 1] <- 1
	all <- cbind(communities[order(rownames(communities)),], site.branch.ones[order(rownames(communities)),])
	all.weight <- weights.mat * t(all)
  
	similarity <- t(all.weight) %*% (all.weight)
	self.similiarity <- diag(similarity)
	denom <- sapply(self.similiarity, function(s) s + self.similiarity)
	sorensen <- 1 - (2*similarity / denom)
	
	return(sorensen)

}
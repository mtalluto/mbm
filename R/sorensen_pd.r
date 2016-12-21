#' Sorensen phylogenetic distance
#' 
#' @param community Community matrix giving sites (rows) by species (columns). Column names should match the
#'                   names in the phylogeny
#' @param phylogeny Phylogeny in the format used by \link{\code{ape}}
#' @param na.species.action Action to take if species are found in the community matrix that are missing in the
#'         phylogeny. See details.
#' @details Compute an analogue of the sorensen distance, using shared pd instead of taxonomic diversity
#' @details If na.species.action is set to 'omit,' missing species will be dropped from the community matrix
#' @value Matrix of pairwise dissimilarity values between sites
#' @export
sorensen_pd <- function(communities, phylogeny, na.species.action = c('error', 'omit')) {
	na.species.action <- match.arg(na.species.action)
	
	# remove unneeded tips from the phylogeny
	drop.tree <- tree$tip.label[ which(! tree$tip.label %in% colnames(sites))] 
	tree <- drop.tip(phylogeny, drop.tree)

	# check for missing species in the phylogeny and drop if necessary
	drop.sites <- which(! colnames(sites) %in% tree$tip.label)
	if(na.species.action == 'error' & length(drop.sites) > 0)
		stop("Missing species found in community matrix; set na.species.action to 'omit' to ignore")
	sites <- sites[,-drop.sites]

	# set up data frames of branches
	branches <- cbind(tree$edge, data.frame(tree$edge.length))
	ext.branches <- branches[ which(branches[,2] <= length(tree$tip.label) ),]
	labs <- data.frame(tree$tip.label)
	ext.branches.n <- merge(ext.branches, labs, by.x='2', by.y=0)
	int.branches <- branches[which (!branches[ ,2] %in% 1:length(tree$tip.label)), ]

	#####Make species by node matrix
	dimnames <- list(as.character(ext.branches.n$"tree.tip.label"), as.character(int.branches$"2")) #changed from '2'
	tip.branch.mat <- matrix(0,nrow(ext.branches.n), nrow(int.branches), dimnames=dimnames) 

	## make site by branch matrix
	idx.first.branch <- min(int.branches$"2")
	for (brch in int.branches$"2") {
	  tips.nums <- internal2tips(tree, brch)
	  tip.branch.mat[tips.nums, brch-idx.first.branch+1] <- 1
	  ###reorder tip.branch to match the site.species  
	  tip.branch.mat.order <- tip.branch.mat[match(colnames(sites),rownames(tip.branch.mat)),]
	  site.branch <- sites %*% tip.branch.mat.order
	  site.branch.ones <- (site.branch&TRUE)*1
	}
	
	#make a three column matrix of branches and branch lengths
	weights.mat=c(ext.branches.n[,3],int.branches[,3])
  
	#bin tips and internal
	sites[sites > 1] <- 1
	all <- cbind(sites[order(rownames(sites)),], site.branch.ones[order(rownames(sites)),])
	all.weight <- weights.mat * t(all)
  
	similarity <- t(all.weight) %*% (all.weight)
	self.similiarity <- diag(similarity)
	denom <- sapply(self.similiarity, function(s) s + self.similiarity)
	sorensen <- 1 - (2*similarity / denom)
	
	return(sorensen)

}
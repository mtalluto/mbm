#' Build a functional dendrogram to be used with tree-based diversity indices
#' 
#' @param x "Tall" dataframe with species and trait information; see 'details'
#' @param vartypes 2-column data frame or character matrix giving (in the first column) the trait names and (in the
#'         second column) the trait types. Values in the second column should be one of 'Q' (quantitative), 'N' (nominal
#'         or categorical), or 'O' (ordinal).
#' @param species column number or name in x containing species names, or a vector of species names
#' @param traits column number or name in x of trait names, or a vector of trait names
#' @param values column or name of trait values, or vector of trait values; see 'details'
#' @details If 'x' is missing, then species, traits, and values must be vectors containing the data
#' @return named list containing a species by species trait distance matrix and a functional dendrogram (in phylo format)
#' @export
dendrogram <- function(x, vartypes, species=1, traits=2, values=3)
{
	fun.aggregate = list(Q=mean, N=function(x) x[1], O=function(x) x[1])
	
	if(missing(x)) {
		x <- data.frame(species=species, traits=traits, values=values)
	} else {
		colnames(x)[c(species, traits, values)] <- c('species', 'traits', 'values')
	}
	spList <- sort(unique(x[,species, drop=FALSE]))
	
	ktab <- ade4::ktab.list.df(lapply(unique(varypes[,2]), function(type) {
		trNames <- varypes[vartypes[,2] == type,1]
		rows <- which(traits %in% trNames)
		trdf <- reshape2::dcast(x[rows,], species ~ traits, value.var='values', fun.aggregate=fun.aggregate[[type]])
		trdf <- merge(spList, trdf, all=TRUE)
		rownames(trdf) <- trdf[,1]
		trdf <- trdf[,-1]
		for(i in 1:ncol(trdf)) trdf[is.nan(trdf[,i]),i] <- NA
		trdf
	}))
	
	tr_dist <- ade4::dist.ktab(guis_ktab, unique(varypes[,2]))
	tmpf <- tempfile(fileext='.nex')
	phlo <- as.phylo(hclust(dist(tr_dist)))
	# write the tree to a temp file and read back in to fix indexing, per LJP
	ape::write.nexus(phlo, file=tmpf)
	phlo <- ape::read.nexus(tmpf)
	list(distance=tr_dist, dendrogram=phlo)
}


#' Process functional trait data, combining multiple protocols and multiple ways of coding
#' 
#' @param x "Tall" dataframe with species and trait information
#' @param traits column number or name in x of trait names
#' @param categories column number or name for trait categories (i.e., the 'value' of text-coded variables)
#' @param values column number or name for trait values (i.e., numeric values)
#' @return dataframe similar to x, but with re-coded traits and with numeric values for text-coded traits
#' @export
combine_traits <- function(x, traits=1, categories=2, values=3)
{
	## combine traits with different protocols
	x[,traits][x[,traits] == "LDMCp" | x[,traits] == "LDMCst"] <- "LDMC"
	x[,traits][x[,traits] == "LNCst" | x[,traits] == "LNCp"] <- "LNC"
	x[,traits][x[,traits] == "PL_REPR_H (ARVES)"] <- "PL_REPR_H"
	x[,traits][x[,traits] == "PL_VEG_H (ARVES)" | x[,traits] == "PL_VEG_H_mean"] <- "PL_VEG_H"
	x[,traits][x[,traits] == "SLAp" | x[,traits] == "SLAst"] <- "SLA"
	x[,traits][x[,traits] == "LHIST_INDK"] <- "LHIST"

	# combine categories for some categorical variables
	x[,categories][x[,categories] == "Chamaephyte" | x[,categories] == "woody chamaephyte"] <- "chamaephyte"
	x[,categories][x[,categories] == "Geophyte"] <- "geophyte"
	x[,categories][x[,categories] == "Hemicryptophyte" | x[,categories] == "hemicryptophyte (long-lived)"] <- 
								"hemicryptophyte"
	x[,categories][x[,categories] == "Phanerophyte"] <- "phanerophyte"
	x[,categories][x[,categories] == "Therophyte"] <- "therophyte"
	
	# numerical coding for categorical variables
	cat_traits <-c("DISP_INDK", "LHIST")
	for(v in cat_traits) {
		rows <- which(x[,traits] == v)
		x[rows,values] <- as.integer(factor(x[rows,categories]))
	}
	
	# numerical coding for ordinal text-based traits
	ord_traits <- c("ELLLITE_INDK", "ELLMOIST_INDK", "ELLNUT_INDK")
	rows <- which(x[,traits] %in% ord_traits)
	pattern <- ".+\\(([0-9.]+)/[0-9]\\)"
	x[rows,values] <- as.numeric(sub(pattern, "\\1", x[rows,categories]))
	
	x
}
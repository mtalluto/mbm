#' Higher-level function for reading predictions
#' 
#' Reads all predictions and site names if available
#' @param fname The name of the input file used to generate the predictions
#' @param nsamp Number of samples, if used
#' @param tfOutput file extension for output files
#' @param namesExt What should be appended to the names filename
#' @keywords internal
#' @return A data frame of predictions and site names
read_mbm_predict <- function(fname, nsamp=NA, tfOutput = '_out.csv', nameExt = '.names')
{
	stop("need to check for bigpr")
	file <- paste0(fname, tfOutput)
	colnames <- if(is.na(nsamp)) c('fit', 'stdev') else paste0('samp', 1:nsamp)
	preds <- data.table::fread(file, sep=',', data.table=FALSE, col.names = colnames)
	nmFile <- paste0(fname, nameExt)
	if(file.exists(nmFile))
	{
		prNames <- data.table::fread(nmFile, sep=',', data.table=FALSE, col.names = c('site1', 'site2'))
		preds <- cbind(prNames, preds)
	}
	preds
}


#' Write prediction data to disk
#' @param x Prediction dataset, including site names in the last two columns
#' @param datname The name of the prediction set
#' @param bigLim Integer, what's the largest number of sites before we start treating this as a big data problem
#' @param namesExt What should be appended to the names filename
#' @keywords internal
#' @return a character vector; the first element is a string for the file name, the second is an argument for passing
#'      to the mbm python program; as a side effect, the file is written; a companion file is written with just site names
write_mbm_predict <- function(x, datname, tfBase = 'mbm_', tfExt = '.csv', bigLim = 200, namesExt = '.names')
{
	if(missing(datname)) 
		datname <- ""

	write_pr <- function(x, fname)
	{
		sitecols <- grep('site', colnames(x))
		if(length(sitecols) > 0)
		{
			nmfile <- paste0(fname, namesExt)
			data.table::fwrite(as.data.frame(x[,-sitecols]), fname)
			data.table::fwrite(as.data.frame(x[,sitecols]), nmfile)
		} else {
			data.table::fwrite(as.data.frame(x), fname)
		}
	}

	if(nrow(x) > bigLim^2)
	{
		file <- tempfile(paste0(tfBase, 'bigpr_', datname))
		arg <- paste0('--bigpr=', file)
		if(!dir.create(file)) 
			stop("Could not create a tempdir for large prediction dataset")
		inds <- seq(1, nrow(x), bigLim^2)
		if(inds[length(inds)] <= nrow(x))
			inds <- c(inds, nrow(x) + 1)
		for(i in 1:(length(inds) - 1))
		{
			st <- inds[i]
			en <- inds[i+1] - 1
			dat <- x[st:en, , drop=FALSE]
			filename <- file.path(file, paste0('bigpr', i, tfExt))
			write_pr(dat, filename)
		}
	} else {
		file <- tempfile(paste0(tfBase, 'pr_', datname, '_'), fileext=tfExt)
		arg <- paste0('--pr=', file)
		write_pr(x, file)
	}
	return(c(file, arg))
}


#' Write MBM data to disk
#' 
#' @param x An MBM object
#' @param tfBase Base name for the temp files
#' @param tfExt Extension for temp files
#' @param namesExt What should be appended to the names filename
#' @keywords internal
#' @return A named character vector with the file names for the \code{response},
#'     \code{covariates}, and \code{params}. As a side effect, the response,
#'     covariate, and (if present) params are written to the files.
write_mbm_dat <- function(x, tfBase = 'mbm_', tfExt = '.csv', namesExt = '.names')
{
	files <- c(response = tempfile(paste0(tfBase, 'y_'), fileext=tfExt),
		    covariates = tempfile(paste0(tfBase, 'x'), fileext=tfExt),
		    params = tempfile(paste0(tfBase, 'par_'), fileext = tfExt))
	data.table::fwrite(as.data.frame(x$response), files['response'])
	data.table::fwrite(x$covariates, files['covariates'])
	data.table::fwrite(x$covar_sites, paste0(files['covariates'], namesExt))
	if("params" %in% names(x))
		data.table::fwrite(as.data.frame(x$params), files['params'])
	return(files)

}
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
#' @return a string for the file name; as a side effect, the file is written; a companion file is written with just site names
write_mbm_predict <- function(x, datname, tfBase = 'mbm_', tfExt = '.csv', bigLim = 200, namesExt = '.names')
{
	if(nrow(x) >= bigLim^2)
		warning("Large (>", bigLim, " sites) dataset for prediction; you may experience memory issues")
	if(missing(datname)) 
		datname <- ""
	file <- tempfile(paste0(tfBase, 'pr_', datname, '_'), fileext=tfExt)
	nmfile <- paste0(file, namesExt)
	sitecols <- grep('site', colnames(x))
	if(length(sitecols) > 0)
	{
		data.table::fwrite(as.data.frame(x[,-sitecols]), file)
		data.table::fwrite(as.data.frame(x[,sitecols]), nmfile)
	} else {
		data.table::fwrite(as.data.frame(x), file)
	}
	return(file)
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
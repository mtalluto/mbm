#' Create an mbm model
#' 
#' @param y Square dissimilarity or distance matrix, can be complete or lower triangular only. Row and column names are
#'          required and must match the site names in the rows of \code{x}.
#' @param x Matrix giving a series of covariates (in columns) for all sites (in rows). Row names are required. All 
#'          variables will be included in the model.
#' @param predictX List of prediction datasets; each list element is a matrix with same number of columns as \code{x}. See details.
#' @param scale Boolean, if true, x values will be centered and scaled before fitting the model.
#' @param y_name A name to give to the y variable
#' @param GPy_location Optional character giving the location of the user's GPy installaion
#' @param pyCmd Where to look for python; the version in use must have GPy installed
#' 
#' @details Prediction datasets can either be supplied when the model is called, or by using the \code{predict} method on the \code{mbm}
#'          object. The former will generally be much faster to run; see \code{\link{predict.mbm}}. Note that predictions are always
#'          generated for the input dataset.
#' @return An S3 object of class mbm. 
#' @export
mbm <- function(y, x, predictX, scale = TRUE, y_name = 'beta', GPy_location, pyCmd = 'python')
{
	# y <- as.matrix(y)
	# x <- as.matrix(x)
	if(scale) 
		x <- scale(x)
	xDF <- env_dissim(x)
	
	yDF <- reshape2::melt(y,varnames=c('site1', 'site2'), value.name = y_name)
	dat <- merge(xDF, yDF, all.x = TRUE)
	x_cols <- which(!colnames(dat) %in% c('site1', 'site2', y_name))
	
	tfBase <- "mbm_"
	tfExt <- ".csv"
	tfOutput <- '_out.csv'
	
	# generate temporary files
	yFile <- tempfile(paste0(tfBase, 'y_'), fileext=tfExt)
	xFile <- tempfile(paste0(tfBase, 'x'), fileext=tfExt)
	data.table::fwrite(as.data.frame(dat[,y_name]), yFile)
	data.table::fwrite(dat[,x_cols], xFile)

	
	# deal with prediction datasets
	if(!missing(predictX))
	{
		stop("Prediction not implemented yet")
		if(!is.list(predictX))
			predictX <- list(predictX)
		if(scale) {
			warning("Prediction datasets will be scaled to the same scale as x")
		}
		predictFiles <- sapply(names(predictX), function(nm) tempfile(paste0(tfBase, 'pr_', nm, '_'), fileext=tfExt))
		
	}

	# run the model
	mbmArgs <- c(system.file('mbm.py', package='mbmtools', mustWork = TRUE), # the file name of the python script
							 paste0('--y=', yFile), paste0(' --x=', xFile), paste0(' --out=', tfOutput))  # additional arguments
	if(!missing(GPy_location))
		mbmArgs <- c(mbmArgs, paste0('--gpy=', GPy_location))
	
	result <- system2('python', args=mbmArgs)
	if(result) 
		stop("MBM returned an error")
		
	# collect results
	fits <- data.table::fread(paste0(xFile, tfOutput), sep=',', data.table=FALSE)
	
	model <- list(
		response = dat[,y_name],
		covariates = dat[,x_cols],
		fits = fits
		## NEED TO SAVE X scaling!
	)
	class(model) <- c('mbm', class(model))
	model
}
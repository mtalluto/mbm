#' MBM Fit Class
#' 
#' Creates and object from an mbm fit (from Python)
#' 
#' @param fit_file File name for the data used to fit the model
#' @param result_dir Directory with results
#' @param holdout_file (optional) File name for the holdout dataset
#' @param predict_file (optional) File name for the prediction dataset
#' @param link_function Link function used to fit the model
#' @param model_name What to call the model
#' @param type Alpha or beta diversity?
#' @return An object of class mbmfit
#' @export
mbmfit <- function(fit_file, result_dir, holdout_file, predict_file, link_function = c("identity", "log", "probit"), 
									 model_name = basename(result_dir), type=c("alpha", "beta"))
{
	cat("1\n")
	fit.colnames <- c('fit.mean', 'fit.sd', 'fit.lower', 'fit.upper')
	link_function <- match.arg(link_function)
	type <- match.arg(type)
	fitDat <- read.csv(fit_file)
	resultFiles <- list.files(result_dir)
	fit_res_file <- resultFiles[grep("fit", resultFiles)]
	fitResult <- read.table(file.path(result_dir, fit_res_file), header=FALSE, sep=' ')
	colnames(fitResult) <- fit.colnames
	cat("2\n")
	
	val <- list(fits = cbind(fitDat, fitResult))

	cat("3\n")
	if(missing(holdout_file))
	{
		val$holdout <- NULL
	} else {
		hd <- read.csv(holdout_file)
		hold_res_file <- resultFiles[grep("holdout", resultFiles)]
		hf <- read.table(file.path(result_dir, hold_res_file), header=FALSE, sep=' ')
		colnames(hf) <- fit.colnames
		val$holdout <- cbind(hd, hf)
	}
	cat("4\n")
	
	if(missing(predict_file))
	{
		val$predict <- NULL
	} else {
		pd <- read.csv(predict_file)
		pred_res_file <- resultFiles[grep("newX", resultFiles)]
		pf <- read.table(file.path(result_dir, pred_res_file), header=FALSE, sep=' ')
		colnames(pf) <- fit.colnames
		val$predict <- cbind(pd, pf)
	}

	attr(val, "name") <- model_name
	attr(val, "link") <- link_function
	attr(val, "type") <- type
	class(val) <- c(class(val), "mbmfit")
	val
}
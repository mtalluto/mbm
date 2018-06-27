#' Check if necessary Python packages can be loaded
#'
check_python <- function()
{
	efun <- function(module, err)
		paste("Loading python module", module, "produced an error:\n", err)

	modules <- c('GPy', 'numpy', 'scipy')
	mTest <- sapply(modules, function(m) {
		tryCatch(sp <- reticulate::import(m), warning = function(w) warning(w), 
			error = function(e) warning(efun(m, e)))
	})
}

#' Convenience function to install python dependencies
#'
#' @details MBM includes requires python dependencies. Assuming a working Python installation,
#' this function will install them for you
#' NOT PRESENTLY IMPLEMENTED; to-do
#' @internal
install_python_dependencies <- function()
{
	stop("Feature not implemented")
}
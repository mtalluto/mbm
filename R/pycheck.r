#' Check if necessary Python packages can be loaded
#' @examples check_python()
#' @export
check_python <- function()
{
	efun <- function(module, err)
		paste("Loading python module", module, "produced an error:\n", err)

	modules <- c('GPy', 'numpy', 'scipy')
	mTest <- sapply(modules, function(m) {
		tryCatch(sp <- reticulate::import(m), warning = function(w) warning(w), 
			error = function(e) warning(efun(m, e)))
	})

	print("check_python done; if no errors/warnings you are good to go!")
}


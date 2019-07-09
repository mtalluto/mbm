.onLoad <- function(libname, pkgname){
	# try to load GPy; if it fails, package loading will fail with an error
	GPy <- reticulate::import("GPy")
}

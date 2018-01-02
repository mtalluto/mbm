#' Predict method for MBM objects
#' 
#' @param x A previously-fit MBM object
#' @param newdata Optional dataset for prediction. If present, it should be either a character vector giving
#'    the name of one of the datasets specified for \code{predictX} when the model was fit, or  a new dataset
#'    in the same format used to fit the model (i.e., a site by covariate matrix). If missing, predictions 
#'    will be for the original data.
#' @param n_samples NA or integer; if NA, analytical predictions with standard deviation are returned, otherwise posterior samples are returned.
#'     Currently ignored if rasterdat is specified.
#' @param GPy_location Optional string giving the location of the user's GPy installaion
#' @param pyMsg boolean, should we print messages from python? Useful for debugging
#' @details Prediction to new data is possible after the fact for mbm models, however there are significant
#'     performance penalties for doing so. Thus, whenever possible, it is preferable to predict during
#'     model fitting via the \code{predictX} argument to the \link{\code{mbm}} function. All prediction
#'     is done on the response scale.
#'     This function caches to disk, thus it is important to ensure that adequate disk space is
#'     available when using large prediction datasets.
#' @return A data frame of predictions and standard deviations (on the link scale); use \code{x$y_rev_transform(x$rev_link(predictions$fit))}
#'    for the response scale
#' @export
predict.mbm <- function(x, newdata, n_samples = NA, GPy_location = NA, pyMsg = FALSE)
{
	tfOutput <- '_out.csv'


	if(missing(newdata))
	{
		preds <- x$fitted.values
	} else if(is.character(newdata))
	{
		preds <- x$predictions[[newdata]]
	} else {

		# 1. parse x to recreate an mbm call similar to the mbm function
		files <- write_mbm_dat(x)
		# 2. set up mbm call to be a re-launch, not a new model
		mbmArgs <- make_args(x, files, GPy_location=GPy_location, n_samples = n_samples)
		mbmArgs <- c(mbmArgs, '--resume')

		# 3. parse newdata into predictX format
		dat <- prep_predict(newdata, x)
		prFile <- write_mbm_predict(dat)
		mbmArgs <- c(mbmArgs, prFile[2])
		# mbmArgs <- c(mbmArgs, sapply(prFile, function(prf) paste0('--pr=', prf)))

		# 4. run model
		result <- system2('python', args=mbmArgs, stdout = TRUE)
		if(pyMsg) print(result)
		if("status" %in% names(attributes(result))) 
			stop("MBM returned an error: ", attr(result, "status"))

		# 5. post-process
		preds <- read_mbm_predict(prFile[1], nsamp = n_samples)
	}

	return(preds)
}


#' Spatial MBM prediction
#' 
#' @param x A previously-fit MBM object
#' @param rasterdat Raster stack containing named layers matching the variable names in x (i.e., colnames(x$covariates)[-1]).
#'      If a layer named 'names' is included, this layer will be used as sitenames, otherwise they will be assigned unique
#'      numbers
#' @param method How to compute the spatial predictions; see details
#' @param ... Other named parameters to pass to \code{\link{predict.mbm}}.
#' @details If \code{method} is "slow", spatial predictions will be computed by first predicting dissimilarity to all 
#'      pairs of raster cells, then performing an ordination on the dissimilarity matrix to produce an RGB raster
#'      of spatial predictions.
#'      For method == 'fast' (currently not implemented), computation is sped up by first performing hierarchical clustering
#'      on the predicted dissimilarity matrix for the calibration data (which will have already been computed when mbm was run)
#'      to produce cell categories. Each raster cell will then be assigned the category of the calibration data point that is
#'      closest environmentally. Then, we compute the dissimilarity matrix of the categories (based on the mean environmental
#'      values). The ordination is performed as with the slow method on this dissimilarity matrix.
#' @return An object of class mbmSP, which is a list with three named items: fits is a 3-band raster giving the first three
#'     prinipal components of predicted pairwise dissimilarity, stdev is a raster giving the mean of pairwise dissimilarities
#'     among all other sites in a given site, and pca is the principal components analysis for the fits.
#' @export
spatial_predict <- function(x, rasterdat, method = c('slow', 'fast'), ...)
{
	method <- match.arg(method)

	if(method == 'fast') {
		stop('Fast method is not implemented, use method="slow"')
	} else {
		preds <- predict_mbm_raster(x, rasterdat)
		fits <- x$y_rev_transform(x$rev_link(preds$fits))
		# for fits, use a PCoA to collapse to 3 axes. For stdev, use rowmeans to collapse to one
		fitPCoA <- ade4::dudi.pco(as.dist(fits), scannf = FALSE, nf = 3)
		fitPCoA_scaled <- as.data.frame(apply(as.matrix(fitPCoA$l1), 2, function(x) {
			x <- x - min(x)
			x / max(x)
		}))
		sdMeans <- data.frame(sd = rowMeans(preds$stdev, na.rm = TRUE))

		# make rasters
		sp::coordinates(fitPCoA_scaled) <- sp::coordinates(sdMeans) <- preds$coords
		sp::gridded(fitPCoA_scaled) <- sp::gridded(sdMeans) <- TRUE
		pcaRas <- raster::stack(fitPCoA_scaled)
		sdRas <- raster::raster(sdMeans)
		ret <- list(fits = pcaRas, stdev = sdRas, pca = fitPCA)
		class(ret) <- c("mbmSP", class(ret))
	}
	ret
}





#' Prediction for MBM from a raster dataset
#' 
#' @param x A previously-fit MBM object
#' @param rasterdat Raster stack containing named layers matching the variable names in x (i.e., colnames(x$covariates)[-1]).
#'      If a layer named 'names' is included, this layer will be used as sitenames, otherwise they will be assigned unique
#'      numbers
#' @param ... Other named parameters to pass to \code{\link{predict.mbm}}.
#' @return A named list; \code{fits} is a cell by cell matrix of predictions (on the link scale; use \code{x$y_rev_transform(x$rev_link(predictions$fit))} 
#'    for the response scale), \code{stdev} is a cell by cell matrix of standard deviations, and \code{coords} is a matrix of coordinates. Row/column
#'    names in \code{fits} and \code{stdev} match the rownames in \code{coords}.
#' @export
predict_mbm_raster <- function(x, rasterdat, ...)
{
	newdata <- raster::getValues(rasterdat[[colnames(x$covariates)[-1]]]) # the -1 is to account for the fact that the first covariate name is always distance
	rows <- complete.cases(newdata)
	newdata <- newdata[rows,]
	coords <- sp::coordinates(rasterdat)[rows,]
	if("names" %in% names(rasterdat))
	{
		names <- raster::getValues(rasterdat[['names']])
		names <- names[rows]
	} else {
		names <- 1:nrow(newdata)
	}
	rownames(newdata) <- rownames(coords) <- names
	preds <- predict(x, newdata, ...)
	diagSites <- unique(c(preds$site1, preds$site2))
	predsDF <- rbind(preds, data.frame(site1 = diagSites, site2 = diagSites, fit = 0, stdev = NA))

	# make site by site matrices and fill in lower triangle
	fitMat <- make_symmetric(predsDF, site1 ~ site2, value.var = "fit")
	sdMat <- make_symmetric(predsDF, site1 ~ site2, value.var = "stdev")
	list(fits = fitMat, stdev = sdMat, coords = coords)
}


#' Turn an MBM prediction dataframe into a symmetric matrix 
#' @param DF MBM predfiction dataframe
#' @param formula A formula to be passed to \code{link{reshape2::acast}}
#' @param value.var Name of value variable for \code{acast}
#' @param ... Additional parameters for \code{acast}
#' @return A symmetric matrix of predictions
#' @keywords internal
make_symmetric <- function(DF, formula, value.var, ...)
{
	mat <- reshape2::acast(DF, formula, value.var = value.var, ...)
	mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
	mat
}
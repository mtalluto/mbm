#' Predict method for MBM objects
#' 
#' @param x A previously-fit MBM object
#' @param newdata Optional dataset for prediction. If present, it should be a new dataset in 
#' 		the same format used to fit the model (i.e., a site by covariate matrix). If missing,
#' 		predictions will be for the original data.
#' @param n_samples NA or integer; if NA, analytical predictions with standard deviation 
#' 		are returned, otherwise posterior samples are returned.
#' @param type Whether to return predictions on the link or response scale
#' @details Prediction to new data is possible after the fact for mbm models, however 
#' 		there can be performance penalties for doing so with large models. Thus, it is 
#'		sometimes preferable to predict during model fitting via the \code{predictX} 
#' 		argument to the \code{\link{mbm}} function. 
#' 
#' 		All prediction is done on the link scale.
#' 
#'		This function caches to disk, thus it is important to ensure that adequate disk 
#' 		space is available when using large prediction datasets.
#' @return A data frame of predictions and standard deviations (on the link scale); use 
#' 		\code{x$y_rev_transform(x$rev_link(predictions$fit))} for the response scale.
#' @export
# predict.mbm <- function(x, newdata, n_samples = NA, GPy_location = NA, pyMsg = FALSE)
predict.mbm <- function(x, newdata, type=c("link", "response")) {
	type <- match.arg(type)
	if(missing(newdata)) {
		newdata <- x$covariates
	} else {
		newdata <- as.matrix(newdata)
		if(ncol(newdata) != ncol(x$x)) {
			stop("newdata must have the same number of variables as the original data")
		}

		# parse newdata into dissimilarity format
		newdata <- x$x_scaling(newdata)
		newdata <- env_dissim(newdata)
		newdata <- as.matrix(newdata[,-which(colnames(newdata) %in% c("site1", "site2"))])
	}

	pr <- x$pyobj$gp$predict_noiseless(newdata)
	pr[[2]] <- sqrt(pr[[2]])
	names(pr) <- c("mean", "sd")
	if(type == "link") {
		as.data.frame(pr)		
	} else {
		pr <- pr[[1]]
		pr <- x$y_rev_transform(x$inv_link(pr))
	}
}


#' Spatial MBM prediction
#' 
#' @param x A previously-fit MBM object
#' @param prdat New dataset to be used for prediction; either a raster stack or data 
#' 		frame. See 'Details'
#' @param coords matrix with 2 columns containing X-Y coordinates for \code{prdat}, 
#' 		required if prdat does not have a \code{coordinates} method.
#' @param method How to compute the spatial predictions; see 'Details'
#' @param ... Other named parameters to pass to \code{\link{predict.mbm}}.
#' @details \code{prdat} can either be a raster stack with new variables (and spatial
#' 		information) for prediction, or a data frame-like object with previous 
#' 		predictions from \code{\link{predict.mbm}} with 4 columns: 1. site1, 2. site2,
#' 		3. mean, and 4. sd. 
#' 
#'		For rasters, if a layer named 'names' is included (recommended), this layer will 
#' 		be used as sitenames, otherwise they will be assigned unique numbers.
#' 
#' 		If \code{method} is "slow", spatial predictions will be computed by first 
#' 		predicting dissimilarity to all pairs of raster cells, then performing an 
#' 		ordination on the dissimilarity matrix to produce an RGB raster of spatial 
#' 		predictions.
#'
#'		For method == 'fast' (currently not implemented), computation is sped up by first 
#'		performing hierarchical clustering on the predicted dissimilarity matrix for the 
#'		calibration data (which will have already been computed when mbm was run) to 
#'		produce cell categories. Each raster cell will then be assigned the category of 
#'		the calibration data point that is closest environmentally. Then, we compute the
#'		dissimilarity matrix of the categories (based on the mean environmental
#'		values). The ordination is performed as with the slow method on this
#' 		dissimilarity matrix.
#' @return An object of class mbmSP, which is a list with three named items: \code{fits}
#'		is a 3-band gridded SpatialPointsDataFrame giving the first three prinipal
#'		components of predicted pairwise dissimilarity, stdev is a SpatialPointsDataFrame
#'		giving the mean of pairwise dissimilarities among all other sites in a given site, 
#'		and pcoa is the principal coordinates analysis for the fits. Both fits and stdev 
#'		can be made into rasters using raster::stack() and raster::raster().
#' @export
spatial_predict <- function(x, prdat, coords, method = c('slow', 'fast'), ...)
{
	method <- match.arg(method)

	if(method == 'fast') {
		stop('Fast method is not implemented, use method="slow"')
	} else {
		if(inherits(prdat, "RasterStack") | inherits(prdat, "RasterLayer") | 
			inherits(prdat, "RasterBrick"))
		{
			preds <- predict_mbm_raster(x, prdat)
		} else {
			if(missing(coords))
				coords <- coordinates(prdat)
			fitMat <- make_symmetric(prdat, site1 ~ site2, value.var = "fit")
			sdMat <- make_symmetric(prdat, site1 ~ site2, value.var = "stdev")
			preds <- list(fits = fitMat, stdev = sdMat, coords = coords)
		}
		
		fits <- x$y_rev_transform(x$rev_link(preds$fits))
		# for fits, use a PCoA to collapse to 3 axes.
		# For stdev, use rowmeans to collapse to one
		fitPCoA <- ade4::dudi.pco(as.dist(fits), scannf = FALSE, nf = 3)
		fitPCoA_scaled <- as.data.frame(apply(as.matrix(fitPCoA$l1), 2, function(x) {
			x <- x - min(x)
			x / max(x)
		}))
		sdMeans <- data.frame(sd = rowMeans(preds$stdev, na.rm = TRUE))

		# make grids
		sp::coordinates(fitPCoA_scaled) <- sp::coordinates(sdMeans) <- preds$coords
		sp::gridded(fitPCoA_scaled) <- sp::gridded(sdMeans) <- TRUE
		ret <- list(fits = fitPCoA_scaled, stdev = sdMeans, pcoa = fitPCoA)
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


## deprecated - old hack, delete if nothing breaks
# #' Turn an MBM prediction dataframe into a symmetric matrix 
# #' @param DF MBM predfiction dataframe
# #' @param formula A formula to be passed to \code{link{reshape2::acast}}
# #' @param value.var Name of value variable for \code{acast}
# #' @param ... Additional parameters for \code{acast}
# #' @return A symmetric matrix of predictions
# #' @keywords internal
# make_symmetric <- function(DF, formula, value.var, ...)
# {
# 	mat <- reshape2::acast(DF, formula, value.var = value.var, ...)
# 	mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
# 	mat
# }


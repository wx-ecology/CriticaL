# ============================================================================
# HaversineLMEfunctions.R
# ----------------------------------------------------------------------------
# Custom spatial-correlation structure (`corHaversine`) for `nlme::lme()` and
# `nlme::gls()` that uses great-circle (Haversine) distances rather than the
# Euclidean distances assumed by the built-in nlme spatial-correlation
# classes (corSpher, corExp, corGaus, corLin, corRatio).
#
# This is required for global-scale analyses with longitude/latitude inputs:
# Euclidean distances on a longitude/latitude grid distort actual spatial
# proximity, especially at high latitudes and across the antimeridian.
#
# Usage:
#   library(nlme)
#   source("HaversineLMEfunctions.R")
#   m <- lme(y ~ x, random = ~ 1 | group,
#            correlation = corHaversine(form = ~ lon + lat, mimic = "corExp"),
#            data = mydata)
#
# Notes:
#   - Coordinate columns must be named exactly "lon" and "lat".
#   - Inputs are assumed in decimal degrees unless `radians = TRUE`.
#   - The `mimic` argument selects which built-in correlation function shape
#     to emulate ("corSpher", "corExp", "corGaus", "corLin", or "corRatio").
#
# Provenance:
#   Adapted from a community helper widely circulated for fitting nlme
#   models with great-circle distances. Test code that originally
#   accompanied this file (smoke tests using simulated data) has been
#   removed so that sourcing this file does not trigger unrelated model
#   fitting.
# ============================================================================


# ---- Distance functions ----------------------------------------------------
# Great-circle distance between two points (in radians) using the Haversine
# formula. Returns distance in km, assuming Earth radius = 6371 km.
haversine <- function(x0, x1, y0, y1) {
  a <- sin((y1 - y0) / 2)^2 + cos(y0) * cos(y1) * sin((x1 - x0) / 2)^2
  v <- 2 * asin(min(1, sqrt(a)))
  6371 * v
}

# Pairwise great-circle distance matrix for a two-column matrix of
# longitude/latitude. Input may be in decimal degrees (default) or radians.
# Note: `fields::rdist.earth()` is more efficient for very large inputs.
haversineDist <- function(xy, radians = FALSE) {
  if (ncol(xy) > 2) stop("Input must have two columns (longitude and latitude)")
  if (!radians) xy <- xy * pi / 180
  hMat <- matrix(NA, ncol = nrow(xy), nrow = nrow(xy))
  for (i in 1:nrow(xy)) {
    for (j in i:nrow(xy)) {
      hMat[j, i] <- haversine(xy[i, 1], xy[j, 1], xy[i, 2], xy[j, 2])
    }
  }
  as.dist(hMat)
}


# ---- Inherit standard methods from corSpatial ------------------------------
# Most of the machinery from nlme's corSpatial parent class works without
# modification; we simply alias these methods to corHaversine.
Initialize.corHaversine <- nlme:::Initialize.corSpatial
recalc.corHaversine     <- nlme:::recalc.corSpatial
Variogram.corHaversine  <- nlme:::Variogram.corSpatial
corFactor.corHaversine  <- nlme:::corFactor.corSpatial
corMatrix.corHaversine  <- nlme:::corMatrix.corSpatial
coef.corHaversine       <- nlme:::coef.corSpatial
"coef<-.corHaversine"   <- nlme:::"coef<-.corSpatial"


# ---- Constructor for the corHaversine class --------------------------------
corHaversine <- function(value  = numeric(0),
                         form   = ~ 1,
                         mimic  = "corExp",
                         nugget = FALSE,
                         fixed  = FALSE) {
  spClass <- "corHaversine"
  attr(value, "formula")  <- form
  attr(value, "nugget")   <- nugget
  attr(value, "fixed")    <- fixed
  attr(value, "function") <- mimic
  class(value) <- c(spClass, "corStruct")
  value
}
environment(corHaversine) <- asNamespace("nlme")


# ---- Dim method -------------------------------------------------------------
# Custom Dim method that uses the third element of the Dim list to store the
# spatial-class index, so that downstream methods know which built-in
# correlation shape to emulate.
Dim.corHaversine <- function(object, groups, ...) {
  if (missing(groups)) return(attr(object, "Dim"))
  val <- Dim.corStruct(object, groups)
  val[["start"]] <- c(0, cumsum(val[["len"]] * (val[["len"]] - 1) / 2)[-val[["M"]]])
  names(val)[3] <- "spClass"
  val[[3]] <- match(attr(object, "function"),
                    c("corSpher", "corExp", "corGaus", "corLin", "corRatio"),
                    0)
  val
}
environment(Dim.corHaversine) <- asNamespace("nlme")


# ---- getCovariate method ---------------------------------------------------
# Extracts the lon/lat covariate from the model frame and converts it to a
# pairwise Haversine-distance representation (grouped or ungrouped, depending
# on whether the formula includes a grouping factor).
getCovariate.corHaversine <- function(object, form = formula(object), data) {
  
  if (is.null(covar <- attr(object, "covariate"))) {
    
    if (missing(data)) stop("need data to calculate covariate")
    
    covForm <- getCovariateFormula(form)
    
    if (length(all.vars(covForm)) > 0) {
      
      # Drop intercept if present
      if (attr(terms(covForm), "intercept") == 1) {
        covForm <- eval(parse(text = paste0("~", deparse(covForm[[2]]), "-1")))
      }
      
      # Validate that the covariates are exactly "lon" and "lat"
      if (length(all.vars(covForm)) > 2) {
        stop("corHaversine can only take two covariates, 'lon' and 'lat'")
      }
      if (!all(all.vars(covForm) %in% c("lon", "lat"))) {
        stop("covariates must be named 'lon' and 'lat'")
      }
      
      covar <- as.data.frame(unclass(
        model.matrix(covForm,
                     model.frame(covForm, data, drop.unused.levels = TRUE))
      ))
      covar <- covar[, order(colnames(covar), decreasing = TRUE)]  # ensure lon, lat order
      
    } else {
      covar <- NULL
    }
    
    # Compute pairwise distances, optionally within groups
    if (!is.null(getGroupsFormula(form))) {
      
      grps <- getGroups(object, data = data)
      
      if (is.null(covar)) {
        covar <- lapply(split(grps, grps),
                        function(x) as.vector(dist(1:length(x))))
      } else {
        giveDist <- function(el) {
          el <- as.matrix(el)
          if (nrow(el) > 1) as.vector(haversineDist(el)) else numeric(0)
        }
        covar <- lapply(split(covar, grps), giveDist)
      }
      covar <- covar[sapply(covar, length) > 0]  # drop singleton groups
      
    } else {
      
      if (is.null(covar)) {
        covar <- as.vector(dist(1:nrow(data)))
      } else {
        covar <- as.vector(haversineDist(as.matrix(covar)))
      }
    }
    
    if (any(unlist(covar) == 0)) {
      stop("cannot have zero distances in \"corHaversine\"")
    }
  }
  
  covar
}
environment(getCovariate.corHaversine) <- asNamespace("nlme")

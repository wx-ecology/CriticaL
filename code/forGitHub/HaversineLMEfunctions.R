#### corHaversine - spatial correlation with haversine distance

# Calculates the geodesic distance between two points specified by radian latitude/longitude using Haversine formula.
# output in km
haversine <- function(x0, x1, y0, y1) {
  a <- sin( (y1 - y0)/2 )^2 + cos(y0) * cos(y1) * sin( (x1 - x0)/2 )^2
  v <- 2 * asin( min(1, sqrt(a) ) )
  6371 * v
}

# function to compute geodesic haversine distance given two-column matrix of longitude/latitude
# input is assumed in form decimal degrees if radians = F
# note fields::rdist.earth is more efficient
haversineDist <- function(xy, radians = F) {
  if (ncol(xy) > 2) stop("Input must have two columns (longitude and latitude)")
  if (radians == F) xy <- xy * pi/180
  hMat <- matrix(NA, ncol = nrow(xy), nrow = nrow(xy))
  for (i in 1:nrow(xy) ) {
    for (j in i:nrow(xy) ) {
      hMat[j,i] <- haversine(xy[i,1], xy[j,1], xy[i,2], xy[j,2]) 
    }
  }
  as.dist(hMat)
}

## for most methods, machinery from corSpatial will work without modification
Initialize.corHaversine <- nlme:::Initialize.corSpatial
recalc.corHaversine <- nlme:::recalc.corSpatial
Variogram.corHaversine <- nlme:::Variogram.corSpatial
corFactor.corHaversine <- nlme:::corFactor.corSpatial
corMatrix.corHaversine <- nlme:::corMatrix.corSpatial
coef.corHaversine <- nlme:::coef.corSpatial
"coef<-.corHaversine" <- nlme:::"coef<-.corSpatial"

## Constructor for the corHaversine class
corHaversine <- function(value = numeric(0), form = ~ 1, mimic = "corExp", nugget = FALSE, fixed = FALSE) {
  spClass <- "corHaversine"
  attr(value, "formula") <- form
  attr(value, "nugget") <- nugget
  attr(value, "fixed") <- fixed
  attr(value, "function") <- mimic
  class(value) <- c(spClass, "corStruct")
  value
}   # end corHaversine class
environment(corHaversine) <- asNamespace("nlme")

Dim.corHaversine <- function(object, groups, ...) {
  if (missing(groups)) return(attr(object, "Dim"))
  val <- Dim.corStruct(object, groups)
  val[["start"]] <- c(0, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
  ## will use third component of Dim list for spClass
  names(val)[3] <- "spClass"
  val[[3]] <- match(attr(object, "function"), c("corSpher", "corExp", "corGaus", "corLin", "corRatio"), 0)
  val
}
environment(Dim.corHaversine) <- asNamespace("nlme")


## getCovariate method for corHaversine class
getCovariate.corHaversine <- function(object, form = formula(object), data) {
  if (is.null(covar <- attr(object, "covariate"))) {          # if object lacks covariate attribute
    if (missing(data)) {                                    # if object lacks data
      stop("need data to calculate covariate")
    }
    covForm <- getCovariateFormula(form)
    if (length(all.vars(covForm)) > 0) {                    # if covariate present
      if (attr(terms(covForm), "intercept") == 1) {       # if formula includes intercept
        covForm <- eval(parse(text = paste("~", deparse(covForm[[2]]),"-1",sep="")))    # remove intercept
      }
      # can only take covariates with correct names
      if (length(all.vars(covForm)) > 2) stop("corHaversine can only take two covariates, 'lon' and 'lat'")
      if ( !all(all.vars(covForm) %in% c("lon", "lat")) ) stop("covariates must be named 'lon' and 'lat'")
      covar <- as.data.frame(unclass(model.matrix(covForm, model.frame(covForm, data, drop.unused.levels = TRUE) ) ) )
      covar <- covar[,order(colnames(covar), decreasing = T)] # order as lon ... lat
    }
    else {
      covar <- NULL
    }
    
    if (!is.null(getGroupsFormula(form))) {                 # if groups in formula extract covar by groups
      grps <- getGroups(object, data = data)
      if (is.null(covar)) {
        covar <- lapply(split(grps, grps), function(x) as.vector(dist(1:length(x) ) ) )     # filler?
      } 
      else {
        giveDist <- function(el) {
          el <- as.matrix(el)
          if (nrow(el) > 1) as.vector(haversineDist(el))                       
          else numeric(0)
        }
        covar <- lapply(split(covar, grps), giveDist )
      }
      covar <- covar[sapply(covar, length) > 0]  # no 1-obs groups
    } 
    else {                                  # if no groups in formula extract distance
      if (is.null(covar)) {
        covar <- as.vector(dist(1:nrow(data) ) )
      } 
      else {
        covar <- as.vector(haversineDist(as.matrix(covar) ) )
      }
    }
    if (any(unlist(covar) == 0)) {          # check that no distances are zero
      stop("cannot have zero distances in \"corHaversine\"")
    }
  }
  covar
}   # end method getCovariate
environment(getCovariate.corHaversine) <- asNamespace("nlme")
############################################################
## test that corHaversine runs with spherical correlation (not testing that it WORKS ...)
library(MASS)
library(nlme)
set.seed(1001)
sample_data <- data.frame(lon = -121:-22, lat = -50:49)
ran <- 1000 # 'range' parameter for spherical correlation
dist_matrix <- as.matrix(haversineDist(sample_data))    # haversine distance matrix
# set up correlation matrix of response
corr_matrix <- 1-1.5*(dist_matrix/ran)+0.5*(dist_matrix/ran)^3
corr_matrix[dist_matrix > ran] = 0
diag(corr_matrix) <- 1
# set up covariance matrix of response
sigma <- 2  # residual standard deviation
cov_matrix <- (diag(100)*sigma) %*% corr_matrix %*% (diag(100)*sigma)   # correlated response
# generate response
sample_data$y <- mvrnorm(1, mu = rep(0, 100), Sigma = cov_matrix)

# fit model
gls_haversine <- gls(y ~ 1, correlation = corHaversine(form=~lon+lat, mimic="corExp"), data = sample_data)
summary(gls_haversine)
############################################################
## test that corHaversine runs with Gaussian correlation
ran = 100  # parameter for Gaussian correlation
corr_matrix_gauss <- exp(-(dist_matrix/ran)^2)
diag(corr_matrix_gauss) <- 1
# set up covariance matrix of response
cov_matrix_gauss <- (diag(100)*sigma) %*% corr_matrix_gauss %*% (diag(100)*sigma)   # correlated response
# generate response
sample_data$y_gauss <- mvrnorm(1, mu = rep(0, 100), Sigma = cov_matrix_gauss)

# fit model
gls_haversine_gauss <- gls(y_gauss ~ 1, correlation = corHaversine(form=~lon+lat, mimic = "corExp"), data = sample_data)
summary(gls_haversine_gauss)

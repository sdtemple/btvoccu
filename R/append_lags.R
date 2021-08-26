#' Append Lagged Covariates
#'
#' Append lagged variables to a 4D array.
#'
#' @param x numeric array
#' @param covariate character name
#' @param lags integers
#'
#' @return numeric array
#'
#' @export
append_lags <- function(x, covariate, lags){
  covind <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  for(lag in lags){
    xx <- array(NA, dim = c(dim(x)[1:3], 1))
    for(i in 1:dim(x)[1]){
      for(j in 1:dim(x)[2]){
        for(k in (lag + 1):(dim(x)[3])){
          xx[i,j,k,1] <- x[i,j,(k-lag),covind]
        }
      }
    }
    x <- abind::abind(x, xx)
    dimnames(x)[[4]][length(dimnames(x)[[4]])] <- paste('lag', lag, covariate, sep ='')
  }
  return(x)
}

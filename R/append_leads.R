#' Append Leading Covariates
#'
#' Append leading variables to a 4D array.
#'
#' @param x numeric array
#' @param covariate character name
#' @param leads integers
#'
#' @return numeric array
#'
#' @export
append_leads <- function(x, covariate, leads){
  covind <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  for(lead in leads){
    xx <- array(NA, dim = c(dim(x)[1:3], 1))
    for(i in 1:dim(x)[1]){
      for(j in 1:dim(x)[2]){
        for(k in 1:(dim(x)[3] - lead)){
          xx[i,j,k,1] <- x[i,j,(k + lead),covind]
        }
      }
    }
    x <- abind::abind(x, xx)
    dimnames(x)[[4]][length(dimnames(x)[[4]])] <- paste('lead', lead, covariate, sep ='')
  }
  return(x)
}

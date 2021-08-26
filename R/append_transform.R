#' Append Transformed Variable
#'
#' Append a transformed variable to a 4D array.
#'
#' @param x numeric array
#' @param covariate character name
#' @param fn simple function
#' @param prefix character prefix for column name
#'
#' @return numeric array
#'
#' @export
append_transform <- function(x, covariate, fn, prefix){
  covind <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  x <- abind::abind(x, fn(x[,,,covind]))
  dimnames(x)[[4]][dim(x)[4]] <- paste(prefix, covariate, sep = '')
  return(x)
}

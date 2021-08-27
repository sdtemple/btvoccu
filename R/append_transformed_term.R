#' Append transformed variable
#'
#' Append a transformed variable to a 4D array.
#'
#' @param x numeric array
#' @param covariate character name
#' @param fn function
#' @param prefix character prefix for column name
#'
#' @return numeric array
#' 
#' @examples 
#' x <- append_transformed_term(x, "precipitation", sqrt, "sqrt.")
#' x <- append_transformed_term(x, "popdensity", log10, "log10.")
#'
#' @export
append_transformed_term <- function(x, covariate, fn, prefix){
  covind <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  x <- abind::abind(x, fn(x[,,,covind]))
  dimnames(x)[[4]][dim(x)[4]] <- paste(prefix, covariate, sep = '')
  return(x)
}

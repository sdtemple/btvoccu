#' Append interaction term
#'
#' Append an interaction variable to a 4D array.
#'
#' @param x numeric array
#' @param covariate character name
#' @param predictor character name
#'
#' @return numeric array
#' 
#' @examples 
#' x <- append_interaction_term(x, "temperature", "precipitation")
#'
#' @export
append_interaction_term <- function(x, covariate, predictor){
  covind <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  predind <- which(dimnames(x)[[4]] == predictor, arr.ind = T)
  x <- abind::abind(x, x[,,,covind] * x[,,,predind])
  dimnames(x)[[4]][dim(x)[4]] <- paste(covariate, '*', predictor, sep = '')
  return(x)
}

#' Subset 4D Array
#'
#' Select predictors for design arrays.
#'
#' @param x numeric array
#' @param d integer: array dimension
#' @param labels character vector: predictor names
#'
#' @return numeric array
#'
#' @export
subset_4darray <- function(x, d, labels){
  indices <- which(dimnames(x)[[d]] %in% labels, arr.ind = T)
  if(d == 1){
    return(x[indices,,,, drop = F])
  } else if(d == 2){
    return(x[,indices,,, drop = F])
  } else if(d == 3){
    return(x[,,indices,, drop = F])
  } else if(d == 4){
    return(x[,,,indices, drop = F])
  }
}

#' Subset 4D array
#'
#' Subset 4-dimensional arrays along any major dimension.
#'
#' @param x numeric array
#' @param d integer: array dimension
#' @param labels character vector: predictor names
#'
#' @return numeric array
#' 
#' @examples 
#' w <- subset_4darray(x, 1, c("TOR", "YRK", "PEE"))
#' w <- subset_4darray(x, 2, 2005:2010)
#' w <- subset_4darray(x, 3, 20:40)
#' w <- subset_4darray(x, 4, c("precipitation", "temperature", "popdensity"))
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
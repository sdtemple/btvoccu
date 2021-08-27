#' Q matrix for ICAR component
#'
#' Compute the Q matrix for an ICAR component.
#'
#' @param A adjacency matrix
#'
#' @return matrix
icarQ <- function(A){
  D <- diag(apply(A, 2, sum))
  return(D - A)
}
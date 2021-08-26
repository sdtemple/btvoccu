#' Append Spatial Random Effects
#'
#' Append eigenvectors from Moran operator decomposition
#' to a 4D array.
#'
#' @param x numeric array
#' @param arealeffs character vector: areal effects
#' @param A adjacency matrix
#'
#' @return numeric array
#'
#' @export
append_sp <- function(x, arealeffs, A){

  # local function
  spmatrix <- function(X, A){
    P <- X %*% tcrossprod(solve(crossprod(X, X)), X)
    I <- diag(nrow(P))
    M <- (I - P) %*% A %*% (I - P)
    return(eigen(M))
  }

  xx <- subset_4darray(x, 4, arealeffs)
  eM <- spmatrix(xx[,1,1,], A)
  print(eM$values)

  sparray <- array(NA,
                   dim = c(dim(x)[1:3],
                               ncol(eM$vectors)),
                   dimnames = list(dimnames(x)[[1]],
                                   dimnames(x)[[2]],
                                   dimnames(x)[[3]],
                                   paste('sp',
                                         1:ncol(eM$vectors),
                                         sep = '')
                                   )
                   )
  for(j in 1:dim(x)[2]){
    for(k in 1:dim(x)[3]){
      sparray[,j,k,] <- eM$vectors
    }
  }
  return(abind::abind(x, sparray))
}

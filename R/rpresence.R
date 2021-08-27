#' Simulate presence-absence data
#'
#' Simulate presence-absence data
#' based on presence probabilities.
#'
#' @param pre numeric array: presence probabilities
#'
#' @return binary array
#' 
#' @examples 
#' sim <- rpresence(posterior_simulate(mdl, x, c("TOR", "YRK", "PEE"), 2004:2005, 20:40))
#' 
#' @export
rpresence <- function(pre){
  yrep <- array(NA, dim = c(dim(pre)[1:3], 1),
                dimnames = list(dimnames(pre)[[1]],
                                dimnames(pre)[[2]],
                                dimnames(pre)[[3]],
                                "Response")
  )
  d <- sample(dim(pre)[4], 1)
  for(i in 1:dim(pre)[1]){
    for(j in 1:dim(pre)[2]){
      for(k in 1:dim(pre)[3]){
        yrep[i,j,k,1] <- rbinom(1, 1, pre[i,j,k,d])
      }
    }
  }
  
  return(yrep)
}
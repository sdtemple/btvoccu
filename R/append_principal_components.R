#' Append principal components
#'
#' Append principal components to a 4D array.
#' Plots from \code{factoextra} describe the PCA.
#'
#' @param x numeric array
#' @param covariates vector of covariates
#' @param npc integer: number of principal components
#' @param prefix character prefix for column name
#'
#' @return numeric array
#' 
#' @examples 
#' x <- append_principal_components(x, dimnames(x)[[4]][12:23], 2, "climate")
#'
#' @export
append_principal_components <- function(x, covariates, npc, prefix){
  
  z <- apply(x, 4, as.numeric)
  indices <- which(dimnames(z)[[2]] %in% covariates, arr.ind = T)
  components <- prcomp(~., data=data.frame(z[,indices]), na.action = na.exclude)
  print(factoextra::fviz_eig(components))
  print(factoextra::fviz_pca_var(components, repel = T))
  z.pca <- factoextra::get_pca_ind(components)
  x.pca <- array(dim = c(dim(x)[1:3], npc),
                 dimnames = list(dimnames(x)[[1]], dimnames(x)[[2]], dimnames(x)[[3]],
                                 paste(prefix, '.pc', 1:npc, sep = '')))
  
  for(l in 1:npc){
    n <- 1
    for(k in 1:dim(x)[3]){
      for(j in 1:dim(x)[2]){
        for(i in 1:dim(x)[1]){
          x.pca[i,j,k,l] <- z.pca$coord[n,l]
          n <- n + 1
        }
      }
    }
  }
  
  return(abind::abind(x, x.pca))
}
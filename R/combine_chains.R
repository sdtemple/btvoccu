#' Combine MCMC chains
#'
#' Combine chains for mixing diagnostics
#' and posterior inference.
#'
#' @param model \code{btvoccu} model object
#' @param effect "betas", "alphas", or "thetas"
#' @param mnridx integer index
#' @param burnin percent of posterior samples to burn
#' @param thin keep every nth posterior draw
#'
#' @return vector of posterior samples
#' 
#' @examples 
#' vec <- combine_chains(mdl, "betas", 1, burnin = 0.10)
#' 
#' @export
combine_chains <- function(model, effect, mnridx, burnin = .5, thin = 1){
  mjridx <- which(names(model) == effect, arr.ind = T)
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  vec <- as.vector(model[[mjridx]][1, keep, mnridx])
  for(n in 2:model$nchains){vec <- c(vec, as.vector(model[[mjridx]][n, keep, mnridx]))}
  return(vec)
}

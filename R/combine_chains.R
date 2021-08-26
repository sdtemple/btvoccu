#' Combine MCMC Chains
#'
#' Combine chains for mixing diagnostics
#' and posterior inference.
#'
#' @param model \code{btvoccu} output
#' @param effect "betas", "alphas", or "thetas"
#' @param mnridx integer index
#' @param burnin percent of posterior samples to burn
#' @param thin integer: default \code{1} implies no thinning, \code{2} means keep every other draw, etc.
#'
#' @return vector of posterior samples
combine_chains <- function(model,
                           effect,
                           mnridx,
                           burnin = .5,
                           thin = 1
                           ){

  mjridx <- which(names(model) == effect, arr.ind = T)

  keep <- seq(ceiling(burnin * model$niter),
              model$niter,
              by = thin)

  vec <- as.vector(model[[mjridx]][1, keep, mnridx])
  for(n in 2:model$nchains){
    vec <- c(vec,
             as.vector(model[[mjridx]][n, keep, mnridx]))
  }

  return(vec)
}

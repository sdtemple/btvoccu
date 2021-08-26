#' Posterior Correlations
#'
#' Summarize posterior correlations.
#'
#' @param model \code{btvoccu} output
#' @param Sigmas "Sigmabetas" or "Sigmaalphas"
#' @param burnin percent of posterior samples to burn
#' @param thin integer: default \code{1} implies no thinning, \code{2} means keep every other draw, etc.
#' @param credible vector of quantiles
#'
#' @return matrix
#'
#' @export
posterior_corrs <- function(model,
                            Sigmas,
                            burnin = .5,
                            thin = 1,
                            credible = c(.025,.5,.975)
                            ){

  mjridx <- which(names(model) == Sigmas, arr.ind = T)
  out <- model[[mjridx]]
  keep <- seq(ceiling(burnin * model$niter),
              model$niter, by = thin)

  n <- 0
  table <- matrix(nrow = 0, ncol = (length(credible) + 2))
  for(i in 1:(dim(out)[3] - 1)){
    for(j in (i+1):dim(out)[3]){

      # format as coda mcmc object
      mcmc <- coda::mcmc.list()
      chains <- apply(out[,keep,i,j] /
                        sqrt(out[,keep,i,i] *
                               out[,keep,j,j]),
                      1,
                      as.numeric)

      for(k in 1:dim(chains)[2]){
        mcmc[[k]] <- coda::as.mcmc(chains[,k])
      }

      # check chain convergence
      rhat <- coda::gelman.diag(mcmc, autoburnin = F)[[1]][1]

      # combine chains
      vec <- as.vector(mcmc[[1]])
      for(k in 2:dim(chains)[2]){
        vec <- c(vec,
                 as.vector(mcmc[[k]]))
      }

      # store in table
      n <- n + 1
      table <- rbind(table,
                     c(quantile(vec, credible),
                       mean(vec),
                       rhat
                      )
                    )
      rownames(table)[n] <- paste(substr(names(model)[mjridx],
                                         1,
                                         (nchar(names(model)[mjridx])-1)),
                                  i,
                                  j,
                                  sep = '')
    }
  }
  colnames(table) <- c(credible, 'mean', 'rhat')

  return(table)
}

#' Posterior Effects
#'
#' Summarize posterior effects.
#'
#' @param model \code{btvoccu} output
#' @param burnin percent of posterior samples to burn
#' @param thin integer: default \code{1} implies no thinning, \code{2} means keep every other draw, etc.
#' @param credible quantiles
#'
#' @return matrix
#'
#' @export
posterior_effects <- function(model,
                              credible = c(.025, .50, .975),
                              burnin = .5,
                              thin = 1
                              ){

  if(is.null(model$speffs)){
    effects <- c('betas','alphas')
  } else{
    effects <- c('betas','alphas','thetas')
  }

  keep <- seq(ceiling(burnin * model$niter),
              model$niter,
              by = thin)
  nchains <- model$nchains

  n <- 0
  table <- matrix(nrow = 0, ncol = (length(credible) + 2))
  indices <- which(names(model) %in% effects, arr.ind = T)
  for(i in indices){
    for(j in 1:dim(model[[i]])[3]){

      # combine chains
      combined <- combine_chains(model,
                                 names(model)[i],
                                 j,
                                 burnin,
                                 thin)

      # convert to mcmc object for coda::gelman.diag
      mcmc <- coda::mcmc.list()
      for(k in 1:nchains){
        mcmc[[k]] <- coda::as.mcmc(model[[i]][k, keep, j])
      }

      # store in table
      table <- rbind(table,
                     c(quantile(combined, credible),
                       mean(combined),
                       coda::gelman.diag(mcmc, autoburnin = F)[[1]][1]
                      )
                    )
      n <- n + 1
      rownames(table)[n] <- paste(substr(names(model)[i],
                                         1,
                                         (nchar(names(model)[i])-1)),
                                  j,
                                  sep = '')
    }
  }
  colnames(table) <- c(credible, 'mean', 'rhat')

  return(table)
}

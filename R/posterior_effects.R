#' Posterior effects
#'
#' Summarize posterior beta, alpha, and theta coefficients.
#'
#' @param model \code{btvoccu} model object
#' @param burnin percent of posterior samples to burn
#' @param thin keep every nth posterior draw
#' @param credible quantiles
#'
#' @return matrix
#' 
#' @examples
#' View(posterior_effects(mdl, 0.10, 2, c(0.25, 0.5, 0.75)))
#' xtable(posterior_effects(mdl))
#'
#' @export
posterior_effects <- function(model, burnin = .5, thin = 1, credible = c(.025, .50, .975)){
  
  if(is.null(model$speffs)){
    effects <- c("betas","alphas")
  } else{
    effects <- c("betas","alphas","thetas")
  }
  
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  nchains <- model$nchains
  
  n <- 0
  table <- matrix(nrow = 0, ncol = (length(credible) + 3))
  indices <- which(names(model) %in% effects, arr.ind = T)
  for(i in indices){
    nameidx <- names(model)[i+6]
    for(j in 1:dim(model[[i]])[3]){
      
      # combine chains
      combined <- combine_chains(model, names(model)[i], j, burnin, thin)
      
      if(nchains > 1){
        # convert to mcmc object for coda::gelman.diag
        mcmc <- coda::mcmc.list()
        for(k in 1:nchains){mcmc[[k]] <- coda::as.mcmc(model[[i]][k, keep, j])}
        
        # store in table
        table <- rbind(table, c(getElement(model, nameidx)[j],
                                quantile(combined, credible),mean(combined),
                                coda::gelman.diag(mcmc, autoburnin = F)[[1]][1]))
      } else{
        # store in table
        table <- rbind(table, c(getElement(model, nameidx)[j],
                                quantile(combined, credible),mean(combined),
                                NA))
      }
      
      n <- n + 1
      rownames(table)[n] <- paste(substr(names(model)[i],
                                         1, (nchar(names(model)[i])-1)),
                                  j, sep = "")
      
    }
  }
  
  colnames(table) <- c("covariate", credible, "mean", "rhat")
  
  return(table)
}

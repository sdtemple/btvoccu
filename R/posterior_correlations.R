#' Posterior correlations
#'
#' Summarize posterior correlations for beta and alpha coefficients.
#'
#' @param model \code{btvoccu} model object
#' @param Sigmas "Sigmabetas" or "Sigmaalphas"
#' @param burnin percent of posterior samples to burn
#' @param thin keep every nth posterior draw
#' @param credible quantiles
#'
#' @return matrix
#' 
#' @examples 
#' View(posterior_correlations(mdl, "Sigmabetas"))
#' xtable(posterior_correlations(mdl, "Sigmaalphas", burnin = 0.10))
#'
#' @export
posterior_correlations <- function(model, Sigmas, burnin = .5, thin = 1, credible = c(.025,.5,.975)){
  
  mjridx <- which(names(model) == Sigmas, arr.ind = T)
  nameidx <- names(model)[mjridx+3]
  out <- model[[mjridx]]
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  
  n <- 0
  table <- matrix(nrow = 0, ncol = (length(credible) + 4))
  for(i in 1:(dim(out)[3] - 1)){
    namei <- getElement(model, nameidx)[i]
    for(j in (i+1):dim(out)[3]){
      namej <- getElement(model, nameidx)[j]
      
      # format as coda mcmc object
      mcmc <- coda::mcmc.list()
      chains <- apply(out[,keep,i,j, drop=F] / sqrt(out[,keep,i,i, drop=F] * out[,keep,j,j, drop=F]),
                      1, as.numeric)
      nchains <- dim(chains)[2]
      for(k in 1:nchains){mcmc[[k]] <- coda::as.mcmc(chains[,k])}
      
      # check chain convergence
      rhat <-  NA
      if(nchains > 1){rhat <- coda::gelman.diag(mcmc, autoburnin = F)[[1]][1]}
      
      # combine chains
      vec <- as.vector(mcmc[[1]])
      if(nchains > 1){for(k in 2:nchains){vec <- c(vec,as.vector(mcmc[[k]]))}}
      
      
      # store in table
      n <- n + 1
      table <- rbind(table, c(namei, namej, quantile(vec, credible), mean(vec), rhat))
      rownames(table)[n] <- paste(substr(names(model)[mjridx],
                                         1, (nchar(names(model)[mjridx])-1)),
                                  i, j, sep = "")
    }
  }
  colnames(table) <- c("covariate1", "covariate2", credible, "mean", "rhat")
  
  return(table)
}

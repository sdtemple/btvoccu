#' Posterior variances
#'
#' Summarize posterior variances.
#'
#' @param model \code{btvoccu} model object
#' @param Sigmas "Sigmabetas", "Sigmaalphas", or "sigmathetas"
#' @param burnin percent of posterior samples to burn
#' @param thin keep every nth posterior draw
#' @param credible quantiles
#'
#' @return matrix
#' 
#' @examples 
#' posterior_variances(mdl, "sigmathetas")
#' View(posterior_variances(mdl, "Sigmabetas"))
#' xtable(posterior_variances(mdl, "Sigmaalphas", burnin = 0.10))
#'
#' @export
posterior_variances <- function(model, Sigmas, burnin = .5, thin = 1, credible = c(.025,.5,.975)){
  
  mjridx <- which(names(model) == Sigmas, arr.ind = T)
  nameidx <- names(model)[mjridx+3]
  out <- model[[mjridx]]
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  nchains <- model$nchains
  
  if(Sigmas == "sigmathetas"){
    if(is.null(model$sigmathetas)){
      stop("This model is not a spatial model.")
    }
    
    # format as coda mcmc object
    mcmc <- coda::mcmc.list()
    for(k in 1:nchains){
      mcmc[[k]] <- out[k, keep, 1]
    }
    
    # check chain convergence
    rhat <- coda::gelman.diag(mcmc, autoburnin = F)[[1]][1]
    
    # combine chains
    vec <- as.vector(out[1, keep, 1])
    for(k in 2:nchains){vec <- c(vec, as.vector(out[k, keep, 1]))}
    
    # store in table
    table <- matrix(nrow = 0, ncol = (length(credible) + 2))
    table <- rbind(table, c(quantile(vec, credible), mean(vec), rhat))
    rownames(table)[1] <- "sigmathetasq"
    colnames(table) <- c(credible, "mean", "rhat")
    
    return(table)
  }
  
  n <- 0
  table <- matrix(nrow = 0, ncol = (length(credible) + 3))
  for(i in 1:(dim(out)[3])){
    name <- getElement(model, nameidx)[i]
    
    # format as coda mcmc object
    mcmc <- coda::mcmc.list()
    chains <- apply(out[,keep,i,i, drop=F], 1, as.numeric)
    for(k in 1:nchains){mcmc[[k]] <- coda::as.mcmc(chains[,k])}
    
    # check chain convergence
    rhat <- NA
    if(nchains > 1){rhat <- coda::gelman.diag(mcmc, autoburnin = F)[[1]][1]}
    
    # combine chains
    vec <- as.vector(mcmc[[1]])
    if(nchains > 1){for(k in 2:nchains){vec <- c(vec, as.vector(mcmc[[k]]))}}
    
    # store in table
    n <- n + 1
    table <- rbind(table, c(name, quantile(vec, credible), mean(vec), rhat))
    rownames(table)[n] <- paste(substr(names(model)[mjridx],
                                       1, (nchar(names(model)[mjridx])-1)),
                                i, i, sep = "")
  }
  colnames(table) <- c("covariate", credible, "mean", "rhat")
  
  return(table)
}

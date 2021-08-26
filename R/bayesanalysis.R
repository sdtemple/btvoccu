# Analyze fitted Bayesian model
# August 25, 2021

# Drawing posteriors ----------------------------------------------------------

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
combine_chains <- function(model, effect, mnridx, burnin = .5, thin = 1){
  mjridx <- which(names(model) == effect, arr.ind = T)
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  vec <- as.vector(model[[mjridx]][1, keep, mnridx])
  for(n in 2:model$nchains){vec <- c(vec, as.vector(model[[mjridx]][n, keep, mnridx]))}
  return(vec)
}

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
#' posterior_effects(mdl, 0.10, 2, c(0.25, 0.5, 0.75))
#'
#' @export
posterior_effects <- function(model, burnin = .5, thin = 1, credible = c(.025, .50, .975)){
  
  if(is.null(model$speffs)){
    effects <- c('betas','alphas')
  } else{
    effects <- c('betas','alphas','thetas')
  }
  
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  nchains <- model$nchains
  
  n <- 0
  table <- matrix(nrow = 0, ncol = (length(credible) + 2))
  indices <- which(names(model) %in% effects, arr.ind = T)
  for(i in indices){
    for(j in 1:dim(model[[i]])[3]){
      
      # combine chains
      combined <- combine_chains(model, names(model)[i], j, burnin, thin)
      
      # convert to mcmc object for coda::gelman.diag
      mcmc <- coda::mcmc.list()
      for(k in 1:nchains){mcmc[[k]] <- coda::as.mcmc(model[[i]][k, keep, j])}
      
      # store in table
      table <- rbind(table, c(quantile(combined, credible),mean(combined),
                              coda::gelman.diag(mcmc, autoburnin = F)[[1]][1]))
      n <- n + 1
      rownames(table)[n] <- paste(substr(names(model)[i],
                                         1, (nchar(names(model)[i])-1)),
                                  j, sep = '')
    }
  }
  colnames(table) <- c(credible, 'mean', 'rhat')
  
  return(table)
}

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
#' @export
posterior_correlations <- function(model, Sigmas, burnin = .5, thin = 1, credible = c(.025,.5,.975)){
  
  mjridx <- which(names(model) == Sigmas, arr.ind = T)
  out <- model[[mjridx]]
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  
  n <- 0
  table <- matrix(nrow = 0, ncol = (length(credible) + 2))
  for(i in 1:(dim(out)[3] - 1)){
    for(j in (i+1):dim(out)[3]){
      
      # format as coda mcmc object
      mcmc <- coda::mcmc.list()
      chains <- apply(out[,keep,i,j] / sqrt(out[,keep,i,i] * out[,keep,j,j]),
                      1, as.numeric)
      
      for(k in 1:dim(chains)[2]){mcmc[[k]] <- coda::as.mcmc(chains[,k])}
      
      # check chain convergence
      rhat <- coda::gelman.diag(mcmc, autoburnin = F)[[1]][1]
      
      # combine chains
      vec <- as.vector(mcmc[[1]])
      for(k in 2:dim(chains)[2]){vec <- c(vec,as.vector(mcmc[[k]]))}
      
      # store in table
      n <- n + 1
      table <- rbind(table, c(quantile(vec, credible), mean(vec), rhat))
      rownames(table)[n] <- paste(substr(names(model)[mjridx],
                                         1, (nchar(names(model)[mjridx])-1)),
                                  i, j, sep = '')
    }
  }
  colnames(table) <- c(credible, 'mean', 'rhat')
  
  return(table)
}


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
#' @export
posterior_variances <- function(model, Sigmas, burnin = .5, thin = 1, credible = c(.025,.5,.975)){
  
  mjridx <- which(names(model) == Sigmas, arr.ind = T)
  out <- model[[mjridx]]
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  nchains <- model$nchains
  
  if(Sigmas == 'sigmathetas'){
    if(is.null(model$sigmathetas)){
      stop('This model is not a spatial model.')
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
    rownames(table)[1] <- 'sigmathetasq'
    colnames(table) <- c(credible, 'mean', 'rhat')
    
    return(table)
  }
  
  n <- 0
  table <- matrix(nrow = 0, ncol = (length(credible) + 2))
  for(i in 1:(dim(out)[3])){
    
    # format as coda mcmc object
    mcmc <- coda::mcmc.list()
    chains <- apply(out[,keep,i,i], 1, as.numeric)
    for(k in 1:nchains){mcmc[[k]] <- coda::as.mcmc(chains[,k])}
    
    # check chain convergence
    rhat <- coda::gelman.diag(mcmc, autoburnin = F)[[1]][1]
    
    # combine chains
    vec <- as.vector(mcmc[[1]])
    for(k in 2:nchains){vec <- c(vec, as.vector(mcmc[[k]]))}
    
    # store in table
    n <- n + 1
    table <- rbind(table, c(quantile(vec, credible), mean(vec), rhat))
    rownames(table)[n] <- paste(substr(names(model)[mjridx],
                                       1, (nchar(names(model)[mjridx])-1)),
                                i, i, sep = '')
  }
  colnames(table) <- c(credible, 'mean', 'rhat')
  
  return(table)
}

#' Trace plots
#'
#' Plot trace of posterior samples.
#'
#' @param model \code{btvoccu} model object
#' @param effects "betas", "alphas", or "thetas"
#' @param mnridx integer index
#' @param burnin percent of posterior samples to burn
#' 
#' @examples 
#' trace_plot(mdl, "betas", 1)
#' trace_plot(mdl, "alphas", 1, 0)
#'
#' @export
plot_trace <- function(model, effects, mnridx, burnin = .5){
  
  mjridx <- which(names(model) == effects, arr.ind = T)
  colors <- rainbow(model$nchains)
  burn <- -(1:(floor(burnin * model$niter)))
  
  # ylim
  mx <- max(model[[mjridx]][1, burn, mnridx])
  mn <- min(model[[mjridx]][1, burn, mnridx])
  for(n in 2:model$nchains){
    mx <- max(mx, model[[mjridx]][n, burn, mnridx])
    mn <- min(mn, model[[mjridx]][n, burn, mnridx])
  }
  
  # plot
  plot(model[[mjridx]][1, burn, mnridx],
       ylab = paste(substr(effects, 1, (nchar(effects)-1)),
                    mnridx,
                    sep = ''),
       xlab = 'index',
       type = 'l',
       col = colors[1],
       ylim = c(mn, mx))
  for(n in 2:model$nchains){lines(model[[mjridx]][n, burn, mnridx], col = colors[n])}
}

# Evaluating models --------------------------------------------------------

#' Watanabe Akaike information criterion
#'
#' Compute Watanabe-Akaike information criterion.
#' When indices are set to \code{NULL},
#' WAIC is calculated using the training data.
#'
#' @param model \code{btvoccu} model object
#' @param x numeric array: covariate array
#' @param y binary array: response array
#' @param sites vector
#' @param seasons vector
#' @param periods vector
#' @param burnin percent of posterior samples to burn
#' @param thin keep every nth posterior draw
#'
#' @return list
#'
#' @export
waic_score <- function(model, x, y,
                       sites = NULL,
                       seasons = NULL,
                       periods = NULL,
                       burnin = .5,
                       thin = 1){
  
  # setup data
  if(is.null(sites)){
    sites <- model$sites
  }
  if(is.null(seasons)){
    seasons <- model$seasons
  }
  if(is.null(periods)){
    periods <- model$periods
  }
  
  sites <- sort(sites)
  seasons <- sort(seasons)
  periods <- sort(periods)
  
  y <- subset_4darray(y, 1, sites)
  y <- subset_4darray(y, 2, seasons)
  y <- subset_4darray(y, 3, periods)
  
  x <- subset_4darray(x, 1, sites)
  x <- subset_4darray(x, 2, seasons)
  x <- subset_4darray(x, 3, periods)
  
  occueffs <- model$occueffs
  deteffs <- model$deteffs
  speffs <- model$speffs
  
  X <- subset_4darray(x, 4, occueffs)
  W <- subset_4darray(x, 4, deteffs)
  if(is.null(speffs)){
    M <- NULL
  } else{
    M <- subset_4darray(x, 4, speffs)
  }
  
  # access posterior samples
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  nsample <- model$nchains * length(keep)
  if(is.null(speffs)){
    effects <- c('betas','alphas') # no spatial effects
    indices <- which(effects %in% names(model), arr.ind = T)
    samples <- list()
    for(m in 1:length(indices)){
      neffects <- dim(model[[m]])[3]
      combined <- array(NA, dim = c(nsample, neffects))
      for(n in 1:neffects){
        combined[,n] <- combine_chains(model, effects[indices[m]], n, burnin, thin)
      }
      samples[[m]] <- combined
      names(samples)[m] <- c('betas','alphas')[indices[m]]
    }
  } else{
    effects <- c('betas','alphas','thetas') # spatial effects included
    indices <- which(effects %in% names(model), arr.ind = T)
    samples <- list()
    for(m in 1:length(indices)){
      neffects <- dim(model[[m]])[3]
      combined <- array(NA, dim = c(nsample, neffects))
      for(n in 1:neffects){
        combined[,n] <- combine_chains(model, effects[indices[m]], n, burnin, thin)
      }
      samples[[m]] <- combined
      names(samples)[m] <- c('betas','alphas','thetas')[indices[m]]
    }
  }
  
  # compute and store posterior glm output
  occu <- array(NA, dim = c(dim(x)[1:3], nsample),
                dimnames = list(sites, seasons, periods, 1:nsample))
  det <- array(NA, dim = c(dim(x)[1:3], nsample),
               dimnames = list(sites, seasons, periods, 1:nsample))
  yll <- array(NA, dim = c(dim(y)[1:3], nsample))
  if(is.null(M)){ # without spatial random effects
    for(i in 1:length(sites)){
      for(j in 1:length(seasons)){
        for(k in length(periods):1){
          if(is.na(y[i,j,k,1])){
            yll[i,j,k,] <- NA
          } else{
            for(d in 1:nsample){
              occueffs <- samples[[1]][d,]
              deteffs <- samples[[2]][d,]
              occu[i,j,k,d] <- model$link(X[i,j,k,] %*% occueffs)
              det[i,j,k,d] <- model$link(W[i,j,k,] %*% deteffs)
              if(y[i,j,k,1]){
                yll[i,j,k,d] <- occu[i,j,k,d] * det[i,j,k,d]
              } else{
                yll[i,j,k,d] <- occu[i,j,k,d] * (1 - det[i,j,k,d]) + (1 - occu[i,j,k,d])
              }
            }
          }
        }
      }
    }
  } else{ # with spatial random effects
    for(i in 1:length(sites)){
      for(j in 1:length(seasons)){
        for(k in length(periods):1){
          if(is.na(y[i,j,k,1])){
            yll[i,j,k,] <- NA
          } else{
            for(d in 1:nsample){
              occueffs <- samples[[1]][d,]
              deteffs <- samples[[2]][d,]
              speffs <- samples[[3]][d,]
              occu[i,j,k,d] <- model$link(X[i,j,k,] %*% occueffs + M[i,j,k,] %*% speffs)
              det[i,j,k,d] <- model$link(W[i,j,k,] %*% deteffs)
              if(y[i,j,k,1]){
                yll[i,j,k,d] <- occu[i,j,k,d] * det[i,j,k,d]
              } else{
                yll[i,j,k,d] <- occu[i,j,k,d] * (1 - det[i,j,k,d]) + (1 - occu[i,j,k,d])
              }
            }
          }
        }
      }
    }
  }
  
  m <- 0
  llarray <- array(NA, dim = c(nsample, sum(!is.na(y))))
  for(i in 1:length(sites)){
    for(j in 1:length(seasons)){
      for(k in 1:length(periods)){
        if(!is.na(y[i,j,k,1])){
          m <- m + 1
          llarray[,m] <- yll[i,j,k,]
        }
      }
    }
  }
  a <- log(apply(llarray, 2, mean))
  b <- apply(apply(llarray, c(1,2), log), 2, mean)
  pwaic <- 2 * sum(a - b)
  llpd <- sum(a)
  ellpd <- llpd - pwaic
  
  output <- list()
  output[[1]] <- -2 * ellpd
  output[[2]] <- ellpd
  output[[3]] <- llpd
  output[[4]] <- pwaic
  names(output) <- c('waic','ellpd','llpd','pwaic')
  
  return(output)
}

#' Posterior sigmoids
#'
#' Compute posterior occupancy, detection,
#' and presence probabilities.
#'
#' @param model \code{btvoccu} model object
#' @param x numeric covariate array
#' @param sites vector
#' @param seasons vector
#' @param periods vector
#' @param ndraws number of posterior draws
#' @param burnin percent of posterior samples to burn
#' @param thin keep every nth draw
#'
#' @return numeric array
#'
#' @export
posterior_sigmoids <- function(model, x,
                               sites,
                               seasons,
                               periods,
                               ndraws = 1000,
                               burnin = .5,
                               thin = 1){
  
  # setup data
  sites <- sort(sites)
  seasons <- sort(seasons)
  periods <- sort(periods)
  
  y <- subset_4darray(y, 2, seasons)
  y <- subset_4darray(y, 3, periods)
  
  x <- subset_4darray(x, 1, sites)
  x <- subset_4darray(x, 2, seasons)
  x <- subset_4darray(x, 3, periods)
  
  occueffs <- model$occueffs
  deteffs <- model$deteffs
  speffs <- model$speffs
  
  X <- subset_4darray(x, 4, occueffs)
  W <- subset_4darray(x, 4, deteffs)
  if(is.null(speffs)){
    M <- NULL
  } else{
    M <- subset_4darray(x, 4, speffs)
  }
  
  # access posterior samples
  keep <- seq(ceiling(burnin * model$niter), model$niter, by = thin)
  nsample <- model$nchains * length(keep)
  if(is.null(speffs)){
    effects <- c('betas','alphas') # no spatial effects
    indices <- which(effects %in% names(model), arr.ind = T)
    samples <- list()
    for(m in 1:length(indices)){
      neffects <- dim(model[[m]])[3]
      combined <- array(NA, dim = c(nsample, neffects))
      for(n in 1:neffects){
        combined[,n] <- combine_chains(model, effects[indices[m]], n, burnin, thin)
      }
      samples[[m]] <- combined
      names(samples)[m] <- c('betas','alphas')[indices[m]]
    }
  } else{
    effects <- c('betas','alphas','thetas') # spatial effects included
    indices <- which(effects %in% names(model), arr.ind = T)
    samples <- list()
    for(m in 1:length(indices)){
      neffects <- dim(model[[m]])[3]
      combined <- array(NA, dim = c(nsample, neffects))
      for(n in 1:neffects){
        combined[,n] <- combine_chains(model, effects[indices[m]], n, burnin, thin)
      }
      samples[[m]] <- combined
      names(samples)[m] <- c('betas','alphas','thetas')[indices[m]]
    }
  }
  
  # compute and store posterior glm output
  occu <- array(NA, dim = c(dim(x)[1:3], ndraws),
                dimnames = list(sites, seasons, periods, 1:ndraws))
  det <- array(NA, dim = c(dim(x)[1:3], ndraws),
               dimnames = list(sites, seasons, periods, 1:ndraws))
  pre <- array(NA, dim = c(dim(x)[1:3], ndraws),
               dimnames = list(sites, seasons, periods, 1:ndraws))
  if(is.null(M)){ # without spatial random effects
    for(i in 1:length(sites)){
      for(j in 1:length(seasons)){
        for(k in length(periods):1){
          for(d in 1:ndraws){
            sampleid <- sample(nsample, 1)
            occueffs <- samples[[1]][sampleid,]
            deteffs <- samples[[2]][sampleid,]
            occu[i,j,k,d] <- model$link(X[i,j,k,] %*% occueffs)
            det[i,j,k,d] <- model$link(W[i,j,k,] %*% deteffs)
            pre[i,j,k,d] <- occu[i,j,k,d] * det[i,j,k,d]
          }
        }
      }
    }
  } else{ # with spatial random effects
    for(i in 1:length(sites)){
      for(j in 1:length(seasons)){
        for(k in length(periods):1){
          for(d in 1:ndraws){
            sampleid <- sample(nsample, 1)
            occueffs <- samples[[1]][sampleid,]
            deteffs <- samples[[2]][sampleid,]
            speffs <- samples[[3]][sampleid,]
            occu[i,j,k,d] <- model$link(X[i,j,k,] %*% occueffs +
                                          M[i,j,k,] %*% speffs)
            det[i,j,k,d] <- model$link(W[i,j,k,] %*% deteffs)
            pre[i,j,k,d] <- occu[i,j,k,d] * det[i,j,k,d]
          }
        }
      }
    }
  }
  
  # hold detection probability when detection covariate is missing
  for(i in 1:length(sites)){
    for(j in 1:length(seasons)){
      for(k in length(periods):2){
        for(d in 1:ndraws){
          if(is.na(det[i,j,k,d])){
            if(!is.na(det[i,j,(k-1),d])){
              det[i,j,k,d] <- det[i,j,(k-1),d]
              pre[i,j,k,d] <- occu[i,j,k,d] * det[i,j,k,d]
            }
          }
        }
      }
    }
  }
  
  sigmoids <- list()
  sigmoids[[1]] <- occu
  sigmoids[[2]] <- det
  sigmoids[[3]] <- pre
  names(sigmoids) <- c('Occupancy','Detection','Presence')
  
  return(sigmoids)
}

#' Simulate presence-absence data
#'
#' Simulate presence-absence data
#' based on presence probabilities.
#'
#' @param pre numeric array: presence probabilities
#'
#' @return binary array
rpresence <- function(pre){
  yrep <- array(NA, dim = c(dim(pre)[1:3], 1),
                dimnames = list(dimnames(pre)[[1]],
                                dimnames(pre)[[2]],
                                dimnames(pre)[[3]],
                                'Response')
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

#' Posterior simulation
#'
#' Simulate presence-absences based on a \code{btvoccu} model.
#' Output some predictive checks and statistics,
#' including a histogram.
#'
#' @param model \code{btvoccu} model object
#' @param x numeric array: covariate array
#' @param y binary array: response array (default \code{NULL} does not provide posterior predictive checks)
#' @param sites vector (default \code{NULL} uses training sites)
#' @param seasons vector (default \code{NULL} uses training seasons)
#' @param periods vector (default \code{NULL} uses training periods)
#' @param nrep number of replicates
#' @param ndraws number of posterior draws
#' @param burnin percent of posterior samples to burn
#' @param thin keep every nth draw
#' 
#' @return list
#'
#' @export
posterior_check <- function(model, 
                            x, y = NULL,
                            sites = NULL,
                            seasons = NULL,
                            periods = NULL,
                            nrep = 1000,
                            ndraws = 1000,
                            burnin = .5,
                            thin = 1){
                               
  # get data indices
  if(is.null(sites)){
    sites <- model$sites
  }
  if(is.null(seasons)){
    seasons <- model$seasons
  }
  if(is.null(periods)){
    periods <- model$periods
  }
  
  # sort because posterior_sigmoids sorts
  sites <- sort(sites)
  seasons <- sort(seasons)
  periods <- sort(periods)
  
  # subset data
  x <- subset_4darray(x, 1, sites)
  x <- subset_4darray(x, 2, seasons)
  x <- subset_4darray(x, 3, periods)
  
  occueffs <- model$occueffs
  deteffs <- model$deteffs
  speffs <- model$speffs
  
  X <- subset_4darray(x, 4, occueffs)
  W <- subset_4darray(x, 4, deteffs)
  if(is.null(speffs)){
    M <- NULL
  } else{
    M <- subset_4darray(x, 4, speffs)
  }
  
  if(is.null(M)){
    x <- abind::abind(X, W)
  } else{
    x <- abind::abind(X, W, M)
  }
  
  # occupancy and detection designs can have overlapping covariates
  x <- x[,,, which(!duplicated(dimnames(x)[[4]]), arr.ind = T), drop = F]
  
  # access posterior sigmoid probabilities
  sigmoids <- posterior_sigmoids(model, x, sites, seasons, periods, ndraws, burnin, thin)
  if(!is.null(y)){
    y <- subset_4darray(y, 1, sites)
    y <- subset_4darray(y, 2, seasons)
    y <- subset_4darray(y, 3, periods)
    for(i in 1:length(sites)){
      for(j in 1:length(seasons)){
        for(k in 1:length(periods)){
          if(is.na(y[i,j,k,1])){
            sigmoids$Presence[i,j,k,] <- NA
          }
        }
      }
    }
  }
  
  # simulate presence-absence data
  yreppre <- rep(NA, nrep)
  for(n in 1:nrep){
    yrep <- rpresence(sigmoids$Presence)
    yreppre[n] <- sum(apply(yrep, 1, sum, na.rm = T))
  }
  
  # outputs
  output <- list()
  if(!is.null(y)){
    ypre <- sum(apply(y, 1, sum, na.rm = T))
    
    # plotting
    output[[1]] <- hist(yreppre,
                        main = NULL,
                        ylab = 'Density',
                        xlab = 'Presences',
                        probability = T,
                        xlim = c(min(ypre, min(yreppre)),
                                 max(ypre, max(yreppre))),
                        breaks <- seq(min(ypre, min(yreppre)) - .5,
                                      max(ypre, max(yreppre)) + .5),
                        col = 'gray')
    abline(v = ypre, col = 'red')
    
    # checks
    prepval <- sum(yreppre <  (ypre + 0.5)) / nrep
    prerel <- 1 + (mean(yreppre) - ypre) / ypre
    premse <- sum((ypre - yreppre) ** 2) / length(yreppre)
    
    output[[2]] <- prepval
    output[[3]] <- prerel
    output[[4]] <- premse
    output[[5]] <- mean(yreppre)
    output[[6]] <- sd(yreppre)
    output[[7]] <- ypre
    names(output) <- c('histogram','pval','relativepresence','mse','avgpresences', 'sdpresences','actualpresences')
  } else{
    output[[1]] <- hist(yreppre,
                        main = NULL,
                        ylab = 'Density',
                        xlab = 'Presences',
                        probability = T,
                        xlim = c(min(yreppre),
                                 max(yreppre)),
                        breaks <- seq(min(yreppre) - .5,
                                      max(yreppre) + .5),
                        col = 'gray'
    )
    output[[2]] <- NULL
    output[[3]] <- NULL
    output[[4]] <- NULL
    output[[5]] <- mean(yreppre)
    output[[6]] <- sd(yreppre)
    output[[7]] <- NULL
    names(output) <- c('histogram',
                       'pval',
                       'relativepresences',
                       'mse',
                       'avgpresences',
                       'sdpresences',
                       'actualpresences')
  }
  
  return(output)
}

# Plotting and predicting -------------------------------------------------
 
#' \code{btvoccu} predictions
#'
#' Predict occupancy, detection, and presence
#' probabilities with uncertainty over time.
#'
#' @param model \code{btvoccu} model object
#' @param x numeric array: predictors array
#' @param sites character vector
#' @param season integer
#' @param periods integer vector
#' @param value 'Occupancy', 'Detection', or 'Presence'
#' @param burnin percent of posterior samples to burn
#' @param thin keep every nth posterior draw
#' @param credible quantiles (length \code{3})
#' @param ndraws number of posterior draws
#' @param spline logical: default \code{TRUE} to fit smoothing splines
#' @param spline.predict logical: default \code{FALSE} to not predict missing values using splines
#' @param nknots see \code{smooth.spline()}
#' 
#' @return data frame
#'
#' @export
btvoccu_predict <- function(model, x,
                            sites,
                            season,
                            periods,
                            value = 'Occupancy',
                            burnin = .5,
                            thin = 1,
                            credible = c(.025,.5,.975),
                            ndraws = 1000,
                            spline = TRUE,
                            spline.predict = FALSE,
                            nknots = 5){
  
  # sort because posterior_sigmoids sorts
  sites <- sort(sites)
  periods <- sort(periods)
  
  # access posterior sigmoid probabilities
  sigmoids <- posterior_sigmoids(model, x, sites, season, periods, ndraws, burnin, thin)
  sigmoid <- sigmoids[[which(names(sigmoids) %in% value, arr.ind = T)]]

  # set up data frame
  n <- 0
  df <- data.frame(matrix(NA, ncol = 6, nrow = 0))
  colnames(df) <- c('Site','Season','Period','Lower','Middle','Upper')
  for(i in 1:length(sites)){
    for(k in length(periods):1){ # iterate in reverse for weeks
      n <- n + 1
      df[n,] <- c(sites[i], season, periods[k],
                  quantile(sigmoid[i,1,k,], probs = credible, na.rm = T))
    }
  }
  df$Period <- as.numeric(df$Period)
  df$Lower <- as.numeric(df$Lower)
  df$Middle <- as.numeric(df$Middle)
  df$Upper <- as.numeric(df$Upper)
  
  # splining
  if(spline){
    for(i in 1:length(sites)){
      sub <- subset(df, Site == sites[i])
      sub <- sub[order(sub$Period),]
      indices <- as.numeric(row.names(sub))
      sub.nna <- sub[!is.na(sub$Middle),] # nna for not NA
      indices.nna <- as.numeric(row.names(sub.nna))
      
      # splining requires at least 4 points
      if(nrow(sub.nna) >= 4){
        # Lower
        spl <- smooth.spline(sub.nna$Period, sub.nna$Lower, nknots = 5)
        if(spline.predict){
          df$Lower[indices] <- predict(spl, sub$Period)$y
        } else{
          df$Lower[indices.nna] <- spl$y
        }
        # Middle
        spl <- smooth.spline(sub.nna$Period, sub.nna$Middle, nknots = 5)
        if(spline.predict){
          df$Middle[indices] <- predict(spl, sub$Period)$y
        } else{
          df$Middle[indices.nna] <- spl$y
        }
        # Upper
        spl <- smooth.spline(sub.nna$Period, sub.nna$Upper, nknots = 5)
        if(spline.predict){
          df$Upper[indices] <- predict(spl, sub$Period)$y
        } else{
          df$Upper[indices.nna] <- spl$y
        }
      }
    }
  }
  
  # correct splining
  # splines sometimes give probabilities > 1 or < 0
  gindices <- which(df$Lower > 1, arr.ind = T)
  df$Lower[gindices] <- 1
  lindices <- which(df$Lower < 0, arr.ind = T)
  df$Lower[lindices] <- 0
  gindices <- which(df$Middle > 1, arr.ind = T)
  df$Middle[gindices] <- 1
  lindices <- which(df$Middle < 0, arr.ind = T)
  df$Middle[lindices] <- 0
  gindices <- which(df$Upper > 1, arr.ind = T)
  df$Upper[gindices] <- 1
  lindices <- which(df$Upper < 0, arr.ind = T)
  df$Upper[lindices] <- 0
  
  return(df)
}

#' \code{btvoccu} predictions visualized
#'
#' Visualize occupancy, detection, and presence
#' probabilities with uncertainty over time.
#'
#' @param model \code{btvoccu} model object
#' @param x numeric array: covariate array
#' @param sites vector
#' @param season length \code{1} vector
#' @param periods vector
#' @param value 'Occupancy', 'Detection', or 'Presence'
#' @param burnin percent of posterior samples to burn
#' @param thin keep every nth posterior draw
#' @param credible quantiles (length \code{3})
#' @param ndraws number of posterior draws
#' @param ribbons logical: default \code{TRUE} to show uncertainty ribbons
#' @param dots logical: default \code{FALSE} to not show actual datapoints
#' @param spline logical: default \code{FALSE} to not fit smoothing splines
#' @param spline.predict logical: default \code{FALSE} to not predict missing values using splines
#' @param nknots see \code{smooth.spline()}
#' @param xaxis character: x-axis label
#' @param legendtile character: legend label
#'
#' @return \code{ggplot} object
#'
#' @export
plot_btvoccu <- function(model,
                         x,
                         sites,
                         season,
                         periods,
                         value = 'Occupancy',
                         burnin = .5,
                         thin = 1,
                         credible = c(.025,.5,.975),
                         ndraws = 1000,
                         ribbons = TRUE,
                         dots = FALSE,
                         spline = FALSE,
                         spline.predict = FALSE,
                         nknots = 5,
                         xaxis = 'Period',
                         legendtitle = 'Site'){
  
  df <- btvoccu_predict(model,
                        x,
                        sites,
                        season,
                        periods,
                        value,
                        burnin,
                        thin,
                        credible,
                        ndraws,
                        spline,
                        spline.predict,
                        nknots)
  # plotting
  plt <- ggplot2::ggplot(df, ggplot2::aes(Period, Middle, colour = factor(Site), fill = factor(Site)))
  if(dots){
    if(ribbons){
      return(
        plt +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 1, show.legend = F) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper), alpha = .1, show.legend = F) +
          ggplot2::theme_classic() +
          ggplot2::labs(x = xaxis, y = value, colour = legendtitle) +
          ggplot2::scale_y_continuous(limits = c(0, 1)) # robust to splining
      )
    } else{
      return(
        plt +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 1, show.legend = F) +
          ggplot2::theme_classic() +
          ggplot2::labs(x = xaxis, y = value, colour = legendtitle) +
          ggplot2::scale_y_continuous(limits = c(0, 1)) # robust to splining
      )
    }
  } else{
    if(ribbons){
      return(
        plt +
          ggplot2::geom_line() +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper), alpha = .1, show.legend = F) +
          ggplot2::theme_classic() +
          ggplot2::labs(x = xaxis, y = value, colour = legendtitle) +
          ggplot2::scale_y_continuous(limits = c(0, 1)) # robust to splining
      )
    } else{
      return(
        plt +
          ggplot2::geom_line() +
          ggplot2::theme_classic() +
          ggplot2::labs(x = xaxis, y = value, colour = legendtitle) +
          ggplot2::scale_y_continuous(limits = c(0, 1)) # robust to splining
      )
    }
  }
}

#' Covariate time series plot
#'
#' Visualize covariates over time.
#'
#' @param x numeric array: covariate array
#' @param sites vector
#' @param season length \code{1} vector
#' @param periods vector
#' @param covariate character
#' @param yaxis characer: y-axis label
#' @param xaxis character: x-axis label
#' @param legendtile character: legend label
#'
#' @return \code{ggplot} object
#'
#' @export
plot_covariate <- function(x,
                           sites,
                           season,
                           periods,
                           covariate,
                           yaxis,
                           xaxis = 'Period',
                           legendtitle = 'Site'){
                              
  sites <- sort(sites)
  periods <- sort(periods)
  
  x <- subset_4darray(x, 3, periods)
  x <- subset_4darray(x, 4, covariate)
  mx <- max(x, na.rm = T)
  mn <- min(x, na.rm = F)
  
  x <- subset_4darray(x, 1, sites)
  x <- subset_4darray(x, 2, season)
  
  n <- 0
  df <- data.frame(matrix(NA, ncol = 4, nrow = 0))
  colnames(df) <- c('Site','Season','Period','Covariate')
  for(i in 1:length(sites)){
    for(k in length(periods):1){ # iterate in reverse for weeks
      n <- n + 1
      df[n,] <- c(sites[i], season, periods[k], x[i,1,k,1])
    }
  }
  
  df$Covariate <- as.numeric(as.character(df$Covariate))
  df$Period <- as.numeric(as.character(df$Period))
  
  # plotting
  plt <- ggplot2::ggplot(df, ggplot2::aes(Period, Covariate, colour = factor(Site)))
  return(
    plt +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 1) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = xaxis, y = yaxis, colour = legendtitle) +
      ggplot2::scale_y_continuous(limits = c(mn - .001, mx + .001))
  )
}










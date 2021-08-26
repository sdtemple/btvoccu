#' Posterior Sigmoids
#'
#' Compute posterior occupancy, detection,
#' and presence probabilities.
#'
#' @param model \code{btvoccu} output
#' @param x numeric array: predictors array
#' @param sites vector
#' @param seasons vector
#' @param periods vector
#' @param ndraws integer: # of posterior draws
#' @param burnin percent of posterior samples to burn
#' @param thin integer: default \code{1} implies no thinning, \code{2} means keep every other draw, etc.
#'
#' @return numeric array
#'
#' @export
posterior_sigmoids <- function(model,
                               x,
                               sites,
                               seasons,
                               periods,
                               ndraws = 1000,
                               burnin = .5,
                               thin = 1
                               ){

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

  keep <- seq(ceiling(burnin * model$niter),
              model$niter,
              by = thin)

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

  occu <- array(NA,
                dim = c(dim(x)[1:3], ndraws),
                dimnames = list(sites, seasons, periods, 1:ndraws))

  det <- array(NA,
               dim = c(dim(x)[1:3], ndraws),
               dimnames = list(sites, seasons, periods, 1:ndraws))

  pre <- array(NA,
               dim = c(dim(x)[1:3], ndraws),
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

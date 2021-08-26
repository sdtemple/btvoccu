#' Watanabe Akaike Information Criterion
#'
#' Compute Watanabe-Akaike Information Criterion.
#' When indices are set to \code{NULL},
#' WAIC is calculated using the training data.
#'
#' @param model \code{btvoccu} output
#' @param x numeric array: predictors array
#' @param y binary array: response array
#' @param sites vector (default \code{NULL} uses training sites)
#' @param seasons vector (default \code{NULL} uses training seasons)
#' @param periods vector (default \code{NULL} uses training periods)
#' @param burnin percent of posterior samples to burn
#' @param thin integer: default \code{1} implies no thinning, \code{2} means keep every other draw, etc.
#'
#' @return list
#'
#' @export
waic_score <- function(model,
                       x,
                       y,
                       sites = NULL,
                       seasons = NULL,
                       periods = NULL,
                       burnin = .5,
                       thin = 1
                       ){

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
                dim = c(dim(x)[1:3], nsample),
                dimnames = list(sites, seasons, periods, 1:nsample))

  det <- array(NA,
               dim = c(dim(x)[1:3], nsample),
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

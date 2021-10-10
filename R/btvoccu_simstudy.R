#' Simulate time-varying occupancy model
#'
#' This function simulates data for a time-varying occupancy model.
#'
#' @param nsites number of sites
#' @param nseasons number of seasons
#' @param nperiods number of periods
#' @param beta vector of occupancy effects
#' @param alpha vector of detection effects
#' @param occusite means for site-specific occupancy covariates
#' @param detsite means for site-specific detection covariates
#' @param occutime array of time-varying occupancy covariates
#' @param dettime array of time-varying detection covariates
#' @param naprob percent of NAs in simulated data
#' @param noise numeric parameter controlling covariate noise
#' @param probit logical: default \code{F} for logit link and \code{T} for probit link
#'
#' @return list including simulated data and more
#'
#' @export
btvoccu_simstudy <- function(nsites, 
                             nseasons, 
                             nperiods,
                             beta, 
                             alpha,
                             occusite,
                             detsite,
                             occutime,
                             dettime,
                             naprob = 0, 
                             noise = 1/4,
                             probit = F){
  
  # setup
  rbt <- nrow(occutime)
  rat <- nrow(dettime)
  lb <- length(beta) - rbt - 1
  la <- length(alpha) - rat - 1
  
  # data input checks
  if(lb != length(occusite)){print("occusite has wrong number of means")}
  if(la != length(detsite)){print("detsite has wrong number of means")}
  if(nperiods != ncol(occutime)){print("occutime must have as many columns as nperiods"); return()}
  if(nperiods != ncol(dettime)){print("dettime must have as many columns as nperiods"); return()}
  if(lb < 0){"occutime has too many rows"; return()}
  if(la < 0){"dettime has too many rows"; return()}
  if((naprob > 1) | (naprob < 0)){"naprob must be in [0,1]"; return()}

  # site-occupancy
  psi <- array(dim = c(nsites, nseasons, nperiods))
  X <- array(dim = c(nsites, nseasons, nperiods, length(beta)))
  X[,,,1] <- 1 # intercept
  for(i in 1:nsites){
    X[i,,,2:(lb+1)] <- mvtnorm::rmvnorm(1, occusite, abs(occusite) / 4 * diag(lb)) # site-specific
    for(j in 1:nseasons){
      for(l in 1:rbt){X[i,j,,(lb+1+l)] <- occutime[l,] + rnorm(nperiods, sd = noise / 10)} # time-varying
      for(k in 1:nperiods){
        if(probit){
          psi[i,j,k] <- pnorm(beta %*% X[i,j,k,] + rnorm(1, sd = noise))
        } else{
          psi[i,j,k] <- faraway::ilogit(beta %*% X[i,j,k,] + rnorm(1, sd = noise))
        }
      }
    }
  }

  # detection
  p <- array(dim = c(nsites, nseasons, nperiods))
  W <- array(dim = c(nsites, nseasons, nperiods, length(alpha)))
  W[,,,1] <- 1 # intercept
  for(i in 1:nsites){
    W[i,,,2:(la+1)] <- mvtnorm::rmvnorm(1, detsite, abs(detsite) / 4 * diag(la)) # site-specific
    for(j in 1:nseasons){
      for(l in 1:rat){W[i,j,,(la+1+l)] <- dettime[l,] + rnorm(nperiods, sd = noise / 10)} # time-varying
      for(k in 1:nperiods){
        if(probit){
          p[i,j,k] <- pnorm(alpha %*% W[i,j,k,] + rnorm(1, sd = noise))
        } else{
          p[i,j,k] <- faraway::ilogit(alpha %*% W[i,j,k,] + rnorm(1, sd = noise))
        }
      }
    }
  }

  # presence-absence
  y <- array(dim=c(nsites, nseasons, nperiods, 1))
  for(i in 1:nsites){
    for(j in 1:nseasons){
      for(k in 1:nperiods){
        y[i,j,k, 1] <- rbinom(1, 1, rbinom(1, 1, psi[i,j,k]) * p[i,j,k])
      }
    }
  }

  # insert NAs
  for(i in 1:nsites){
    for(j in 1:nseasons){
      for(k in 1:nperiods){
        if(runif(1) < naprob){
          y[i,j,k,1] <- NA
        }
      }
    }
  }
  na <- !is.na(y)
  na <- apply(na, c(1,2,3,4), as.numeric)
  yna <- replace(y, is.na(y), 0)
  
  # naming
  sites <- 1:nsites
  dimnames(y)[[1]] <- sites
  dimnames(X)[[1]] <- sites
  dimnames(W)[[1]] <- sites
  
  seasons <- 1:nseasons
  dimnames(y)[[2]] <- seasons
  dimnames(X)[[2]] <- seasons
  dimnames(W)[[2]] <- seasons
  
  periods <- 1:nperiods
  dimnames(y)[[3]] <- periods
  dimnames(X)[[3]] <- periods
  dimnames(W)[[3]] <- periods
  
  occueffs <- paste("bvar", 1:length(beta), sep = "")
  deteffs <- paste("avar", 1:length(alpha), sep = "")
  dimnames(y)[[4]] <- "response"
  dimnames(X)[[4]] <- occueffs
  dimnames(W)[[4]] <- deteffs

  return(list(y=y, 
              yna=yna, 
              na=na, 
              X=X, 
              W=W,
              p=p,
              psi=psi,
              nsites=nsites, 
              nseasons=nseasons, 
              nperiods=nperiods,
              beta=beta, 
              alpha=alpha,
              sites=sites,
              seasons=seasons,
              periods=periods,
              occueffs=occueffs,
              deteffs=deteffs))
}

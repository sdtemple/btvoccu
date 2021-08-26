#' Bayesian time-varying occupancy model
#'
#' Train a time-varying occupancy model
#' using Bayesian generalized linear regression.
#'
#' @param niter integer: # of iterations
#' @param x numeric array: covariate array
#' @param y binary array: response array
#' @param sites vector
#' @param seasons vector
#' @param periods vector
#' @param occueffs character vector: occupancy effects
#' @param deteffs character vector: detection effects
#' @param speffs character vector: spatial random effects
#' @param A adjacency matrix
#' @param logit logical: default \code{TRUE} for logistic regression and \code{FALSE} for probit regression
#' @param nchains number of chains (default \code{1})
#' @param print.interval integer: period to print iteration to console (default \code{1000})
#'
#' @return length \code{15} list of posterior samples and more
#'
#' @export
btvoccu <- function(niter, 
                    x, y,
                    sites, 
                    seasons, 
                    periods,
                    occueffs, 
                    deteffs, 
                    speffs = NULL, 
                    A = NULL,
                    logit = TRUE, 
                    nchains = 1,
                    print.interval = 1000){

  # data processing
  y <- subset_4darray(y, 1, sites)
  y <- subset_4darray(y, 2, seasons)
  y <- subset_4darray(y, 3, periods)
  x <- subset_4darray(x, 1, sites)
  x <- subset_4darray(x, 2, seasons)
  x <- subset_4darray(x, 3, periods)
  X <- subset_4darray(x, 4, occueffs) # occupancy design
  W <- subset_4darray(x, 4, deteffs) # detection design

  # spatial random effects design
  if(is.null(speffs)){
    M <- NULL
  } else{
    M <- subset_4darray(x, 4, speffs)
  }

  # local functions
  iw_update <- function(beta, nu, Psi){return(MCMCpack::riwish(1 + nu, Psi + tcrossprod(beta, beta)))}
  iga_update <- function(theta, Sigma, a, b){
    return(MCMCpack::rinvgamma(1, a + length(theta) / 2, b + crossprod(theta, Sigma) %*% theta / 2))}

  # models
  if(is.null(M)){ # without spatial random effect
    
    # samplers
    gibbs_logit_btvoccu <- function(niter,
                                    y,
                                    X,
                                    W,
                                    mubeta,
                                    nubeta,
                                    Psibeta,
                                    mualpha,
                                    nualpha,
                                    Psialpha,
                                    print.interval = 1000){

      # local functions
      mvn_update <- function(y, X, beta, mubeta, Sigmabeta){
        # sample augmented variable
        ya <- rep(NA, length(y))
        for(i in 1:length(y)){
          ya[i] <- BayesLogit::rpg(1, z = X[i,] %*% beta)
        }

        # matrix algebra
        kappa <- y - 1/2
        V <- solve(solve(Sigmabeta) + crossprod(X, diag(ya) %*% X))
        m <- V %*% (solve(Sigmabeta) %*% mubeta  + crossprod(X, kappa))

        return(mvtnorm::rmvnorm(1, mean = m, sigma = V))
      }

      z_update <- function(y, X, W, beta, alpha){
        z <- rep(NA, length(y))
        for(i in 1:length(y)){
          if(y[i]){
            z[i] <- 1
          } else{
            lognum <- log(faraway::ilogit(crossprod(X[i,], beta))) + log(1 - faraway::ilogit(crossprod(W[i,], alpha)))
            logden <- log(exp(lognum) + (1 - faraway::ilogit(crossprod(X[i,], beta))))
            z[i] <- rbinom(1, 1, exp(lognum - logden))
          }
        }

        return(z)
      }

      # setup matrices to record posterior draws
      betas <- matrix(NA, nrow = niter, ncol = length(mubeta))
      betas[1,] <- mubeta
      alphas <- matrix(NA, nrow = niter, ncol = length(mualpha))
      alphas[1,] <- mualpha
      Sigmabetas <- array(NA, dim = c(niter, length(mubeta), length(mubeta)))
      Sigmabetas[1,,] <- iw_update(mubeta, nubeta, Psibeta)
      Sigmaalphas <- array(NA, dim = c(niter, length(mualpha), length(mualpha)))
      Sigmaalphas[1,,] <- iw_update(mualpha, nualpha, Psialpha)
      
      # flatten arrays, removing NAs from the start
      ynna <- !is.na(y)
      nobs <- sum(ynna)
      yobs <- rep(NA, nobs)
      Xobs <- array(dim = c(nobs, dim(X)[4]))
      Wobs <- array(dim = c(nobs, dim(W)[4]))
      n <- 0
      for(i in 1:dim(y)[1]){
        for(j in 1:dim(y)[2]){
          for(k in 1:dim(y)[3]){
            if(ynna[i,j,k,1]){
              n <- n + 1
              yobs[n] <- y[i,j,k,1]
              Xobs[n,] <- X[i,j,k,]
              Wobs[n,] <- W[i,j,k,]
            }
          }
        }
      }

      # gibbs sampling
      for(n in 2:niter){
        if((n %% print.interval) == 0){print(n)}
        z <- z_update(yobs, Xobs, Wobs, betas[(n-1),], alphas[(n-1),])
        betas[n,] <- mvn_update(z, Xobs, betas[(n-1),], mubeta, Sigmabetas[(n-1),,])
        Sigmabetas[n,,] <- iw_update(betas[n,], nubeta, Psibeta)

        # subset to occupied sites
        z1 <- which(z == 1, arr.ind = T)
        y1 <- yobs[z1]
        W1 <- Wobs[z1,, drop = F]
        alphas[n,] <- mvn_update(y1, W1, alphas[(n-1),], mualpha, Sigmaalphas[(n-1),,])
        Sigmaalphas[n,,] <- iw_update(alphas[n,], nualpha,Psialpha)
      }
      
      return(list(betas = betas, alphas = alphas, Sigmabetas = Sigmabetas, Sigmaalphas = Sigmaalphas))
    }

    gibbs_probit_btvoccu <- function(niter,
                                     y,
                                     X,
                                     W,
                                     mubeta,
                                     nubeta,
                                     Psibeta,
                                     mualpha,
                                     nualpha,
                                     Psialpha,
                                     print.interval = 1000){

      # local functions
      mvn_update <- function(y, X, beta, mubeta, Sigmabeta, XtX){
        # sample augmented variable
        ya <- rep(NA, length(y))
        for(i in 1:length(y)){
          if(y[i]){
            ya[i] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = X[i,] %*% beta, sd = 1)
          } else{
            ya[i] <- truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = X[i,] %*% beta, sd = 1)
          }
        }
        
        # matrix algebra
        V <- solve(solve(Sigmabeta) + XtX)
        m <- V %*% (solve(Sigmabeta) %*% mubeta  + crossprod(X, ya))
        
        return(mvtnorm::rmvnorm(1, mean = m, sigma = V))
      }

      z_update <- function(y, X, W, beta, alpha){
        z <- rep(NA, length(y))
        for(i in 1:length(y)){
          if(y[i]){
            z[i] <- 1
          } else{
            lognum <- pnorm(crossprod(X[i,], beta), log.p = T) + log(1 - pnorm(crossprod(W[i,], alpha)))
            logden <- log(exp(lognum) + (1 - pnorm(crossprod(X[i,], beta))))
            z[i] <- rbinom(1, 1, exp(lognum - logden))
          }
        }
        return(z)
      }

      # setup matrices to record posterior draws
      betas <- matrix(NA, nrow = niter, ncol = length(mubeta))
      betas[1,] <- mubeta
      alphas <- matrix(NA, nrow = niter, ncol = length(mualpha))
      alphas[1,] <- mualpha
      Sigmabetas <- array(NA, dim = c(niter, length(mubeta), length(mubeta)))
      Sigmabetas[1,,] <- iw_update(mubeta, nubeta, Psibeta)
      Sigmaalphas <- array(NA, dim = c(niter, length(mualpha), length(mualpha)))
      Sigmaalphas[1,,] <- iw_update(mualpha, nualpha, Psialpha)

      # flatten arrays, removing NAs from the start
      ynna <- !is.na(y)
      nobs <- sum(ynna)
      yobs <- rep(NA, nobs)
      Xobs <- array(dim = c(nobs, dim(X)[4]))
      Wobs <- array(dim = c(nobs, dim(W)[4]))
      n <- 0
      for(i in 1:dim(y)[1]){
        for(j in 1:dim(y)[2]){
          for(k in 1:dim(y)[3]){
            if(ynna[i,j,k,1]){
              n <- n + 1
              yobs[n] <- y[i,j,k,1]
              Xobs[n,] <- X[i,j,k,]
              Wobs[n,] <- W[i,j,k,]
            }
          }
        }
      }
      XtX <- crossprod(Xobs, Xobs)

      # gibbs sampling
      for(n in 2:niter){
        if((n %% print.interval) == 0){print(n)}
        z <- z_update(yobs, Xobs, Wobs, betas[(n-1),], alphas[(n-1),])
        betas[n,] <- mvn_update(z, Xobs, betas[(n-1),], mubeta, Sigmabetas[(n-1),,], XtX)
        Sigmabetas[n,,] <- iw_update(betas[n,], nubeta, Psibeta)

        # subset to occupied sites
        z1 <- which(z == 1, arr.ind = T)
        y1 <- yobs[z1]
        W1 <- Wobs[z1,, drop = F]
        alphas[n,] <- mvn_update(y1, W1, alphas[(n-1),], mualpha, Sigmaalphas[(n-1),,], crossprod(W1, W1))
        Sigmaalphas[n,,] <- iw_update(alphas[n,], nualpha, Psialpha)
      }
      
      return(list(betas = betas, alphas = alphas, Sigmabetas = Sigmabetas, Sigmaalphas = Sigmaalphas))
    }

    # choose sampler
    if(logit){
      sampler <- gibbs_logit_btvoccu
    } else{
      sampler <- gibbs_probit_btvoccu
    }

    # priors
    mubeta <- rep(0, dim(X)[4])
    nubeta <- dim(X)[4] + 1
    Psibeta <- diag(dim(X)[4])
    mualpha <- rep(0, dim(W)[4])
    nualpha <- dim(W)[4] + 1
    Psialpha <- diag(dim(W)[4])

    # sampling
    output <- list()
    betas <- array(NA, dim = c(nchains, niter, length(mubeta)))
    alphas <- array(NA, dim = c(nchains, niter, length(mualpha)))
    Sigmabetas <- array(NA, dim = c(nchains, niter, length(mubeta), length(mubeta)))
    Sigmaalphas <- array(NA, dim = c(nchains, niter, length(mualpha), length(mualpha)))
    for(n in 1:nchains){
      print(paste('Chain', n))
      sample <- sampler(niter, y, X, W, mubeta, nubeta, Psibeta, mualpha, nualpha, Psialpha, print.interval)
      betas[n,,] <- sample$betas
      alphas[n,,] <- sample$alphas
      Sigmabetas[n,,,] <- sample$Sigmabetas
      Sigmaalphas[n,,,] <- sample$Sigmaalphas
    }

    # formatting
    output[[1]] <- betas
    output[[2]] <- alphas
    output[[3]] <- NULL
    output[[4]] <- Sigmabetas
    output[[5]] <- Sigmaalphas
    output[[6]] <- NULL
    output[[7]] <- dimnames(X)[[4]]
    output[[8]] <- dimnames(W)[[4]]
    output[[9]] <- NULL
    output[[10]] <- sites
    output[[11]] <- seasons
    output[[12]] <- periods
    output[[13]] <- niter
    output[[14]] <- nchains
    if(logit){
      output[[15]] <- faraway::ilogit
    } else{
      output[[15]] <- pnorm
    }
    names(output) <- c('betas',
                       'alphas',
                       'thetas',
                       'Sigmabetas',
                       'Sigmaalphas',
                       'sigmathetas',
                       'occueffs',
                       'deteffs',
                       'speffs',
                       'sites',
                       'seasons',
                       'periods',
                       'niter',
                       'nchains',
                       'link')
    
    return(output)
  } else{ # with spatial random effects
    if(is.null(A)){
      stop('Provide an adjacency matrix.')
    }
    Qs <- crossprod(M[,1,1,], icarQ(A)) %*% M[,1,1,]
    
    # samplers
    gibbs_logit_btvspoccu <- function(niter,
                                      y,
                                      X,
                                      M,
                                      W,
                                      mubeta,
                                      nubeta,
                                      Psibeta,
                                      mutheta,
                                      Qs,
                                      atheta,
                                      btheta,
                                      mualpha,
                                      nualpha,
                                      Psialpha,
                                      print.interval = 1000){

      # local functions
      sp_mvn_update <- function(y, X, beta, mubeta, Sigmabeta, M, theta, mutheta, Sigmatheta){
        # sample augmented variable
        ya <- rep(NA, length(y))
        for(i in 1:length(y)){
          ya[i] <- BayesLogit::rpg(1, z = X[i,] %*% beta + M[i,] %*% theta)
        }

        # matrix algebra
        kappa <- y - 1/2
        Vbeta <- solve(solve(Sigmabeta) + crossprod(X, diag(ya) %*% X))
        mbeta <- Vbeta %*% (solve(Sigmabeta) %*% mubeta  + 
                              crossprod(X, (kappa - (diag(ya) %*% M %*% theta))))
        Vtheta <- solve(solve(Sigmatheta) + crossprod(M, diag(ya) %*% M))
        mtheta <- Vtheta %*% (solve(Sigmatheta) %*% mutheta +
                                crossprod(M, (kappa - (diag(ya) %*% X %*% beta))))

        return(list(beta = mvtnorm::rmvnorm(1, mean = mbeta, sigma = Vbeta),
                    theta = mvtnorm::rmvnorm(1, mean = mtheta, sigma = Vtheta)))
      }

      mvn_update <- function(y, X, beta, mubeta, Sigmabeta){
        # sample augmented variable
        ya <- rep(NA, length(y))
        for(i in 1:length(y)){
          ya[i] <- BayesLogit::rpg(1, z = X[i,] %*% beta)
        }

        # matrix algebra
        kappa <- y - 1/2
        V <- solve(solve(Sigmabeta) +
                     crossprod(X, diag(ya) %*% X))
        m <- V %*% (solve(Sigmabeta) %*% mubeta  +
                      crossprod(X, kappa))

        return(mvtnorm::rmvnorm(1, mean = m, sigma = V))
      }

      z_update <- function(y, X, M, W, beta, theta, alpha){
        z <- rep(NA, length(y))
        for(i in 1:length(y)){
          if(y[i]){
            z[i] <- 1
          } else{
            lognum <- log(faraway::ilogit(crossprod(X[i,], beta)+ crossprod(M[i,] %*% theta))) +
              log(1 - faraway::ilogit(crossprod(W[i,], alpha)))
            logden <- log(exp(lognum) + (1 - faraway::ilogit(crossprod(X[i,], beta) + crossprod(M[i,], theta))))
            z[i] <- rbinom(1, 1, exp(lognum - logden))
          }
        }

        return(z)
      }

      # setup matrices to record posterior draws
      betas <- matrix(NA, nrow = niter, ncol = length(mubeta))
      betas[1,] <- mubeta
      thetas <- matrix(NA, nrow = niter, ncol = length(mutheta))
      thetas[1,] <- mutheta
      alphas <- matrix(NA, nrow = niter, ncol = length(mualpha))
      alphas[1,] <- mualpha
      Sigmabetas <- array(NA, dim = c(niter, length(mubeta), length(mubeta)))
      Sigmabetas[1,,] <- iw_update(mubeta, nubeta, Psibeta)
      sigmathetas <- array(NA, dim = c(niter, 1))
      sigmathetas[1,] <- 1
      Sigmaalphas <- array(NA, dim = c(niter, length(mualpha), length(mualpha)))
      Sigmaalphas[1,,] <- iw_update(mualpha, nualpha, Psialpha)

      # flatten arrays, removing NAs from the start
      ynna <- !is.na(y)
      nobs <- sum(ynna)
      yobs <- rep(NA, nobs)
      Xobs <- array(dim = c(nobs, dim(X)[4]))
      Mobs <- array(dim = c(nobs, dim(M)[4]))
      Wobs <- array(dim = c(nobs, dim(W)[4]))
      n <- 0
      for(i in 1:dim(y)[1]){
        for(j in 1:dim(y)[2]){
          for(k in 1:dim(y)[3]){
            if(ynna[i,j,k,1]){
              n <- n + 1
              yobs[n] <- y[i,j,k,1]
              Xobs[n,] <- X[i,j,k,]
              Mobs[n,] <- M[i,j,k,]
              Wobs[n,] <- W[i,j,k,]
            }
          }
        }
      }

      # gibbs sampling
      for(n in 2:niter){
        if((n %% print.interval) == 0){
          print(n)
        }

        z <- z_update(yobs, Xobs, Mobs, Wobs, betas[(n-1),], thetas[(n-1),], alphas[(n-1),])
        spdraw <- sp_mvn_update(z, Xobs, betas[(n-1),], mubeta, Sigmabetas[(n-1),,], Mobs, thetas[(n-1),], mutheta, sigmathetas[(n-1),] * solve(Qs))
        betas[n,] <- spdraw$beta
        Sigmabetas[n,,] <- iw_update(betas[n,], nubeta, Psibeta)
        thetas[n,] <- spdraw$theta
        sigmathetas[n,] <- iga_update(thetas[n,], solve(Qs), atheta, btheta)

        # subset to occupied sites
        z1 <- which(z == 1, arr.ind = T)
        y1 <- yobs[z1]
        W1 <- Wobs[z1,, drop = F]
        alphas[n,] <- mvn_update(y1, W1, alphas[(n-1),], mualpha, Sigmaalphas[(n-1),,])
        Sigmaalphas[n,,] <- iw_update(alphas[n,], nualpha, Psialpha)
      }
      return(list(betas = betas, thetas = thetas, alphas = alphas, Sigmabetas = Sigmabetas, sigmathetas = sigmathetas, Sigmaalphas = Sigmaalphas))
    }

    gibbs_probit_btvspoccu <- function(niter,
                                       y,
                                       X,
                                       M,
                                       W,
                                       mubeta,
                                       nubeta,
                                       Psibeta,
                                       mutheta,
                                       Qs,
                                       atheta,
                                       btheta,
                                       mualpha,
                                       nualpha,
                                       Psialpha,
                                       print.interval = 1000){

      # local functions
      sp_mvn_update <- function(y, X, beta, mubeta, Sigmabeta, XtX, M, theta, mutheta, Sigmatheta, MtM){
        # sample augmented variable
        ya <- rep(NA, length(y))
        for(i in 1:length(y)){
          if(y[i]){
            ya[i] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = X[i,] %*% beta + M[i,] %*% theta, sd = 1)
          } else{
            ya[i] <- truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = X[i,] %*% beta + M[i,] %*% theta, sd = 1)
          }
        }

        # matrix algebra
        Vbeta <- solve(solve(Sigmabeta) + XtX)
        mbeta <- Vbeta %*% (solve(Sigmabeta) %*% mubeta  + crossprod(X, (ya - (M %*% theta))))
        Vtheta <- solve(solve(Sigmatheta) + MtM)
        mtheta <- Vtheta %*% (solve(Sigmatheta) %*% mutheta + crossprod(M, (ya - (X %*% beta))))

        return(list(beta = mvtnorm::rmvnorm(1, mean = mbeta, sigma = Vbeta),
                    theta = mvtnorm::rmvnorm(1, mean = mtheta, sigma = Vtheta)))
      }


      mvn_update <- function(y, X, beta, mubeta, Sigmabeta, XtX){
        # sample augmented variable
        ya <- rep(NA, length(y))
        for(i in 1:length(y)){
          if(y[i]){
            ya[i] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = X[i,] %*% beta, sd = 1)
          } else{
            ya[i] <- truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = X[i,] %*% beta, sd = 1)
          }
        }

        # matrix algebra
        V <- solve(solve(Sigmabeta) + XtX)
        m <- V %*% (solve(Sigmabeta) %*% mubeta + crossprod(X, ya))

        return(mvtnorm::rmvnorm(1, mean = m, sigma = V))
      }

      z_update <- function(y, X, M, W, beta, theta, alpha){
        z <- rep(NA, length(y))
        for(i in 1:length(y)){
          if(y[i]){
            z[i] <- 1
          } else{
            lognum <- pnorm(crossprod(X[i,], beta) + crossprod(M[i,] %*% theta), log.p = T) +
              log(1 - pnorm(crossprod(W[i,], alpha)))
            logden <- log(exp(lognum) + (1 - pnorm(crossprod(X[i,], beta) + crossprod(M[i,], theta))))
            z[i] <- rbinom(1, 1, exp(lognum - logden))
          }
        }
        
        return(z)
      }

      # setup matrices to record posterior draws
      betas <- matrix(NA, nrow = niter, ncol = length(mubeta))
      betas[1,] <- mubeta
      thetas <- matrix(NA, nrow = niter, ncol = length(mutheta))
      thetas[1,] <- mutheta
      alphas <- matrix(NA, nrow = niter, ncol = length(mualpha))
      alphas[1,] <- mualpha
      Sigmabetas <- array(NA, dim = c(niter, length(mubeta), length(mubeta)))
      Sigmabetas[1,,] <- iw_update(mubeta, nubeta, Psibeta)
      sigmathetas <- matrix(NA, nrow = niter, ncol = 1)
      sigmathetas[1,] <- 1
      Sigmaalphas <- array(NA, dim = c(niter, length(mualpha), length(mualpha)))
      Sigmaalphas[1,,] <- iw_update(mualpha, nualpha, Psialpha)

      # flatten arrays, removing NAs from the start
      ynna <- !is.na(y)
      nobs <- sum(ynna)
      yobs <- rep(NA, nobs)
      Xobs <- array(dim = c(nobs, dim(X)[4]))
      Mobs <- array(dim = c(nobs, dim(M)[4]))
      Wobs <- array(dim = c(nobs, dim(W)[4]))
      n <- 0
      for(i in 1:dim(y)[1]){
        for(j in 1:dim(y)[2]){
          for(k in 1:dim(y)[3]){
            if(ynna[i,j,k,1]){
              n <- n + 1
              yobs[n] <- y[i,j,k,1]
              Xobs[n,] <- X[i,j,k,]
              Mobs[n,] <- M[i,j,k,]
              Wobs[n,] <- W[i,j,k,]
            }
          }
        }
      }
      XtX <- crossprod(Xobs, Xobs)
      MtM <- crossprod(Mobs, Mobs)

      # gibbs sampling
      for(n in 2:niter){
        if((n %% print.interval) == 0){print(n)}
        z <- z_update(yobs, Xobs, Mobs, Wobs, betas[(n-1),], thetas[(n-1),], alphas[(n-1),])
        spdraw <- sp_mvn_update(z, Xobs, betas[(n-1),], mubeta, Sigmabetas[(n-1),,], XtX, Mobs, thetas[(n-1),], mutheta, sigmathetas[(n-1),] * solve(Qs), MtM)
        betas[n,] <- spdraw$beta
        Sigmabetas[n,,] <- iw_update(betas[n,], nubeta, Psibeta)
        thetas[n,] <- spdraw$theta
        sigmathetas[n,] <- iga_update(thetas[n,], solve(Qs), atheta, btheta)

        # subset to occupied sites
        z1 <- which(z == 1, arr.ind = T)
        y1 <- yobs[z1]
        W1 <- Wobs[z1,, drop = F]
        alphas[n,] <- mvn_update(y1, W1, alphas[(n-1),], mualpha, Sigmaalphas[(n-1),,], crossprod(W1, W1))
        Sigmaalphas[n,,] <- iw_update(alphas[n,], nualpha, Psialpha)
      }

      return(list(betas = betas, thetas = thetas, alphas = alphas, Sigmabetas = Sigmabetas, sigmathetas = sigmathetas, Sigmaalphas = Sigmaalphas))
    }

    # choose sampler
    if(logit){
      sampler <- gibbs_logit_btvspoccu
    } else{
      sampler <- gibbs_probit_btvspoccu
    }

    # priors
    mubeta <- rep(0, dim(X)[4])
    nubeta <- dim(X)[4] + 1
    Psibeta <- diag(dim(X)[4])
    mutheta <- rep(0, dim(M)[4])
    atheta <- .001
    btheta <- .001
    mualpha <- rep(0, dim(W)[4])
    nualpha <- dim(W)[4] + 1
    Psialpha <- diag(dim(W)[4])


    # sampling
    output <- list()
    betas <- array(NA, dim = c(nchains, niter, length(mubeta)))
    thetas <- array(NA, dim = c(nchains, niter, length(mutheta)))
    alphas <- array(NA, dim = c(nchains, niter, length(mualpha)))
    Sigmabetas <- array(NA, dim = c(nchains, niter, length(mubeta), length(mubeta)))
    sigmathetas <- array(NA, dim = c(nchains, niter, 1))
    Sigmaalphas <- array(NA, dim = c(nchains, niter, length(mualpha), length(mualpha)))
    for(n in 1:nchains){
      print(paste('Chain', n))
      sample <- sampler(niter,
                        y,
                        X,
                        M,
                        W,
                        mubeta,
                        nubeta,
                        Psibeta,
                        mutheta,
                        Qs,
                        atheta,
                        btheta,
                        mualpha,
                        nualpha,
                        Psialpha,
                        print.interval)
      betas[n,,] <- sample$betas
      thetas[n,,] <- sample$thetas
      alphas[n,,] <- sample$alphas
      Sigmabetas[n,,,] <- sample$Sigmabetas
      sigmathetas[n,,] <- sample$sigmathetas
      Sigmaalphas[n,,,] <- sample$Sigmaalphas
    }

    # formatting
    output[[1]] <- betas
    output[[2]] <- alphas
    output[[3]] <- thetas
    output[[4]] <- Sigmabetas
    output[[5]] <- Sigmaalphas
    output[[6]] <- sigmathetas
    output[[7]] <- dimnames(X)[[4]]
    output[[8]] <- dimnames(W)[[4]]
    output[[9]] <- dimnames(M)[[4]]
    output[[10]] <- sites
    output[[11]] <- seasons
    output[[12]] <- periods
    output[[13]] <- niter
    output[[14]] <- nchains
    if(logit){
      output[[15]] <- faraway::ilogit
    } else{
      output[[15]] <- pnorm
    }
    names(output) <- c('betas',
                       'alphas',
                       'thetas',
                       'Sigmabetas',
                       'Sigmaalphas',
                       'sigmathetas',
                       'occueffs',
                       'deteffs',
                       'speffs',
                       'sites',
                       'seasons',
                       'periods',
                       'niter',
                       'nchains',
                       'link'
                       )

    return(output)
  }
}

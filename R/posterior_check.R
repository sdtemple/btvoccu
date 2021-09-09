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
#' @examples 
#' output <- posterior_check(mdl, x)
#' output <- posterior_check(mdl, x, NULL, c("TOR", "YRK", "PEE"), 2004:2005, 20:40)
#' output <- posterior_check(mdl, x, y, c("TOR", "YRK", "PEE"), 2004:2005, 20:40)
#'
#' @export
posterior_check <- function(model, 
                            x, y = NULL,
                            sites = NULL,
                            seasons = NULL,
                            periods = NULL,
                            nrep = 100,
                            ndraws = 200,
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
                        ylab = "Density",
                        xlab = "Presences",
                        probability = T,
                        xlim = c(min(ypre, min(yreppre)),
                                 max(ypre, max(yreppre))),
                        breaks <- seq(min(ypre, min(yreppre)) - .5,
                                      max(ypre, max(yreppre)) + .5),
                        col = "gray")
    abline(v = ypre, col = "red")
    
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
    names(output) <- c("histogram","pval","relativepresence","mse","avgpresences", "sdpresences","actualpresences")
  } else{
    output[[1]] <- hist(yreppre,
                        main = NULL,
                        ylab = "Density",
                        xlab = "Presences",
                        probability = T,
                        xlim = c(min(yreppre),
                                 max(yreppre)),
                        breaks <- seq(min(yreppre) - .5,
                                      max(yreppre) + .5),
                        col = "gray"
    )
    output[[2]] <- NULL
    output[[3]] <- NULL
    output[[4]] <- NULL
    output[[5]] <- mean(yreppre)
    output[[6]] <- sd(yreppre)
    output[[7]] <- NULL
    names(output) <- c("histogram",
                       "pValue",
                       "relativePresence",
                       "mse",
                       "avgPresence",
                       "sdPresence",
                       "actualPresence")
  }
  
  return(output)
}

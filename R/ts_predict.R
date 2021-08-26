#' Predictive Time Series
#'
#' Predict occupancy, detection, and presence
#' probabilities with uncertainty over time.
#'
#' @param model \code{btvoccu} output
#' @param x numeric array: predictors array
#' @param sites character vector
#' @param season integer
#' @param periods integer vector
#' @param credible length \code{3} vector of quantiles
#' @param value 'Occupancy', 'Detection', or 'Presence'
#' @param spline logical: default \code{TRUE} to fit smoothing splines
#' @param spline.predict logical: default \code{FALSE} to not predict missing values using splines
#' @param nknots see \code{smooth.spline()}
#' @param ndraws integer: # of posterior draws
#' @param burnin percent of posterior samples to burn
#' @param thin integer: default \code{1} implies no thinning, \code{2} means keep every other draw, etc.
#'
#' @return data frame
#'
#' @export
ts_predict <- function(model,
                       x,
                       sites,
                       season,
                       periods,
                       credible = c(.025,.5,.975),
                       value = 'Occupancy',
                       spline = TRUE,
                       spline.predict = FALSE,
                       nknots = 5,
                       ndraws = 1000,
                       burnin = .5,
                       thin = 1
                       ){

  # sort because posterior_sigmoids sorts
  sites <- sort(sites)
  periods <- sort(periods)

  # access posterior sigmoid probabilities

  sigmoids <- posterior_sigmoids(model,
                                 x,
                                 sites,
                                 season,
                                 periods,
                                 ndraws,
                                 burnin,
                                 thin
                                 )
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

        spl <- smooth.spline(sub.nna$Period, sub.nna$Lower,
                             nknots = 5)

        if(spline.predict){
          df$Lower[indices] <- predict(spl, sub$Period)$y
        } else{
          df$Lower[indices.nna] <- spl$y
        }

        # Middle

        spl <- smooth.spline(sub.nna$Period, sub.nna$Middle,
                             nknots = 5)

        if(spline.predict){
          df$Middle[indices] <- predict(spl, sub$Period)$y
        } else{
          df$Middle[indices.nna] <- spl$y
        }

        # Upper

        spl <- smooth.spline(sub.nna$Period, sub.nna$Upper,
                             nknots = 5)

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

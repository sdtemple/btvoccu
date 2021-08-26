#' Predictive Time Series Plots
#'
#' Visualize occupancy, detection, and presence
#' probabilities with uncertainty over time.
#'
#' @param model \code{btvoccu} output
#' @param x numeric array: predictors array
#' @param sites vector
#' @param season length \code{1} vector
#' @param periods vector
#' @param credible length \code{3} vector of quantiles
#' @param value 'Occupancy', 'Detection', or 'Presence'
#' @param ribbons logical: default \code{TRUE} to show uncertainty ribbons
#' @param spline logical: default \code{FALSE} to not fit smoothing splines
#' @param spline.predict logical: default \code{FALSE} to not predict missing values using splines
#' @param nknots see \code{smooth.spline()}
#' @param xaxis character: x-axis label
#' @param legendtile character: legend label
#' @param dots logical: default \code{FALSE} to not show actual datapoints
#' @param ndraws integer: # of posterior draws
#' @param burnin percent of posterior samples to burn
#' @param thin integer: default \code{1} implies no thinning, \code{2} means keep every other draw, etc.
#'
#' @return \code{ggplot} object
#'
#' @export
ts_predict_plot <- function(model,
                            x,
                            sites,
                            season,
                            periods,
                            credible = c(.025,.5,.975),
                            value = 'Occupancy',
                            ribbons = TRUE,
                            spline = FALSE,
                            spline.predict = FALSE,
                            nknots = 5,
                            xaxis = 'Epiweek',
                            legendtitle = 'Locality',
                            dots = FALSE,
                            ndraws = 1000,
                            burnin = .5,
                            thin = 1
                            ){

  df <- ts_predict(model,
                   x,
                   sites,
                   season,
                   periods,
                   credible,
                   value,
                   spline,
                   spline.predict,
                   nknots,
                   ndraws,
                   burnin,
                   thin
                   )


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

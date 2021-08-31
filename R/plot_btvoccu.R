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
#' @param value "Occupancy", "Detection", or "Presence"
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
#' @examples 
#' plot_btvoccu(mdl, x, c("TOR", "YRK", "PEE"), 2005, 20:40)
#' plot_btvoccu(mdl, x, c("TOR", "YRK", "PEE"), 2005, 20:40, spline = T)
#' plot_btvoccu(mdl, x, c("TOR", "YRK", "PEE", "THB", "ALG", "OTT"), 2005, 20:40, spline = T)
#' plot_btvoccu(mdl, x, c("TOR", "YRK", "PEE"), 2005, 20:40, "Detection")
#' plot_btvoccu(mdl, x, c("TOR", "YRK", "PEE"), 2005, 20:40, "Presence")
#'
#' @export
plot_btvoccu <- function(model,
                         x,
                         sites,
                         season,
                         periods,
                         value = "Occupancy",
                         burnin = .5,
                         thin = 1,
                         credible = c(.025,.5,.975),
                         ndraws = 1000,
                         ribbons = TRUE,
                         dots = FALSE,
                         spline = FALSE,
                         spline.predict = FALSE,
                         nknots = 5,
                         xaxis = "Period",
                         legendtitle = "Site"){
  
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
          ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper), alpha = .1, linetype = 2, show.legend = F) +
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
          ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper), alpha = .1, linetype = 2, show.legend = F) +
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
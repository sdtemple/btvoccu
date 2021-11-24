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
#' @param font.size see \code{theme_classic()}
#' @param legend.key.width see \code{theme()}
#' @param legend.key.height see \code{theme()}
#'
#' @return \code{ggplot} object
#' 
#' @examples 
#' plot_covariate(x, c("TOR", "YRK", "PEE"), 2005, 20:40, "popdensity", "Population Density")
#' plot_covariate(x, c("TOR", "YRK", "PEE"), 2005, 20:40, "temperature", "Temperature (Celsius)")
#'
#' @export
plot_covariate <- function(x,
                           sites,
                           season,
                           periods,
                           covariate,
                           yaxis,
                           xaxis = "Period",
                           legendtitle = "Site",
                           font.size = 8,
                           legend.key.width = 2,
                           legend.key.height = 0){
  
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
  colnames(df) <- c("Site","Season","Period","Covariate")
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
      ggplot2::theme_classic(base_size = font.size) +
      ggplot2::theme(legend.key.height = ggplot2::unit(legend.key.height, "mm"), 
                     legend.key.width = ggplot2::unit(legend.key.width, "mm")) + 
      ggplot2::labs(x = xaxis, y = yaxis, colour = legendtitle) +
      ggplot2::scale_y_continuous(limits = c(mn - .001, mx + .001))
  )
}
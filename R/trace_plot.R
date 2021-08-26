#' Trace Plots
#'
#' Plot trace of posterior samples.
#'
#' @param model \code{btvoccu} output
#' @param effects "betas", "alphas", or "thetas"
#' @param mnridx integer index
#' @param burnin percent of posterior samples to burn
#'
#' @export
trace_plot <- function(model,
                       effects,
                       mnridx,
                       burnin = .5
                       ){

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
  for(n in 2:model$nchains){
    lines(model[[mjridx]][n, burn, mnridx],
          col = colors[n])
  }
}

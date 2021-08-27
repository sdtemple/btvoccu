#' Trace plots
#'
#' Plot trace of posterior samples.
#'
#' @param model \code{btvoccu} model object
#' @param effects "betas", "alphas", or "thetas"
#' @param mnridx integer index
#' @param burnin percent of posterior samples to burn
#' 
#' @examples 
#' trace_plot(mdl, "betas", 1)
#' trace_plot(mdl, "alphas", 1, 0.10)
#' for(i in 1:length(mdl$occueffs)){plot_trace(mdl, "betaS", i, burnin = 0)}
#'
#' @export
plot_trace <- function(model, effects, mnridx, burnin = 0){
  
  mjridx <- which(names(model) == effects, arr.ind = T)
  colors <- grDevices::rainbow(model$nchains)
  burn <- -(1:(floor(burnin * model$niter)))
  
  # ylim
  mx <- max(model[[mjridx]][1, burn, mnridx])
  mn <- min(model[[mjridx]][1, burn, mnridx])
  for(n in 2:model$nchains){
    mx <- max(mx, model[[mjridx]][n, burn, mnridx])
    mn <- min(mn, model[[mjridx]][n, burn, mnridx])
  }
  
  # plot
  if(effects == "betas"){
    name <- "occueffs"
  }
  else if(effects == "alphas"){
    name <- "deteffs"
  }
  else {
    name <- "speffs"
  }
  plot(model[[mjridx]][1, burn, mnridx],
       ylab = getElement(model, name)[mnridx],
       xlab = "index",
       type = "l",
       col = colors[1],
       ylim = c(mn, mx),
       main = effects)
  for(n in 2:model$nchains){lines(model[[mjridx]][n, burn, mnridx], col = colors[n])}
}

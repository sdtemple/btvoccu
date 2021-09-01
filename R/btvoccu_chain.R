#' Chaining 'btvoccu' models
#' 
#' Concatenate independent MCMC chains for \code{btvoccu} model objects
#' 
#' @param model1 \code{btvoccu} model object
#' @param model2 \code{btvoccu} model object
#' 
#' @return \code{btvoccu} model object
#' 
#' @export
btvoccu_chain <- function(model1, model2){
  o1 <- model1$occueffs; o2 <- model2$occueffs
  d1 <- model1$deteffs; d2 <- model2$deteffs
  s1 <- model1$speffs; s2 <- model2$speffs
  if(identical(model1$sites, model2$sites)){
    if(identical(model1$seasons, model2$seasons)){
      if(identical(model1$periods, model2$periods)){
        if(identical(model1$niter, model2$niter)){
          if(identical(model1$link, model2$link)){
            if(identical(o1, o2)){
              if(identical(d1, d2)){
                if(identical(s1, s2)){
                  output <- list()
                  output[[1]] <- abind::abind(model1$betas, model2$betas, along = 1)
                  output[[2]] <- abind::abind(model1$alphas, model2$alphas, along = 1)
                  if(is.null(s1)){
                    output[[3]] <- NULL
                  } else{
                    output[[3]] <- abind::abind(model1$thetas, model2$thetas, along = 1)
                  }
                  output[[4]] <- abind::abind(model1$Sigmabetas, model2$Sigmabetas, along = 1)
                  output[[5]] <- abind::abind(model1$Sigmaalphas, model2$Sigmaalphas, along = 1)
                  if(is.null(s1)){
                    output[[6]] <- NULL
                  } else{
                    output[[6]] <- abind::abind(model1$sigmathetas, model2$sigmathetas, along = 1)
                  }
                  output[[7]] <- o1
                  output[[8]] <- d1
                  output[[9]] <- s1
                  output[[10]] <- model1$sites
                  output[[11]] <- model1$seasons
                  output[[12]] <- model1$periods
                  output[[13]] <- model1$niter
                  output[[14]] <- model1$nchains + model2$nchains
                  output[[15]] <- model1$link
                  names(output) <- names(model1)
                  return(output)
                }
              }
            }
          }
        }
      }
    }
  }
  print("Chains differ: these are not the same models.")
  return()
}

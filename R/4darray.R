# Manage covariate and response data in 4-dimensional array
# August 25, 2021

# Append to 4d array ------------------------------------------------------

#' Append transformed variable
#'
#' Append a transformed variable to a 4D array.
#'
#' @param x numeric array
#' @param covariate character name
#' @param fn function
#' @param prefix character prefix for column name
#'
#' @return numeric array
#' 
#' @examples 
#' x <- append_transformed_term(x, "precipitation", sqrt, "sqrt.")
#' x <- append_transformed_term(x, "popdensity", log10, "log10.")
#'
#' @export
append_transformed_term <- function(x, covariate, fn, prefix){
  covind <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  x <- abind::abind(x, fn(x[,,,covind]))
  dimnames(x)[[4]][dim(x)[4]] <- paste(prefix, covariate, sep = '')
  return(x)
}

#' Append interaction term
#'
#' Append an interaction variable to a 4D array.
#'
#' @param x numeric array
#' @param covariate character name
#' @param predictor character name
#'
#' @return numeric array
#' 
#' @examples 
#' x <- append_interaction_term(x, "temperature", "precipitation")
#'
#' @export
append_interaction_term <- function(x, covariate, predictor){
  covind <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  predind <- which(dimnames(x)[[4]] == predictor, arr.ind = T)
  x <- abind::abind(x, x[,,,covind] * x[,,,predind])
  dimnames(x)[[4]][dim(x)[4]] <- paste(covariate, '*', predictor, sep = '')
  return(x)
}

#' Append lagged covariates
#'
#' Append lagged variables to a 4D array.
#'
#' @param x numeric array
#' @param covariate character name
#' @param lags integers
#'
#' @return numeric array
#' 
#' @examples
#' x <- append_lagged_terms(x, "precipitation", 2)
#' x <- append_lagged_terms(x, "precipitation", c(2, 4, 6))
#'
#' @export
append_lagged_terms <- function(x, covariate, lags){
  
  covind <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  for(lag in lags){
    xx <- array(NA, dim = c(dim(x)[1:3], 1))
    for(i in 1:dim(x)[1]){
      for(j in 1:dim(x)[2]){
        for(k in (lag + 1):(dim(x)[3])){
          xx[i,j,k,1] <- x[i,j,(k-lag),covind]
        }
      }
    }
    x <- abind::abind(x, xx)
    dimnames(x)[[4]][length(dimnames(x)[[4]])] <- paste('lag', lag, covariate, sep ='')
  }
  
  return(x)
}

#' Append leading covariates
#'
#' Append leading variables to a 4D array.
#'
#' @param x numeric array
#' @param covariate character name
#' @param leads integers
#'
#' @return numeric array
#'
#' @export
append_leading_terms <- function(x, covariate, leads){
  
  covind <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  for(lead in leads){
    xx <- array(NA, dim = c(dim(x)[1:3], 1))
    for(i in 1:dim(x)[1]){
      for(j in 1:dim(x)[2]){
        for(k in 1:(dim(x)[3] - lead)){
          xx[i,j,k,1] <- x[i,j,(k + lead),covind]
        }
      }
    }
    x <- abind::abind(x, xx)
    dimnames(x)[[4]][length(dimnames(x)[[4]])] <- paste('lead', lead, covariate, sep ='')
  }
  
  return(x)
}

#' Append principal components
#'
#' Append principal components to a 4D array.
#' Plots from \code{factoextra} describe the PCA.
#'
#' @param x numeric array
#' @param covariates vector of covariates
#' @param npc integer: number of principal components
#' @param prefix character prefix for column name
#'
#' @return numeric array
#' 
#' @examples 
#' x <- append_principal_components(x, dimnames(x)[[4]][12:23], 2, "climate")
#'
#' @export
append_principal_components <- function(x, covariates, npc, prefix){
  
  z <- apply(x, 4, as.numeric)
  indices <- which(dimnames(z)[[2]] %in% covariates, arr.ind = T)
  components <- prcomp(~., data=data.frame(z[,indices]), na.action = na.exclude)
  print(factoextra::fviz_eig(components))
  print(factoextra::fviz_pca_var(components, repel = T))
  z.pca <- factoextra::get_pca_ind(components)
  x.pca <- array(dim = c(dim(x)[1:3], npc),
                 dimnames = list(dimnames(x)[[1]], dimnames(x)[[2]], dimnames(x)[[3]],
                                 paste(prefix, '.pc', 1:npc, sep = '')))
  
  for(l in 1:npc){
    n <- 1
    for(k in 1:dim(x)[3]){
      for(j in 1:dim(x)[2]){
        for(i in 1:dim(x)[1]){
          x.pca[i,j,k,l] <- z.pca$coord[n,l]
          n <- n + 1
        }
      }
    }
  }
  
  return(abind::abind(x, x.pca))
}

#' Append spatial random effects
#'
#' Append eigenvectors from Moran operator decomposition
#' to a 4D array.
#'
#' @param x numeric array
#' @param arealeffs character vector: areal effects
#' @param A adjacency matrix
#'
#' @return numeric array
#' 
#' @examples 
#' x <- append_spatial_effects(x, c("urban", "agri", "forest"), adjacency_matrix)
#'
#' @export
append_spatial_effects <- function(x, arealeffs, A){
  
  # local function
  spmatrix <- function(X, A){
    P <- X %*% tcrossprod(solve(crossprod(X, X)), X)
    I <- diag(nrow(P))
    M <- (I - P) %*% A %*% (I - P)
    return(eigen(M))
  }
  
  xx <- subset_4darray(x, 4, arealeffs)
  eM <- spmatrix(xx[,1,1,], A)
  print(eM$values)
  
  sparray <- array(NA,
                   dim = c(dim(x)[1:3],
                           ncol(eM$vectors)),
                   dimnames = list(dimnames(x)[[1]],
                                   dimnames(x)[[2]],
                                   dimnames(x)[[3]],
                                   paste('sp',
                                         1:ncol(eM$vectors),
                                         sep = '')
                   )
  )
  
  for(j in 1:dim(x)[2]){
    for(k in 1:dim(x)[3]){
      sparray[,j,k,] <- eM$vectors
    }
  }
  
  return(abind::abind(x, sparray))
}

#' Q matrix for ICAR component
#'
#' Compute the Q matrix for an ICAR component.
#'
#' @param A adjacency matrix
#'
#' @return matrix
icarQ <- function(A){
  D <- diag(apply(A, 2, sum))
  return(D - A)
}


# Split or subset 4d array -------------------------------------------------

#' Subset 4D array
#'
#' Subset 4-dimensional arrays along any major dimension.
#'
#' @param x numeric array
#' @param d integer: array dimension
#' @param labels character vector: predictor names
#'
#' @return numeric array
#' 
#' @examples 
#' w <- subset_4darray(x, 1, c("TOR", "YRK", "PEE"))
#' w <- subset_4darray(x, 2, 2005:2010)
#' w <- subset_4darray(x, 3, 20:40)
#' w <- subset_4darray(x, 4, c("precipitation", "temperature", "popdensity"))
#'
#' @export
subset_4darray <- function(x, d, labels){
  indices <- which(dimnames(x)[[d]] %in% labels, arr.ind = T)
  if(d == 1){
    return(x[indices,,,, drop = F])
  } else if(d == 2){
    return(x[,indices,,, drop = F])
  } else if(d == 3){
    return(x[,,indices,, drop = F])
  } else if(d == 4){
    return(x[,,,indices, drop = F])
  }
}

#' Train-validate-test for 4D array
#'
#' Split 4D array into training,
#' validation, and testing arrays.
#' Splitting is done by setting
#' values to NA because main function
#' \code{btvoccu()} ignores NA responses
#' when fitting.
#'
#' @param y binary array: response array
#' @param p proportions (vector length \code{1} or \code{2})
#' @param sites vector
#' @param seasons vector
#' @param periods vector
#'
#' @return list
#' 
#' @examples
#' ytest <- subset_4darray(y, 2, c(2008, 2012, 2016))
#' ysplit <- subset_4darray(y, 2, c(2002:2007, 2009:2011, 2013:2015, 2017))
#' ysplit <- split_4darray(y, .75, periods = 20:40)
#' ytrain <- ysplit$train
#' yvalid <- ysplit$test
#' 
#' ysplit <- split_4darray(y, c(6, .2), periods = 20:40)
#'
#' @export
split_4darray <- function(y, p, sites = NULL, seasons = NULL, periods = NULL){
  
  `%notin%` <- Negate(`%in%`)
  
  # when NULL
  if(is.null(sites)){
    sites <- dimnames(y)[[1]]
  }
  if(is.null(seasons)){
    seasons <- dimnames(y)[[2]]
  }
  if(is.null(periods)){
    periods <- dimnames(y)[[3]]
  }
  
  # some error handling
  if(sum(p) > 1){
    stop('Invalid splitting proportions')
  }
  if(sum(p < 0) > 0){
    stop('Invalid splitting proportions')
  }
  if(length(p) > 2){
    stop('p must be of length 1 or 2')
  }
  
  # concatenate indices
  splitter <- c()
  for(i in as.character(sites)){
    for(j in as.character(seasons)){
      for(k in as.character(periods)){
        if(!is.na(y[i,j,k,1])){
          splitter <- c(paste(i, '.', j, '.', k, sep = ''), splitter)
        }
      }
    }
  }
  len <- length(splitter)
  
  if(length(p) == 1){ # train-test split
    
    # sample indices
    train_indices <- sample(splitter, floor(p * len))
    test_indices <- splitter[which(splitter %notin% train_indices, arr.ind = T)]
    
    # split based on periods
    train_arr <- array(NA, dim = c(length(train_indices), 3))
    test_arr <- array(NA, dim = c(length(test_indices), 3))
    
    for(n in 1:length(train_indices)){
      str <- train_indices[n]
      locs <- gregexpr('.', str, fixed = T)[[1]][1:2]
      train_arr[n,] <- c(substr(str, 1, locs[1] - 1),
                         substr(str, locs[1] + 1, locs[2] - 1),
                         substr(str, locs[2] + 1, nchar(str))
      )
    }
    
    for(n in 1:length(test_indices)){
      str <- test_indices[n]
      locs <- gregexpr('.', str, fixed = T)[[1]][1:2]
      test_arr[n,] <- c(substr(str, 1, locs[1] - 1),
                        substr(str, locs[1] + 1, locs[2] - 1),
                        substr(str, locs[2] + 1, nchar(str))
      )
    }
    
    # set values to NA
    train_array <- y
    test_array <- y
    
    # set non-training values to NA
    for(n in 1:dim(test_arr)[1]){
      i <- test_arr[n,1]
      j <- test_arr[n,2]
      k <- test_arr[n,3]
      train_array[i,j,k,1] <- NA
    }
    
    # set non-testing values to NA
    for(n in 1:dim(train_arr)[1]){
      i <- train_arr[n,1]
      j <- train_arr[n,2]
      k <- train_arr[n,3]
      test_array[i,j,k,1] <- NA
    }
    
    # formatting
    output <- list()
    output[[1]] <- train_array
    output[[2]] <- test_array
    names(output) <- c('train','test')
    
  } else if(length(p) == 2){ # train-test-validate split
    
    # sample indices
    prop <- 1 - sum(p)
    test_indices <- sample(splitter, floor(prop * len))
    splitter <- splitter[which(splitter %notin% test_indices, arr.ind = T)]
    
    prop <- p[1] / sum(p)
    len <- length(splitter)
    train_indices <- sample(splitter, floor(prop * len))
    validate_indices <- splitter[which(splitter %notin% train_indices, arr.ind = T)]
    
    # split based on periods
    train_arr <- array(NA, dim = c(length(train_indices), 3))
    validate_arr <- array(NA, dim = c(length(validate_indices), 3))
    test_arr <- array(NA, dim = c(length(test_indices), 3))
    
    for(n in 1:length(train_indices)){
      str <- train_indices[n]
      locs <- gregexpr('.', str, fixed = T)[[1]][1:2]
      train_arr[n,] <- c(substr(str, 1, locs[1] - 1),
                         substr(str, locs[1] + 1, locs[2] - 1),
                         substr(str, locs[2] + 1, nchar(str))
      )
    }
    
    for(n in 1:length(validate_indices)){
      str <- validate_indices[n]
      locs <- gregexpr('.', str, fixed = T)[[1]][1:2]
      validate_arr[n,] <- c(substr(str, 1, locs[1] - 1),
                            substr(str, locs[1] + 1, locs[2] - 1),
                            substr(str, locs[2] + 1, nchar(str))
      )
    }
    
    for(n in 1:length(test_indices)){
      str <- test_indices[n]
      locs <- gregexpr('.', str, fixed = T)[[1]][1:2]
      test_arr[n,] <- c(substr(str, 1, locs[1] - 1),
                        substr(str, locs[1] + 1, locs[2] - 1),
                        substr(str, locs[2] + 1, nchar(str))
      )
    }
    
    # set values to NA
    train_array <- y
    validate_array <- y
    test_array <- y
    
    # set non-training values to NA
    for(n in 1:dim(test_arr)[1]){
      i <- test_arr[n,1]
      j <- test_arr[n,2]
      k <- test_arr[n,3]
      train_array[i,j,k,1] <- NA
    }
    
    for(n in 1:dim(validate_arr)[1]){
      i <- validate_arr[n,1]
      j <- validate_arr[n,2]
      k <- validate_arr[n,3]
      train_array[i,j,k,1] <- NA
    }
    
    # set non-validating values to NA
    for(n in 1:dim(test_arr)[1]){
      i <- test_arr[n,1]
      j <- test_arr[n,2]
      k <- test_arr[n,3]
      validate_array[i,j,k,1] <- NA
    }
    
    for(n in 1:dim(train_arr)[1]){
      i <- train_arr[n,1]
      j <- train_arr[n,2]
      k <- train_arr[n,3]
      validate_array[i,j,k,1] <- NA
    }
    
    # set non-testing values to NA
    for(n in 1:dim(validate_arr)[1]){
      i <- validate_arr[n,1]
      j <- validate_arr[n,2]
      k <- validate_arr[n,3]
      test_array[i,j,k,1] <- NA
    }
    
    for(n in 1:dim(train_arr)[1]){
      i <- train_arr[n,1]
      j <- train_arr[n,2]
      k <- train_arr[n,3]
      test_array[i,j,k,1] <- NA
    }
    
    # formatting
    output <- list()
    output[[1]] <- train_array
    output[[2]] <- test_array
    output[[3]] <- validate_array
    names(output) <- c('train','test','validate')
  }
  
  return(output)
}


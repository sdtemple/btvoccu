#' Train-Test Split for 4D Array
#'
#' Split 4D array into train,
#' test, and validate arrays.
#' Splitting is done by setting
#' values to NA because main function
#' \code{btvoccu()} ignores NA responses
#' when fitting.
#'
#' @param y binary array: response array
#' @param props length \code{1} or \code{2} proportions
#' @param sites vector
#' @param seasons vector
#' @param periods vector
#'
#' @return list
#'
#' @export
split_4darray <- function(y, 
                          props,
                          sites = NULL,
                          seasons = NULL, 
                          periods = NULL){
  
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
  
  if(sum(props) > 1){
    stop('Invalid splitting proportions')
  }
  
  if(sum(props < 0) > 0){
    stop('Invalid splitting proportions')
  }
  
  if(length(props) > 2){
    stop('props must be of length 1 or 2')
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
  
  
  `%notin%` <- Negate(`%in%`)
  if(length(props) == 1){
    
    # train-test split
    
    
    # sample indices
    
    train_indices <- sample(splitter, floor(props * len))
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
    
  } else if(length(props) == 2){
    
    # train-test-validate split
    
    # sample indices
    
    prop <- 1 - sum(props)
    test_indices <- sample(splitter, floor(prop * len))
    splitter <- splitter[which(splitter %notin% test_indices, arr.ind = T)]
    
    prop <- props[1] / sum(props)
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

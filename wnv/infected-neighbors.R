# Infected Neighbors
# Seth Temple, sdtemple@lanl.gov
# September 11, 2021

library(btvoccu)

y <- readRDS("data/response.rds")
y <- apply(y, 1:4, as.numeric)
A <- readRDS("data/adjacency.rds")

# Neighbor Covariate ------------------------------------------------------

w <- array(NA, dim = c(34, 16, 52, 2))
for(i in 1:34){
  adj <- c(which(A[i,] == 1, arr.ind = T), i)
  for(j in 1:16){
    for(k in 18:44){
      slice <- y[adj,j,k,]
      nna <- sum(!is.na(slice))
      sm <- sum(slice, na.rm = T)
      w[i,j,k,1] <- sm
      w[i,j,k,2] <- sm / (nna + 1)
    }
  }
}
for(i in 1:34){
  for(j in 1:16){
    for(k in c(1:17, 45:52)){
      w[i,j,k,1:2] <- 0
    }
  }
}

dimnames(w)[[1]] <- dimnames(y)[[1]]
dimnames(w)[[2]] <- dimnames(y)[[2]]
dimnames(w)[[3]] <- dimnames(y)[[3]]
dimnames(w)[[4]] <- c("infected.neighbors", "infected.neighbor.rate")
x <- readRDS("data/fewer-covariates.rds")
x <- apply(x, 1:4, as.numeric)
x <- abind::abind(x, w)
x <- append_lagged_terms(x, "infected.neighbors", 1:2)
x <- append_lagged_terms(x, "infected.neighbor.rate", 1:2)
saveRDS(x, "data/final-covariates.rds")

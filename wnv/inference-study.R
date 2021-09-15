# btvoccu Inference Study
# Seth Temple, sdtemple@lanl.gov
# September 9, 2021

library(btvoccu)
set.seed(992021)

# Common Inputs -----------------------------------------------------------

# study design
nsites <- 10
nseasons <- 5
nperiods <- 20

# MCMC parameters
niter <- 1000
nchains <- 2
print.interval <- 100

# effects
betas <- c(2, 1, 1)
alphas <- c(-1/2, 0, 1)

# site-specific means
occusite <- 0
detsite <- -1

# time-varying covariates
occutime <- array(NA, dim = c(1, 20))
occutime[1,] <- sin(2 * pi * 1:20 / 40 + pi)
dettime <- array(NA, dim = c(1, 20))
dettime[1,] <- sin(2 * pi * 1:20 / 40)

# Simulation Study --------------------------------------------------------

nrep <- 20
storage <- array(NA, dim = c(length(betas) + length(alphas),
                             nrep, 3))
rp <- rep(NA, nrep)
for(i in 1:nrep){
  print(paste("simulation", i))
  # simulate data
  d <- btvoccu_simstudy(nsites,
                        nseasons,
                        nperiods,
                        betas,
                        alphas,
                        occusite,
                        detsite,
                        occutime,
                        dettime,
                        noise = 0.25)
  y <- d$y
  sites <- d$sites
  seasons <- d$seasons
  periods <- d$periods
  bvar <- d$occueffs
  avar <- d$deteffs
  XW <- abind::abind(d$X, d$W)
  
  # modeling
  m <- btvoccu(niter, XW, y, sites, seasons, periods, bvar, avar, 
               nchains = nchains, print.interval = print.interval)
  smry <- posterior_effects(m)
  smry <- apply(smry[,2:6], 1:2, as.numeric)
  storage[,i,1] <- smry[,2]
  storage[,i,2] <- smry[,3] - smry[,1]
  storage[,i,3] <- smry[,5]
  rp[i] <- posterior_check(m, XW, y)$relativePresence
}
saveRDS(storage, "simstudy/table1.rds")
saveRDS(rp, "simstudy/vector1.rds")

# with NAs (missingness)
nrep <- 20
storage <- array(NA, dim = c(length(betas) + length(alphas),
                             nrep, 3))
rp <- rep(NA, nrep)
for(i in 1:nrep){
  print(paste("simulation", i))
  # simulate data
  d <- btvoccu_simstudy(nsites,
                        nseasons,
                        nperiods,
                        betas,
                        alphas,
                        occusite,
                        detsite,
                        occutime,
                        dettime,
                        noise = 0.25,
                        naprob = 1/2) # only change here
  y <- d$y
  sites <- d$sites
  seasons <- d$seasons
  periods <- d$periods
  bvar <- d$occueffs
  avar <- d$deteffs
  XW <- abind::abind(d$X, d$W)
  
  # modeling
  m <- btvoccu(niter, XW, y, sites, seasons, periods, bvar, avar, 
               nchains = nchains, print.interval = print.interval)
  smry <- posterior_effects(m)
  smry <- apply(smry[,2:6], 1:2, as.numeric)
  storage[,i,1] <- smry[,2]
  storage[,i,2] <- smry[,3] - smry[,1]
  storage[,i,3] <- smry[,5]
  rp[i] <- posterior_check(m, XW, y)$relativePresence
}
saveRDS(storage, "simstudy/table2.rds")
saveRDS(rp, "simstudy/vector2.rds")

# Results -----------------------------------------------------------------

# inference table
t1 <- readRDS("simstudy/table1.rds")
t2 <- readRDS("simstudy/table2.rds")
t0 <- matrix(c(betas, alphas), ncol = 1)
t1 <- apply(t1, c(1,3), mean)
t2 <- apply(t2, c(1,3), mean)
out <- cbind(t0, t1, t2)
library(xtable)
xtable(out, digits = 3)

# predictive performance
v1 <- readRDS("simstudy/vector1.rds")
v2 <- readRDS("simstudy/vector2.rds")
mean(v1); sd(v1)
mean(v2); sd(v2)


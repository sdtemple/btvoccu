# btvoccu Inference Study
# Seth Temple, sdtemple@lanl.gov
# September 29, 2021

library(btvoccu)
set.seed(9292021)

# Common Inputs -----------------------------------------------------------

# study design
nsites <- 35
nseasons <- 3
nperiods <- 20

# MCMC parameters
niter <- 5000
nchains <- 1
print.interval <- 500
nrep <- 40

# site-specific means
occusite <- 1/4
detsite <- -1/2

# time-varying covariates
occutime <- array(NA, dim = c(1, 20))
occutime[1,] <- sin(2 * pi * 1:20 / 40) / 2
dettime <- array(NA, dim = c(1, 20))
dettime[1,] <- sin(2 * pi * 1:20 / 40 - pi / 4) / 2

# Simulation Study 1 --------------------------------------------------------

# effects
betas <- c(3/4, 1, 1)
alphas <- c(1/2, -1, 1/2)

storage <- array(NA, dim = c(length(betas) + length(alphas),
                             nrep, 2))
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
                        dettime)
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
  rp[i] <- posterior_check(m, XW, y)$relativePresence
}
saveRDS(storage, "simstudy/table1.rds")
saveRDS(rp, "simstudy/vector1.rds")

# Simulation Study 2 --------------------------------------------------------

# effects
betas <- c(3/4, 1, 1)
alphas <- c(-1, -1, 1/2)

storage <- array(NA, dim = c(length(betas) + length(alphas),
                             nrep, 2))
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
                        dettime)
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
  rp[i] <- posterior_check(m, XW, y)$relativePresence
}
saveRDS(storage, "simstudy/table2.rds")
saveRDS(rp, "simstudy/vector2.rds")

# Simulation Study 3 --------------------------------------------------------

# effects
betas <- c(-1/2, 1, 1)
alphas <- c(1/2, -1, 1/2)

storage <- array(NA, dim = c(length(betas) + length(alphas),
                             nrep, 2))
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
                        dettime)
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
  rp[i] <- posterior_check(m, XW, y)$relativePresence
}
saveRDS(storage, "simstudy/table3.rds")
saveRDS(rp, "simstudy/vector3.rds")

# Simulation Study 4 --------------------------------------------------------

# effects
betas <- c(-1/2, 1, 1)
alphas <- c(-1, -1, 1/2)

storage <- array(NA, dim = c(length(betas) + length(alphas),
                             nrep, 2))
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
                        dettime)
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
  rp[i] <- posterior_check(m, XW, y)$relativePresence
}
saveRDS(storage, "simstudy/table4.rds")
saveRDS(rp, "simstudy/vector4.rds")

# Results -----------------------------------------------------------------

# inference table
library(xtable)
t1 <- readRDS("simstudy/table1.rds")
t2 <- readRDS("simstudy/table2.rds")
t3 <- readRDS("simstudy/table3.rds")
t4 <- readRDS("simstudy/table4.rds")
t1 <- apply(t1, c(1,3), mean)
t2 <- apply(t2, c(1,3), mean)
t3 <- apply(t3, c(1,3), mean)
t4 <- apply(t4, c(1,3), mean)
out <- rbind(t1, t2, t3, t4)
t0 <- c(c(3/4,1,1),
        c(1/2,-1,1/2),
        c(3/4,1,1),
        c(-1,-1,1/2),
        c(-1/2,1,1),
        c(1/2,-1,1/2),
        c(-1/2,1,1),
        c(-1,-1,1/2))
out <- cbind(t0, out)
xtable(out[,1:3], digits = 3)

# predictive performance
v1 <- readRDS("simstudy/vector1.rds")
v2 <- readRDS("simstudy/vector2.rds")
v3 <- readRDS("simstudy/vector3.rds")
v4 <- readRDS("simstudy/vector4.rds")
mean(v1); sd(v1)
mean(v2); sd(v2)
mean(v3); sd(v3)
mean(v4); sd(v4)


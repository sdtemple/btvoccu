# btvoccu Time Study
# Seth Temple, sdtemple@lanl.gov
# September 11, 2021

library(btvoccu)
library(rstan)

# Common Inputs -----------------------------------------------------------

# study design
nsites <- 10
nseasons <- 5
nperiods <- 20

# MCMC parameters
niter <- 2000
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

d <- btvoccu_simstudy(nsites,
                      nseasons,
                      nperiods,
                      betas,
                      alphas,
                      occusite,
                      detsite,
                      occutime,
                      dettime,
                      naprob = 1/2)
y <- d$y
sites <- d$sites
seasons <- d$seasons
periods <- d$periods
bvar <- d$occueffs
avar <- d$deteffs
XW <- abind::abind(d$X, d$W)

# btvoccu -----------------------------------------------------------------

m <- btvoccu(niter, XW, y, sites, seasons, periods, bvar, avar, 
             nchains = nchains, print.interval = print.interval)


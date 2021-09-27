# btvoccu Time Study
# Seth Temple, sdtemple@lanl.gov
# September 27, 2021

library(btvoccu)
library(rstan)
set.seed(9272021)

# Setup -------------------------------------------------------------------

y <- readRDS("data/response.rds")
x <- readRDS("data/final-covariates.rds")
y <- apply(y, 1:4, as.numeric)
x <- apply(x, 1:4, as.numeric)

sites <- dimnames(x)[[1]]
seasons <- dimnames(y)[[2]]
periods <- 18:44

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")

niter <- 2000
nchains <- 2
print.interval <- niter / 10


# btvoccu -----------------------------------------------------------------

ptm <- proc.time()
m <- btvoccu(niter, x, y, sites, seasons, periods, 
             occueffs, deteffs,
             nchains = nchains, print.interval = print.interval)
proc.time() - ptm


# rstan -------------------------------------------------------------------

X <- subset_4darray(x, 4, occueffs)
X <- subset_4darray(X, 1, sites)
X <- subset_4darray(X, 2, seasons)
X <- subset_4darray(X, 3, periods)
W <- subset_4darray(x, 4, deteffs)
W <- subset_4darray(W, 1, sites)
W <- subset_4darray(W, 2, seasons)
W <- subset_4darray(W, 3, periods)
Y <- subset_4darray(y, 1, sites)
Y <- subset_4darray(Y, 2, seasons)
Y <- subset_4darray(Y, 3, periods)
na <- !is.na(Y)
na <- apply(na, 1:4, as.numeric)
yna <- replace(Y, is.na(Y), 0)
wna <- replace(W, is.na(Y), 0)

dd <- list(nsites = length(sites),
           nseasons = length(seasons),
           nperiods = length(periods),
           nbeta = length(occueffs),
           nalpha = length(deteffs),
           yna = yna,
           na = na,
           X = X,
           W = wna)

stanfile <- "nuts_logit_btvoccu.stan"
s <- stan(stanfile,
          data = dd,
          chains = nchains,
          warmup = niter / 2,
          iter = niter,
          core = 1)


# checking ----------------------------------------------------------------

# btvoccu
posterior_effects(m)
par(mfrow=c(2,3))
for(i in 1:5){plot_trace(m, "betas", i, burnin = 1/2)}
par(mfrow=c(2,2))
for(i in 1:4){plot_trace(m, "alphas", i, burnin = 1/2)}

# rstan
par(mfrow=c(1,1))
plot(s, pars = c("beta", "alpha"))
print(s, pars = c('beta','alpha'), probs = c(.05,.5,.95))
traceplot(s, "beta")
traceplot(s, "alpha")

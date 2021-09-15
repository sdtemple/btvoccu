# btvoccu Models A*
# Seth Temple, sdtemple@lanl.gov
# August 26, 2021

library(btvoccu)
set.seed(8262021)

# Setup -------------------------------------------------------------------

y <- readRDS("data/response-train.rds")
x <- readRDS("data/final-covariates.rds")
A <- readRDS("data/adjacency.rds")

y <- apply(y, 1:4, as.numeric)
x <- apply(x, 1:4, as.numeric)

# spatial effects
arealeffs <- c('sqrt.perc.agri', 'sqrt.perc.urban')
x <- append_spatial_effects(x, arealeffs, A)

sites <- dimnames(x)[[1]]
seasons <- dimnames(y)[[2]]
periods <- 18:44

niter <- 6000
nchains <- 2
print.interval <- 500

# AA ----------------------------------------------------------------------

occueffs <- c("intercept")
deteffs <- c("intercept")
AA <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(AA, "models/A/AA.rds")

# AB ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
AB <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(AB, "models/A/AB.rds")

# AC ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
AC <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(AC, "models/A/AC.rds")

# AD ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
AD <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(AD, "models/A/AD.rds")

# AE ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag6sqrt.level.avg.avg",
              "lag6gini.bird.simpson")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
AE <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(AE, "models/A/AE.rds")

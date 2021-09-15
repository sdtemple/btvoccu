# btvoccu Models A*
# Seth Temple, sdtemple@lanl.gov
# September 11, 2021

library(btvoccu)
set.seed(9112021)

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

# AF ----------------------------------------------------------------------

# occueffs <- c("intercept",
#               "sqrt.perc.agri",
#               "sqrt.perc.urban",
#               "sqrt.wks.meantemp.below.freezing",
#               "lag6meantemp.wk",
#               "lag6sqrt.level.avg.avg",
#               "lag6gini.bird.simpson",
#               "lag1infected.neighbors")
# deteffs <- c("intercept",
#              "sqrt.X.surveys",
#              "sqrt.Cx.prop",
#              "cbrt.catch.rate")
# AF <- btvoccu(niter,
#               x, y,
#               sites,
#               seasons,
#               periods,
#               occueffs,
#               deteffs,
#               nchains = nchains,
#               print.interval = print.interval)
# saveRDS(AF, "models/A/AF.rds")

# AF2 ---------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag6sqrt.level.avg.avg",
              "lag6gini.bird.simpson",
              "lag2infected.neighbors")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
AF2 <- btvoccu(niter,
               x, y,
               sites,
               seasons,
               periods,
               occueffs,
               deteffs,
               nchains = nchains,
               print.interval = print.interval)       
saveRDS(AF2, "models/A/AF2.rds")



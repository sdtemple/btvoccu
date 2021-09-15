# btvoccu Models B*
# Seth Temple, sdtemple@lanl.gov
# August 26, 2021

library(btvoccu)
set.seed(8262021)

# Setup -------------------------------------------------------------------

y <- readRDS("data/response-train.rds")
x <- readRDS("data/final-covariates.rds")

y <- apply(y, 1:4, as.numeric)
x <- apply(x, 1:4, as.numeric)

sites <- dimnames(x)[[1]]
seasons <- dimnames(y)[[2]]
periods <- 18:44

niter <- 2000
nchains <- 2
print.interval <- 200

# BA ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag2meantemp.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
BA <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(BA, "models/B/BA.rds")

# BB ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag4meantemp.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
BB <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(BB, "models/B/BB.rds")

# BC ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
BC <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(BC, "models/B/BC.rds")

# BD ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag8meantemp.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
BD <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(BD, "models/B/BD.rds")

# Model evaluation --------------------------------------------------------

z <- readRDS("data/response-validate.rds")
z <- apply(z, 1:4, as.numeric)

BA <- readRDS("models/B/BA.rds")
BB <- readRDS("models/B/BB.rds")
BC <- readRDS("models/B/BC.rds")
BD <- readRDS("models/B/BD.rds")

table <- matrix(nrow = 4, ncol = 5)

# BA
table[1,1] <- waic_score(BA, x, y)$waic
table[1,2] <- waic_score(BA, x, z)$waic
table[1,3:5] <- posterior_effects(BA)[5,2:4]

# BB
table[2,1] <- waic_score(BB, x, y)$waic
table[2,2] <- waic_score(BB, x, z)$waic
table[2,3:5] <- posterior_effects(BB)[5,2:4]

# BC
table[3,1] <- waic_score(BC, x, y)$waic
table[3,2] <- waic_score(BC, x, z)$waic
table[3,3:5] <- posterior_effects(BC)[5,2:4]

# BD
table[4,1] <- waic_score(BD, x, y)$waic
table[4,2] <- waic_score(BD, x, z)$waic
table[4,3:5] <- posterior_effects(BD)[5,2:4]

colnames(table) <- c("trWAIC", "valWAIC", "0.025", "0.50", "0.975")
rownames(table) <- c("2","4","6","8")
saveRDS(table, "models/B/Btable.rds")

# share
table <- readRDS("models/B/Btable.rds")
table <- apply(table, 1:2, as.numeric)
table <- apply(table, 1:2, round, digits = 3)
library(xtable)
xtable(table, digits = 3)

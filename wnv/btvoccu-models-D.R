# btvoccu Models D*
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

# DA ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag6sqrt.level.avg.avg",
              "lag2gini.bird.simpson")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
DA <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(DA, "models/D/DA.rds")

# DB ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag6sqrt.level.avg.avg",
              "lag4gini.bird.simpson")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
DB <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(DB, "models/D/DB.rds")

# DC ----------------------------------------------------------------------

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
DC <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(DC, "models/D/DC.rds")

# DD ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag6sqrt.level.avg.avg",
              "lag8gini.bird.simpson")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
DD <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(DD, "models/D/DD.rds")

# Model evaluation --------------------------------------------------------

z <- readRDS("data/response-validate.rds")
z <- apply(z, 1:4, as.numeric)

DA <- readRDS("models/D/DA.rds")
DB <- readRDS("models/D/DB.rds")
DC <- readRDS("models/D/DC.rds")
DD <- readRDS("models/D/DD.rds")

table <- matrix(nrow = 4, ncol = 5)

# DA
table[1,1] <- waic_score(DA, x, y)$waic
table[1,2] <- waic_score(DA, x, z)$waic
table[1,3:5] <- posterior_effects(DA)[6,2:4]

# DB
table[2,1] <- waic_score(DB, x, y)$waic
table[2,2] <- waic_score(DB, x, z)$waic
table[2,3:5] <- posterior_effects(DB)[6,2:4]

# DC
table[3,1] <- waic_score(DC, x, y)$waic
table[3,2] <- waic_score(DC, x, z)$waic
table[3,3:5] <- posterior_effects(DC)[6,2:4]

# DD
table[4,1] <- waic_score(DD, x, y)$waic
table[4,2] <- waic_score(DD, x, z)$waic
table[4,3:5] <- posterior_effects(DD)[6,2:4]

colnames(table) <- c("trWAIC", "valWAIC", "0.025", "0.50", "0.975")
rownames(table) <- c("2","4","6","8")
saveRDS(table, "models/D/Dtable.rds")

# share
table <- readRDS("models/D/Dtable.rds")
table <- apply(table, 1:2, as.numeric)
table <- apply(table, 1:2, round, digits = 3)
library(xtable)
xtable(table, digits = 3)

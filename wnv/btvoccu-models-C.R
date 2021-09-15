# btvoccu Models C*
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

# CA ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag2precip.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
CA <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(CA, "models/C/CA.rds")

# CB ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag4precip.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
CB <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(CB, "models/C/CB.rds")

# CC ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag6precip.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
CC <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(CC, "models/C/CC.rds")

# CD ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag8precip.wk")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
CD <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(CD, "models/C/CD.rds")

# CE ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag2sqrt.level.avg.avg")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
CE <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(CE, "models/C/CE.rds")

# CF ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag4sqrt.level.avg.avg")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
CF <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(CF, "models/C/CF.rds")

# CG ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag6sqrt.level.avg.avg")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
CG <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(CG, "models/C/CG.rds")

# CH ----------------------------------------------------------------------

occueffs <- c("intercept",
              "sqrt.perc.agri",
              "sqrt.perc.urban",
              "sqrt.wks.meantemp.below.freezing",
              "lag6meantemp.wk",
              "lag8sqrt.level.avg.avg")
deteffs <- c("intercept",
             "sqrt.X.surveys",
             "sqrt.Cx.prop",
             "cbrt.catch.rate")
CH <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(CH, "models/C/CH.rds")

# Model evaluation --------------------------------------------------------

z <- readRDS("data/response-validate.rds")
z <- apply(z, 1:4, as.numeric)

CA <- readRDS("models/C/CA.rds")
CB <- readRDS("models/C/CB.rds")
CC <- readRDS("models/C/CC.rds")
CD <- readRDS("models/C/CD.rds")
CE <- readRDS("models/C/CE.rds")
CF <- readRDS("models/C/CF.rds")
CG <- readRDS("models/C/CG.rds")
CH <- readRDS("models/C/CH.rds")

table <- matrix(nrow = 8, ncol = 5)

# CA
table[1,1] <- waic_score(CA, x, y)$waic
table[1,2] <- waic_score(CA, x, z)$waic
table[1,3:5] <- posterior_effects(CA)[6,2:4]

# CB
table[2,1] <- waic_score(CB, x, y)$waic
table[2,2] <- waic_score(CB, x, z)$waic
table[2,3:5] <- posterior_effects(CB)[6,2:4]

# CC
table[3,1] <- waic_score(CC, x, y)$waic
table[3,2] <- waic_score(CC, x, z)$waic
table[3,3:5] <- posterior_effects(CC)[6,2:4]

# CD
table[4,1] <- waic_score(CD, x, y)$waic
table[4,2] <- waic_score(CD, x, z)$waic
table[4,3:5] <- posterior_effects(CD)[6,2:4]

# CE
table[5,1] <- waic_score(CE, x, y)$waic
table[5,2] <- waic_score(CE, x, z)$waic
table[5,3:5] <- posterior_effects(CE)[6,2:4]

# CF
table[6,1] <- waic_score(CF, x, y)$waic
table[6,2] <- waic_score(CF, x, z)$waic
table[6,3:5] <- posterior_effects(CF)[6,2:4]

# CG
table[7,1] <- waic_score(CG, x, y)$waic
table[7,2] <- waic_score(CG, x, z)$waic
table[7,3:5] <- posterior_effects(CG)[6,2:4]

# CH
table[8,1] <- waic_score(CH, x, y)$waic
table[8,2] <- waic_score(CH, x, z)$waic
table[8,3:5] <- posterior_effects(CH)[6,2:4]

colnames(table) <- c("trWAIC", "valWAIC", "0.025", "0.50", "0.975")
rownames(table) <- c("2 (precip)","4","6","8", "2 (level)", "4", "6", "8")
saveRDS(table, "models/C/Ctable.rds")

# share
table <- readRDS("models/C/old/Ctable.rds")
table <- apply(table, 1:2, as.numeric)
table <- apply(table, 1:2, round, digits = 3)
library(xtable)
xtable(table, digits = 3)

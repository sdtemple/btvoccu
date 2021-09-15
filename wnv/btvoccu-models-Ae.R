# btvoccu Models A* (Evaluation)
# Seth Temple, sdtemple@lanl.gov
# September 13, 2021

library(btvoccu)
set.seed(9132021)

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

# Model evaluation --------------------------------------------------------

z <- readRDS("data/response-validate.rds")
z <- apply(z, 1:4, as.numeric)

AA <- readRDS("models/A/AA.rds")
AB <- readRDS("models/A/AB.rds")
AC <- readRDS("models/A/AC.rds")
AD <- readRDS("models/A/AD.rds")
AE <- readRDS("models/A/AE.rds")
AF <- readRDS("models/A/AF2.rds")
AG <- readRDS("models/A/AG5.rds")

table <- matrix(nrow = 7, ncol = 5)
gta <- c("TOR", "YRK", "PEE", "HAL", "DUR")
nrd <- c("ALG", "NPS", "NWR", "PQP", "SUD", "THB", "TSK")

# AA
table[1,1] <- waic_score(AA, x, y, burnin = 1/6)$waic
table[1,2] <- waic_score(AA, x, z, burnin = 1/6)$waic
table[1,3] <- posterior_check(AA, x, z, burnin = 1/6)$relativePresence
table[1,4] <- posterior_check(AA, x, z, gta, burnin = 1/6)$relativePresence
table[1,5] <- posterior_check(AA, x, z, nrd, burnin = 1/6)$relativePresence

# AB
table[2,1] <- waic_score(AB, x, y, burnin = 1/6)$waic
table[2,2] <- waic_score(AB, x, z, burnin = 1/6)$waic
table[2,3] <- posterior_check(AB, x, z, burnin = 1/6)$relativePresence
table[2,4] <- posterior_check(AB, x, z, gta, burnin = 1/6)$relativePresence
table[2,5] <- posterior_check(AB, x, z, nrd, burnin = 1/6)$relativePresence

# AC
table[3,1] <- waic_score(AC, x, y, burnin = 1/6)$waic
table[3,2] <- waic_score(AC, x, z, burnin = 1/6)$waic
table[3,3] <- posterior_check(AC, x, z, burnin = 1/6)$relativePresence
table[3,4] <- posterior_check(AC, x, z, gta, burnin = 1/6)$relativePresence
table[3,5] <- posterior_check(AC, x, z, nrd, burnin = 1/6)$relativePresence

# AD
table[4,1] <- waic_score(AD, x, y, burnin = 1/6)$waic
table[4,2] <- waic_score(AD, x, z, burnin = 1/6)$waic
table[4,3] <- posterior_check(AD, x, z, burnin = 1/6)$relativePresence
table[4,4] <- posterior_check(AD, x, z, gta, burnin = 1/6)$relativePresence
table[4,5] <- posterior_check(AD, x, z, nrd, burnin = 1/6)$relativePresence

# AE
table[5,1] <- waic_score(AE, x, y, burnin = 1/6)$waic
table[5,2] <- waic_score(AE, x, z, burnin = 1/6)$waic
table[5,3] <- posterior_check(AE, x, z, burnin = 1/6)$relativePresence
table[5,4] <- posterior_check(AE, x, z, gta, burnin = 1/6)$relativePresence
table[5,5] <- posterior_check(AE, x, z, nrd, burnin = 1/6)$relativePresence

# AF
table[6,1] <- waic_score(AF, x, y, burnin = 1/6)$waic
table[6,2] <- waic_score(AF, x, z, burnin = 1/6)$waic
table[6,3] <- posterior_check(AF, x, z, burnin = 1/6)$relativePresence
table[6,4] <- posterior_check(AF, x, z, gta, burnin = 1/6)$relativePresence
table[6,5] <- posterior_check(AF, x, z, nrd, burnin = 1/6)$relativePresence

# AG
table[7,1] <- waic_score(AG, x, y, burnin = 1/2)$waic
table[7,2] <- waic_score(AG, x, z, burnin = 1/2)$waic
table[7,3] <- posterior_check(AG, x, z, burnin = 1/2)$relativePresence
table[7,4] <- posterior_check(AG, x, z, gta, burnin = 1/2)$relativePresence
table[7,5] <- posterior_check(AG, x, z, nrd, burnin = 1/2)$relativePresence


colnames(table) <- c("trWAIC", "valWAIC", "allRelPre", "gtaRelPre", "nrdRelPre")
rownames(table) <- c("AA","AB","AC","AD", "AE", "AF2", "AG5")
saveRDS(table, "models/A/Atable.rds")

# share
table <- readRDS("models/A/Atable.rds")
table <- apply(table, 1:2, as.numeric)
table <- apply(table, 1:2, round, digits = 3)
library(xtable)
xtable(table, digits = 3)

# btvoccu Final Model
# Seth Temple, sdtemple@lanl.gov
# September 14, 2021

library(btvoccu)
set.seed(9142021)

# Setup -------------------------------------------------------------------

y <- readRDS("data/response-split.rds")
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

niter <- 12000
nchains <- 3
print.interval <- 1000

# AF ----------------------------------------------------------------------

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
AF <- btvoccu(niter,
              x, y,
              sites,
              seasons,
              periods,
              occueffs,
              deteffs,
              nchains = nchains,
              print.interval = print.interval)
saveRDS(AF, "models/F for final/AF2-9142021.rds")

# Tabular -----------------------------------------------------------------

AF <- readRDS("models/F for final/AF2-9142021.rds")
y <- readRDS("data/response-split.rds")
z <- readRDS("data/response-test.rds")
x <- readRDS("data/final-covariates.rds")
y <- apply(y, 1:4, as.numeric)
z <- apply(z, 1:4, as.numeric)
x <- apply(x, 1:4, as.numeric)

phu <- dimnames(z)[[1]]
s <- dimnames(z)[[2]]
rp <- array(NA, dim = c(length(phu), 7))
for(h in 1:length(phu)){
  print(h)
  rp[h,1] <- phu[h]
  out <- posterior_check(AF, x, z, sites = phu[h], seasons = s[1], periods = 18:44, burnin = 1/6) # 2008
  rp[h,2] <- out$avgPresence
  rp[h,3] <- out$actualPresence
  out <- posterior_check(AF, x, z, sites = phu[h], seasons = s[2], periods = 18:44, burnin = 1/6) # 2012
  rp[h,4] <- out$avgPresence
  rp[h,5] <- out$actualPresence
  out <- posterior_check(AF, x, z, sites = phu[h], seasons = s[3], periods = 18:44, burnin = 1/6) # 2016
  rp[h,6] <- out$avgPresence
  rp[h,7] <- out$actualPresence
}
rpn <- apply(rp[,2:7], 1:2, as.numeric) # numeric columns only
rsm <- rpn[,1] + rpn[,3] + rpn[,5] # row sum model averages
rsa <- rpn[,2] + rpn[,4] + rpn[,6] # row sum actual presences
cs <- colSums(rpn) # column sum

# append to table
table <- cbind(rp, rsm)
table <- cbind(table, rsa)
table <- rbind(table, c("Total", cs, sum(cs[c(1,3,5)]), sum(cs[c(2,4,6)])))

# table for latex
library(xtable)
table <- table[,2:9]
rownames(table) <- c(phu, "Total")
xtable(table)

# posterior effects
pe <- posterior_effects(AF, burnin = 1/6)
out <- apply(pe[,2:6], 1:2, as.numeric)
out <- apply(out, 1:2, round, digits = 3)
pe[,2:6] <- out 
xtable(pe) # table for latex

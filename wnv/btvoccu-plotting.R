# btvoccu time series and maps
# Seth Temple, sdtemple@lanl.gov
# September 14, 2021

set.seed(9142021)

library(btvoccu)
library(btvoccumaps)
library(rgdal)
library(gridExtra)
library(magick)

final <- readRDS("models/F for final/AF2-9142021.rds")
x <- readRDS("data/final-covariates.rds")
x <- apply(x, 1:4, as.numeric)

sites <- c("TOR","OTT","SUD","WEC")
moresites <- c(sites, "ALG", "THB", "YRK", "PEE")
season <- 2016
periods <- 20:40

# Trace plots -------------------------------------------------------------

plot_trace <- function(model, effects, mnridx, xlab, ylab, v, burnin = 0){
  
  mjridx <- which(names(model) == effects, arr.ind = T)
  colors <- rcartocolor::carto_pal(12, "Safe")
  colors <- colors[1:model$nchains]
  #colors <- grDevices::rainbow(model$nchains)
  burn <- -(1:(floor(burnin * model$niter)))
  
  # ylim
  mx <- max(model[[mjridx]][1, burn, mnridx])
  mn <- min(model[[mjridx]][1, burn, mnridx])
  for(n in 2:model$nchains){
    mx <- max(mx, model[[mjridx]][n, burn, mnridx])
    mn <- min(mn, model[[mjridx]][n, burn, mnridx])
  }
  
  plot(model[[mjridx]][1, burn, mnridx],
       ylab = ylab,
       xlab = xlab,
       type = "l",
       col = colors[1],
       ylim = c(mn, mx))
  abline(v = v, col = "black", lwd = 5)
  for(n in 2:model$nchains){lines(model[[mjridx]][n, burn, mnridx], col = colors[n])}
}

par(mfrow=c(2,2))

plot_trace(final, "betas", 1, "Index", "Intercept (occupancy)", 2000)
plot_trace(final, "betas", 2, "Index", "Agricultural", 2000)
plot_trace(final, "betas", 3, "Index", "Urban", 2000)
plot_trace(final, "betas", 4, "Index", "Freezing Weeks", 2000)

plot_trace(final, "betas", 5, "Index", "Temperature (6 week lag)", 2000)
plot_trace(final, "betas", 6, "Index", "Water Level (6 week lag)", 2000)
plot_trace(final, "betas", 7, "Index", "Gini Simpson Index (6 week lag)", 2000)
plot_trace(final, "betas", 8, "Index", "Infected Persons (2 week lag)", 2000)

plot_trace(final, "alphas", 1, "Index", "Intercept (detection)", 2000)
plot_trace(final, "alphas", 2, "Index", "Survey Count", 2000)
plot_trace(final, "alphas", 3, "Index", "Culex Proportion", 2000)
plot_trace(final, "alphas", 4, "Index", "Catch Rate", 2000) # save these

# Time Series -------------------------------------------------------------

# not smoothed
a <- plot_btvoccu(final, x, sites, season, periods,
                  value = "Occupancy",
                  burnin = 1/6)
b <- plot_btvoccu(final, x, sites, season, periods,
                  value = "Detection",
                  burnin = 1/6)
c <- plot_btvoccu(final, x, sites, season, periods,
                  value = "Presence",
                  burnin = 1/6)
d <- plot_btvoccu(final, x, moresites, season, periods,
                  value = "Occupancy",
                  burnin = 1/6,
                  ribbons = F)
grid.arrange(a, b, c, d, nrow = 2, ncol = 2) # save this


# smoothed
a <- plot_btvoccu(final, x, sites, season, periods,
                  value = "Occupancy",
                  burnin = 1/6,
                  spline = T)
b <- plot_btvoccu(final, x, sites, season, periods,
                  value = "Detection",
                  burnin = 1/6,
                  spline = T)
c <- plot_btvoccu(final, x, sites, season, periods,
                  value = "Presence",
                  burnin = 1/6, 
                  spline = T)
d <- plot_btvoccu(final, x, moresites, season, periods,
                  value = "Occupancy",
                  burnin = 1/6,
                  spline = T,
                  ribbons = F)
grid.arrange(a, b, c, d, nrow = 2, ncol = 2) # save this

# covariates
a <- plot_covariate(x, moresites, season, 20:40, "lag6gini.bird.simpson", "Bird Index (6 week lag)")
b <- plot_covariate(x, moresites, season, 20:40, "lag6meantemp.wk", "Temperature (6 week lag)")
c <- plot_covariate(x, moresites, season, 20:40, "lag2infected.neighbors", "Infected Neighbors (2 week lag)")
d <- plot_covariate(x, moresites, season, 20:40, "lag6sqrt.level.avg.avg", "Water Level (6 week lag)")
grid.arrange(b, c, d, a, nrow = 2, ncol = 2) # save this

a <- plot_covariate(x, moresites, season, 20:40, "sqrt.X.surveys", "Survey Count (sqrt)")
b <- plot_covariate(x, moresites, season, 20:40, "sqrt.Cx.prop", "Percent Culex (sqrt)")
c <- plot_covariate(x, moresites, season, 20:40, "cbrt.catch.rate", "Catch Rate (cbrt)")
grid.arrange(a, b, c, nrow = 2, ncol = 2) # save this

# Maps --------------------------------------------------------------------

# shapefile
ontario <- readOGR('data/Ministry_of_Health_Public_Health_Unit_Boundary-shp/Ministry_of_Health_Public_Health_Unit_Boundary.shp')

# simple features
cities <- readRDS("data/cities.rds")

# in the order of the shapefile
phus <- c('YRK','HUR','WAT','ELG','HAM','THB','PEE','LAM','WDG','BRN','MSL',
          'SUD','HKP','NIA','CHK','KFL','WEC','PTC','GBO','EOH','NPS','OTT',
          'LGL','NWR','HDN','TSK','REN','TOR','HAL','HPE','SMD','DUR','PQP','ALG')

btvoccu_map_inset("figures/occumap.png",
                  final, x, "Occupancy", 0.5, ontario, cities, 
                  phus, season, 34, latmax = 46)
btvoccu_map_inset("figures/detmap.png",
                  final, x, "Detection", 0.5, ontario, cities, 
                  phus, season, 34, latmax = 46)

# covariates
covariate_map_inset("figures/Covariates/maps/agri-inset.png",
                    x, "sqrt.perc.agri", ontario, cities, phus,
                    2017, 30, "Sqrt Percent Agricultural",
                    latmax= 46)
# covariate_map_inset("figures/Covariates/maps/urban-inset.png",
#                     x, "sqrt.perc.urban", ontario, cities, phus,
#                     2017, 30, "Sqrt Percent Urban",
#                     latmax= 46)

# occupancy movie
for(p in 20:40){
  btvoccu_map_inset(paste("figures/occu2016/frame", p, ".png", sep = ""),
                          final, x, "Occupancy", 0.5, ontario, cities, phus, 
                          season, p, latmax = 46)
}
imgs <- list.files("figures/occu2016/", full.names = T)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 1)
image_write(image = img_animated,
            path = "figures/occu2016/occu2016movie.gif")

# detection movie
for(p in 20:40){
  btvoccu_map_inset(paste("figures/det2016/frame", p, ".png", sep = ""),
                    final, x, "Detection", 0.5, ontario, cities, phus, 
                    season, p, latmax = 46)
}
imgs <- list.files("figures/det2016/", full.names = T)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 1)
image_write(image = img_animated,
            path = "figures/det2016/det2016movie.gif")

# presence movie
for(p in 20:40){
  btvoccu_map_inset(paste("figures/pre2016/frame", p, ".png", sep = ""),
                    final, x, "Presence", 0.5, ontario, cities, phus, 
                    season, p, latmax = 46)
}
imgs <- list.files("figures/pre2016/", full.names = T)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 1)
image_write(image = img_animated,
            path = "figures/pre2016/pre2016movie.gif")

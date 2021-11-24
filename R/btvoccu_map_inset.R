#' Thematic \code{btvoccu} Maps with Insets
#'
#' Create zoomed in thematic map demonstrating
#' how output sigmoids vary across space.
#'
#' @param file output file name
#' @param model \code{btvoccu} model object
#' @param x numeric covariate array
#' @param value "Occupancy", "Detection", or "Presence"
#' @param q quantile
#' @param shape \code{SpatialPolygonsDataFrame} object
#' @param cities \code{sf} object, named points
#' @param sites vector (in order of \code{shape} object)
#' @param season integer
#' @param period integer
#' @param latmax numeric
#' @param latmin numeric
#' @param longmax numeric
#' @param longmin numeric
#' @param palette color palette
#' @param inset.x see \code{grid::viewport()}
#' @param inset.y see \code{grid::viewport()}
#' @param inset.width see \code{grid::viewport()}
#' @param legend.position see \code{tmap::tm_layout()}
#' @param title.position see \code{tmap::tm_layout()}
#' @param title.size see \code{tmap::tm_layout()}
#' @param legend.title.size see \code{tmap::tm_layout()}
#' @param legend.text.size see \code{tmap::tm_layout()}
#' @param xlab.size see \code{tmap::tm_xlab()}
#' @param ylab.size see \code{tmap::tm_ylab()}
#' @param dots.size see \code{tmap::tm_dots()}
#' @param cities.text see \code{tmap::tm_text()}
#' @param cities.just see \code{tmap::tm_text()}
#' @param cities.size see \code{tmap::tm_text()}
#' @param cities.col see \code{tmap::tm_text()}
#'
#' @return \code{tmap} object
#'
#' @export
btvoccu_map_inset <- function(file, model, x, 
                              value, q,
                              shape, cities, 
                              sites, season, period,
                              latmax = 90, latmin = -90,
                              longmax = 180, longmin = -180,
                              palette = 'Reds', 
                              inset.x = .95, 
                              inset.y = .1, 
                              inset.width = .35,
                              legend.position = c('left','top'), 
                              title.position  = c('right','top'),
                              title.size = .8, 
                              legend.title.size = .8, 
                              legend.text.size = .6,
                              xlab.size = .8, 
                              ylab.size = .8,
                              dots.size = .04,
                              cities.text = 'name', 
                              cities.just = 'top',
                              cities.size = .8, 
                              cities.col = 'black'
                              ){

  # organize data
  sigmoids <- posterior_sigmoids(model, x, sites, season, (period - 1):period) # model estimates
  shape@data$site <- sites
  arry <- array(NA, dim = c(length(sites), 2), dimnames = list(sites, c('site','value')))
  arry[,1] <- sites
  idx <- which(names(sigmoids) %in% value, arr.ind = T)
  sigmoid <- sigmoids[[idx]]
  sigmoid <- subset_4darray(sigmoid, 3, period)
  sigmoid <- sigmoid[sites,,,, drop = F] # reorder sites
  arry[,2] <- apply(sigmoid, 1:3, quantile, probs = q, na.rm = T)

  # merge shape
  merged.shape <- sp::merge(shape, arry, by.x = 'site', by.y = 'site')
  merged.shape@data$value <- as.numeric(as.character(merged.shape@data$value))

  # zoomed in shape and cities
  centroids <- rgeos::gCentroid(merged.shape, byid = T)
  merged.shape@data$long <-  centroids@coords[,1]
  merged.shape@data$lat <- centroids@coords[,2]
  focused.shape <- merged.shape[merged.shape$long > longmin,]
  focused.shape <- focused.shape[focused.shape$long < longmax,]
  focused.shape <- focused.shape[focused.shape$lat > latmin,]
  focused.shape <- focused.shape[focused.shape$lat < latmax,]
  cities$long <- sf::st_coordinates(cities)[,1]
  cities$lat <- sf::st_coordinates(cities)[,2]
  focused.cities <- cities[cities$long > longmin,]
  focused.cities <- focused.cities[focused.cities$long < longmax,]
  focused.cities <- focused.cities[focused.cities$lat > latmin,]
  focused.cities <- focused.cities[focused.cities$lat < latmax,]

  # main tmap object
  plt <- tmap::tm_shape(focused.shape) +
    tmap::tm_polygons(col = 'value',
                      title = value,
                      breaks = seq(0, 1, .2),
                      style = 'cont',
                      palette = palette) +
    tmap::tm_grid(lines = F) +
    tmap::tm_layout(title = paste('Epiweek ', period, ', ', season, sep = ''),
                    title.size = title.size,
                    title.position = title.position,
                    legend.position = legend.position,
                    legend.title.size = legend.title.size,
                    legend.text.size = legend.text.size,
                    asp = 0) +
    tmap::tm_xlab("Longitude", size = xlab.size) +
    tmap::tm_ylab("Latitude", size = ylab.size) +
    tmap::tm_shape(focused.cities) +
    tmap::tm_dots(size = dots.size, col = cities.col) +
    tmap::tm_text(cities.text, just = cities.just, 
                  size = cities.size, col = cities.col)

  # inset tmap object
  inset <- tmap::tm_shape(merged.shape) +
    tmap::tm_polygons(col = 'value',
                      title = value,
                      breaks = seq(0, 1, .2),
                      style = 'cont',
                      palette = palette) +
    tmap::tm_grid(lines = F) +
    tmap::tm_layout(legend.show = F, asp = 0)

  # define viewport and manage aspect ratio
  w <- inset.width
  h <- tmaptools::get_asp_ratio(inset) * inset.width
  vp <- grid::viewport(x = inset.x,
                       y = inset.y,
                       width = w,
                       height = h,
                       just = c("right", "bottom"))

  # saving
  tmap::tmap_save(plt, file,
                  insets_tm = inset,
                  insets_vp = vp,
                  height = tmaptools::get_asp_ratio(plt) * 5,
                  width = 8.5 - 2 * 1,
                  units = "in")
}

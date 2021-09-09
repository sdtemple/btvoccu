#' Thematic \code{btvoccu} Maps
#'
#' Create thematic map demonstrating how
#' output sigmoids vary across space.
#'
#' @param model \code{btvoccu} model object
#' @param x numeric covariate array
#' @param value "Occupancy", "Detection", or "Presence"
#' @param q quantile
#' @param shape \code{SpatialPolygonsDataFrame} object
#' @param cities \code{sf} object, named points
#' @param sites vector (in order of \code{shape} object)
#' @param season integer
#' @param period integer
#' @param palette color palette
#' @param legend.position see \code{tmap::tm_layout()}
#' @param title.position see \code{tmap::tm_layout()}
#' @param compass.position see \code{tmap::tm_compass()}
#' @param title.size see \code{tmap::tm_layout()}
#' @param legend.title.size see \code{tmap::tm_layout()}
#' @param legend.text.size see \code{tmap::tm_layout()}
#'
#' @return \code{tmap} object
#'
#' @export
btvoccu_map <- function(file, model, x, value, q,
                        shape, cities, sites, season,period,
                        palette = 'Reds', legend.position = c('left','bottom'), 
                        title.position = c('right','top'), compass.position = c('right','bottom'),
                        title.size = .8, legend.title.size = .8, legend.text.size = .6
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

  # making a map
  plt <- tmap::tm_shape(merged.shape) +
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
          tmap::tm_compass(position = compass.position) +
          tmap::tm_shape(cities) +
          tmap::tm_dots(size = .2,
                        col = 'black') +
          tmap::tm_text('name', just = 'top', size = .6)

  # saving
  tmap::tmap_save(plt, file, height = tmaptools::get_asp_ratio(plt) * 5, width = 5, units = "in")
}

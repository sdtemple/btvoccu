#' Thematic Covariate Maps with Insets
#'
#' Create a thematic map demonstrating how
#' a covariate varies across space.
#'
#' @param file output file name
#' @param x numeric covariate array
#' @param covariate character name of covariate
#' @param shape \code{SpatialPolygonsDataFrame} object
#' @param cities \code{sf} object, named points
#' @param sites vector (in order of \code{shape} object)
#' @param season integer
#' @param period integer
#' @param label legend label
#' @param latmax numeric
#' @param latmin numeric
#' @param longmax numeric
#' @param longmin numeric
#' @param palette color palette
#' @param timestamp logical: default \code{FALSE} to not show timestamp
#' @param seasononly logical: default \code{FALSE} to show full timestamp
#' @param inset.x see \code{grid::viewport()}
#' @param inset.y see \code{grid::viewport()}
#' @param inset.width see \code{grid::viewport()}
#' @param legend.position see \code{tmap::tm_layout()}
#' @param title.position see \code{tmap::tm_layout()}
#' @param title.size see \code{tmap::tm_layout()}
#' @param legend.title.size see \code{tmap::tm_layout()}
#' @param legend.text.size see \code{tmap::tm_layout()}
#'
#' @return \code{tmap} object
#'
#' @export
covariate_map_inset <- function(file, x, covariate,
                                shape, cities, sites, 
                                season, period, label,
                                latmax = 90, latmin = -90,
                                longmax = 180, longmin = -180,
                                palette = 'Reds',
                                timestamp = F, seasononly = F,
                                inset.x = .95, inset.y = .05, inset.width = .4,
                                legend.position = c('left','top'),
                                title.position = c('right','top'),
                                title.size = .8, 
                                legend.title.size = .8, 
                                legend.text.size = .6
                                ){

  # organize data
  shape@data$site <- sites
  j <- which(dimnames(x)[[2]] == season, arr.ind = T)
  k <- which(dimnames(x)[[3]] == period, arr.ind = T)
  l <- which(dimnames(x)[[4]] == covariate, arr.ind = T)
  arry <- data.frame(matrix(NA, ncol = 2, nrow = length(sites)))
  colnames(arry) <- c('site', 'value')
  arry[,1] <- sites
  arry[,2] <- x[sites, j, k, l]

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
  mx <- max(arry[,2], na.rm = T)
  mn <- min(arry[,2], na.rm = T)
  if(timestamp){ # control if timestamp is provided
    if(seasononly){ # control if period is listed
      plt <- tmap::tm_shape(focused.shape) +
        tmap::tm_polygons(col = 'value',
                          title = label,
                          breaks = seq(mn, mx, length.out = 5),
                          style = 'cont',
                          palette = palette) +
        tmap::tm_grid(lines = F) +
        tmap::tm_layout(title = season,
                        title.size = title.size,
                        title.position = title.position,
                        legend.position = legend.position,
                        legend.title.size = legend.title.size,
                        legend.text.size = legend.text.size,
                        asp = 0) +
        tmap::tm_shape(focused.cities) +
        tmap::tm_dots(size = .2, col = 'black') +
        tmap::tm_text('name', just  = 'top', size = .6)

    } else{
      plt <- tmap::tm_shape(focused.shape) +
        tmap::tm_polygons(col = 'value',
                          title = label,
                          breaks = seq(mn, mx, length.out = 5),
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
        tmap::tm_shape(focused.cities) +
        tmap::tm_dots(size = .2, col = 'black') +
        tmap::tm_text('name', just  = 'top', size = .6)
    }
  } else{
    plt <- tmap::tm_shape(focused.shape) +
      tmap::tm_polygons(col = 'value',
                        title = label,
                        breaks = seq(mn, mx, length.out = 5),
                        style = 'cont',
                        palette = palette) +
      tmap::tm_grid(lines = F) +
      tmap::tm_layout(legend.position = legend.position,
                      legend.title.size = legend.title.size,
                      legend.text.size = legend.text.size,
                      asp = 0) +
      tmap::tm_shape(focused.cities) +
      tmap::tm_dots(size = .2, col = 'black') +
      tmap::tm_text('name', just = 'top', size = .6)
  }

  # inset tmap object
  inset <- tmap::tm_shape(merged.shape) +
    tmap::tm_polygons(col = 'value',
                title = label,
                breaks = seq(mn, mx, length.out = 5),
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

  tmap::tmap_save(plt, file,
                  insets_tm = inset, 
                  insets_vp = vp,
                  height = tmaptools::get_asp_ratio(plt)*5,
                  width=5,
                  units="in")
}

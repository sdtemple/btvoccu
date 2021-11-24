#' Thematic Covariate Maps
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
#' @param palette color palette
#' @param timestamp logical: default \code{FALSE} to not show timestamp
#' @param seasononly logical: default \code{FALSE} to show full timestamp
#' @param legend.position see \code{tmap::tm_layout()}
#' @param title.position see \code{tmap::tm_layout()}
#' @param compass.position see \code{tmap::tm_compass()}
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
#' @export
covariate_map <- function(file, x, covariate,
                          shape, cities, 
                          sites, season, period,
                          label, 
                          palette = 'Reds',
                          timestamp = FALSE, 
                          seasononly = FALSE,
                          legend.position = c('left','bottom'),
                          title.position = c('right','top'),
                          compass.position = c('right','bottom'),
                          title.size = .8, 
                          legend.title.size = .8,
                          legend.text.size = .6,
                          xlab.size = .8, 
                          ylab.size = .8,
                          dots.size = .04,
                          cities.text = 'name', 
                          cities.just = 'top',
                          cities.size = .7, 
                          cities.col = 'black'
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

  # making a map
  if(timestamp){ # control if timestamp is provided
    if(seasononly){ # control if period is listed
      plt <- tmap::tm_shape(merged.shape) +
              tmap::tm_polygons(col = 'value',
                                title = label,
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
              tmap::tm_xlab("Longitude", size = xlab.size) +
              tmap::tm_ylab("Latitude", size = ylab.size) + 
              tmap::tm_compass(position = compass.position) +
              tmap::tm_shape(cities) +
              tmap::tm_dots(size = dots.size, col = cities.col) +
              tmap::tm_text(cities.text, just = cities.just,
                            size = cities.size, col = cities.col)
    } else{
      plt <- tmap::tm_shape(merged.shape) +
              tmap::tm_polygons(col = 'value',
                                title = label,
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
              tmap::tm_compass(position = compass.position) +
              tmap::tm_shape(cities) +
              tmap::tm_dots(size = dots.size, col = cities.col) +
              tmap::tm_text(cities.text, just = cities.just, 
                            size = cities.size, col = cities.col)
    }
  } else{
    plt <- tmap::tm_shape(merged.shape) +
            tmap::tm_polygons(col = 'value',
                              title = label,
                              style = 'cont',
                              palette = palette) +
            tmap::tm_grid(lines = F) +
            tmap::tm_layout(legend.position = legend.position,
                            legend.title.size = legend.title.size,
                            legend.text.size = legend.text.size,
                            asp = 0) +
            tmap::tm_xlab("Longitude", size = xlab.size) +
            tmap::tm_ylab("Latitude", size = ylab.size) + 
            tmap::tm_compass(position = compass.position) +
            tmap::tm_shape(cities) +
            tmap::tm_dots(size = dots.size, col = cities.col) +
            tmap::tm_text(cities.text, just = cities.just, 
                          size = cities.size, col = cities.col)
  }

  # saving
  tmap::tmap_save(plt, file,
                  height = tmaptools::get_asp_ratio(plt) * 5,
                  width = 8.5 - 2 * 1,
                  units = "in")
}

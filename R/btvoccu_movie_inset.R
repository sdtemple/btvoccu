#' Thematic \code{btvoccu} Maps Animated (w/ Insets)
#'
#' Create zoomed in thematic movie demonstrating
#' how output sigmoids vary across space.
#'
#' @param file output file name
#' @param directory output directory name (must exist already)
#' @param fps integer, frames per second (1 is recommended)
#' @param model \code{btvoccu} model object
#' @param x numeric covariate array
#' @param value "Occupancy", "Detection", or "Presence"
#' @param q quantile
#' @param shape \code{SpatialPolygonsDataFrame} object
#' @param cities \code{sf} object, named points
#' @param sites vector (in order of \code{shape} object)
#' @param season integer
#' @param periods vector
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
#' @return movie of \code{tmap} objects
#'
#' @export
btvoccu_movie_inset <- function(file, directory, fps,
                                model, x, 
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
                              
  # make multiple maps
  for(p in periods){
    btvoccu_map_inset(paste(directory, "period", p, ".png", sep = ""), # store in directory
                      model, x, 
                      value, q, 
                      shape, cities, 
                      sites, season, p, 
                      latmax, latmin, 
                      longmax, longmin,
                      palette, 
                      inset.x, 
                      inset.y, 
                      inset.width, 
                      legend.position, 
                      title.position, 
                      title.size, 
                      legend.title.size, 
                      legend.text.size,
                      xlab.size, 
                      ylab.size,
                      dots.size,
                      cities.text, 
                      cities.just,
                      cities.size, 
                      cities.col)    
  }
                               
  # animate with magick package
  imgs <- list.files(directory, full.names = T)
  img_list <- lapply(imgs, image_read)
  img_joined <- magick::image_join(img_list)
  img_animated <- magick::image_animate(img_joined, fps = fps)
  magick::image_write(image = img_animated,
                      path = paste(directory, file, ".gif", sep = ""))
}

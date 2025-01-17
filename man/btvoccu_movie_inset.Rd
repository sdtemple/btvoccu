% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/btvoccu_movie_inset.R
\name{btvoccu_movie_inset}
\alias{btvoccu_movie_inset}
\title{Thematic \code{btvoccu} Maps Animated (w/ Insets)}
\usage{
btvoccu_movie_inset(
  file,
  directory,
  fps,
  model,
  x,
  value,
  q,
  shape,
  cities,
  sites,
  season,
  period,
  latmax = 90,
  latmin = -90,
  longmax = 180,
  longmin = -180,
  palette = "Reds",
  inset.x = 0.95,
  inset.y = 0.1,
  inset.width = 0.35,
  legend.position = c("left", "top"),
  title.position = c("right", "top"),
  title.size = 0.8,
  legend.title.size = 0.8,
  legend.text.size = 0.6,
  xlab.size = 0.8,
  ylab.size = 0.8,
  dots.size = 0.04,
  cities.text = "name",
  cities.just = "top",
  cities.size = 0.8,
  cities.col = "black"
)
}
\arguments{
\item{file}{output file name}

\item{directory}{output directory name (must exist already)}

\item{fps}{integer, frames per second (1 is recommended)}

\item{model}{\code{btvoccu} model object}

\item{x}{numeric covariate array}

\item{value}{"Occupancy", "Detection", or "Presence"}

\item{q}{quantile}

\item{shape}{\code{SpatialPolygonsDataFrame} object}

\item{cities}{\code{sf} object, named points}

\item{sites}{vector (in order of \code{shape} object)}

\item{season}{integer}

\item{latmax}{numeric}

\item{latmin}{numeric}

\item{longmax}{numeric}

\item{longmin}{numeric}

\item{palette}{color palette}

\item{inset.x}{see \code{grid::viewport()}}

\item{inset.y}{see \code{grid::viewport()}}

\item{inset.width}{see \code{grid::viewport()}}

\item{legend.position}{see \code{tmap::tm_layout()}}

\item{title.position}{see \code{tmap::tm_layout()}}

\item{title.size}{see \code{tmap::tm_layout()}}

\item{legend.title.size}{see \code{tmap::tm_layout()}}

\item{legend.text.size}{see \code{tmap::tm_layout()}}

\item{xlab.size}{see \code{tmap::tm_xlab()}}

\item{ylab.size}{see \code{tmap::tm_ylab()}}

\item{dots.size}{see \code{tmap::tm_dots()}}

\item{cities.text}{see \code{tmap::tm_text()}}

\item{cities.just}{see \code{tmap::tm_text()}}

\item{cities.size}{see \code{tmap::tm_text()}}

\item{cities.col}{see \code{tmap::tm_text()}}

\item{periods}{vector}
}
\value{
movie of \code{tmap} objects
}
\description{
Create zoomed in thematic movie demonstrating
how output sigmoids vary across space.
}

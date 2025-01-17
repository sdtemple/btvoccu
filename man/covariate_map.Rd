% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariate_map.R
\name{covariate_map}
\alias{covariate_map}
\title{Thematic Covariate Maps}
\usage{
covariate_map(
  file,
  x,
  covariate,
  shape,
  cities,
  sites,
  season,
  period,
  label,
  palette = "Reds",
  timestamp = FALSE,
  seasononly = FALSE,
  legend.position = c("left", "bottom"),
  title.position = c("right", "top"),
  compass.position = c("right", "bottom"),
  title.size = 0.8,
  legend.title.size = 0.8,
  legend.text.size = 0.6,
  xlab.size = 0.8,
  ylab.size = 0.8,
  dots.size = 0.04,
  cities.text = "name",
  cities.just = "top",
  cities.size = 0.7,
  cities.col = "black"
)
}
\arguments{
\item{file}{output file name}

\item{x}{numeric covariate array}

\item{covariate}{character name of covariate}

\item{shape}{\code{SpatialPolygonsDataFrame} object}

\item{cities}{\code{sf} object, named points}

\item{sites}{vector (in order of \code{shape} object)}

\item{season}{integer}

\item{period}{integer}

\item{label}{legend label}

\item{palette}{color palette}

\item{timestamp}{logical: default \code{FALSE} to not show timestamp}

\item{seasononly}{logical: default \code{FALSE} to show full timestamp}

\item{legend.position}{see \code{tmap::tm_layout()}}

\item{title.position}{see \code{tmap::tm_layout()}}

\item{compass.position}{see \code{tmap::tm_compass()}}

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
}
\description{
Create a thematic map demonstrating how
a covariate varies across space.
}

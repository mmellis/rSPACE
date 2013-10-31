\name{encounter.history}
\alias{encounter.history}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce a single encounter history file based on population parameters.  
}
\description{
  Given an underlying landscape, populates the landscape with simulated individuals based on spacing rules in \code{Parameters}.  Calculates probability of at least one individual at every location on the landscape and accumulates probability over grid cells to determine presence/absence.
}
\usage{
encounter.history(Parameters, map, grid_layer, n.cells)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Parameters}{
list of parameter values for simulation
}
  \item{map}{
a raster layer with underlying habitat information
}
  \item{grid_layer}{
a vector the length of the number of pixels in \code{map} indicating a grid cell index for each pixel 
}
  \item{n.cells}{
total number of grid cells in sampling frame
}
}
\details{
  Called through \code{create.landscapes}, but can be used as a standalone function to produce a single encounter history for testing.
}
\value{
a matrix with \code{n.cells} rows and \code{Parameters$n_yrs*Parameters$n_visits} columns with 1/0 indicating presence/absence of at least one individual in each cell at each visit.
}

\seealso{
  \code{\link{create.landscapes}}
}
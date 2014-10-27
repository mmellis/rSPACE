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
encounter.history(map, Parameters, ...)
}

\arguments{
  \item{Parameters}{
list of parameter values for simulation
}
  \item{map}{
a raster layer with underlying habitat information
}
  \item{...}{
additional options
 }
}
\details{
  Called through \code{create.landscapes}, but can be used as a standalone function to produce a single encounter history for testing.
  
  Use \code{showSteps = T} as an additional option to display plots and other output as error checking.
}
\value{
Vector with one element for each grid cell representing the encounter history at that cell. Encounter histories consist of 1/0 values indicating presence/absence of at least one individual in each cell at each visit in each year of the study period.
}

\seealso{
  \code{\link{create.landscapes}}
}
\name{find.power}
\alias{find.power}

\title{
Predict power from an rSPACE simulation and analysis
}
\description{
Uses a loess smoother to predict either the number of grid cells needed to achieve a given level of power or the expected power from a specific sample size based on rSPACE output.
}
\usage{
find.power(folder, pwr = 0.8, CI = 0.95, n_grid = NULL)
}

\arguments{
  \item{folder}{
 Full path for rSPACE folder containing both simulation and analysis output
}
  \item{pwr}{
  If specified, desired power
}
  \item{CI}{
  Significance level
}
  \item{n_grid}{
  If specified, number of grid cells used
}
}
\details{
Depending on whether \code{pwr} or \code{n_grid} are specified, \code{find.power} will return the opposite.  If both \code{pwr} and \code{n_grid}, \code{n_grid} will be ignored
}
\value{
Either the number of cells required to reach a given power or the power expected from a sample of a given number of cells.
}

\author{
Martha Ellis
}


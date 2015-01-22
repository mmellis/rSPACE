\name{findPower}
\alias{findPower}

\title{
Predict power from an rSPACE simulation and analysis
}
\description{
Uses a loess smoother to predict either the number of grid cells needed to achieve a given level of power or the expected power from a specific sample size based on rSPACE output.
}
\usage{
findPower(folder, data, CI = 0.95, pwr = 0.8, n_grid=NULL)
}

\arguments{
  \item{folder}{
 Full path for rSPACE folder containing both simulation and analysis output
}
  \item{data}{
 Dataframe with previously loaded results
 }
  \item{CI}{
  Significance level (one-tailed test)
} 
  \item{pwr}{
  If specified, desired power
}
  \item{n_grid}{
  If specified, number of grid cells used
}
}
\details{
Data must be provided either by specifying a folder within which to look for a results file or providing previously loaded results.

Depending on whether \code{pwr} or \code{n_grid} are specified, \code{findPower} will return the opposite.  If both \code{pwr} and \code{n_grid}, \code{n_grid} will be ignored
}
\value{
Either the number of cells required to reach a given power or the power expected from a sample of a given number of cells.
}

\author{
Martha Ellis
}


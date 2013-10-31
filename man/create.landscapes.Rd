\name{create.landscapes}
\alias{create.landscapes}
\title{Wrapper function to create replicate simulated landscapes}
\description{
Creates a series of files with simulated encounter histories for using in a mark-recapture analysis based on population parameters.
}
\usage{
create.landscapes(n_runs, map, filter.map = NULL, Parameters = NULL, base.name = "tmp", folder.dir=getwd())
}

\arguments{
  \item{n_runs}{
number of replicate encounter histories to produce
}
  \item{map}{
a raster layer with underlying habitat information
}
  \item{filter.map}{
secondary raster layer indicating areas to include in sample
}
  \item{Parameters}{
list of parameter values for simulation (see \code{\link{enter.parameters}})
}
  \item{base.name}{
prefix to use for encounter history files (e.g., the default "tmp" creates files named "tmpXX.txt" where XX indexes the number of the run)
}
  \item{folder.dir}{
full path of where to create output encounter history files
}
}
\details{
This function will produce \code{n_run} text files with simulated encounter histories.  Each individual file does not contain identifying information about the parameters used to create it, so it will be good practice to create a folder structure to keep separate simulations separate.  \code{folder.dir} must already exist.

\code{map} will define the landscape for the population simulation.  Values of \code{map} are used as weights for both selecting home range centers and scaling bivariate normal home ranges.

The optional \code{filter.map} defines the sampling frame for the population.  It must have either equal or smaller extent than \code{map}.  Areas to exclude from the sample should be marked as NA.    
}
\value{
list with directory location and list of output file names
}
\author{
  Martha Ellis
}
\references{ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).}

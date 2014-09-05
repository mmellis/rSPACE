\name{create.landscapes}
\alias{create.landscapes}
\title{Wrapper function to create replicate simulated landscapes}
\description{
Creates a series of files with simulated encounter histories for using in a mark-recapture analysis based on population parameters.
}
\usage{
create.landscapes(n_runs, map, Paramters, ... )
}

\arguments{
  \item{n_runs}{
number of replicate encounter histories to produce
}
  \item{map}{
a raster layer with underlying habitat information
}
  \item{Parameters}{
list of parameter values for simulation 
(see \code{\link{enter.parameters}})
}
  \item{...}{
additional options
}
}
\details{
This function will produce \code{n_run} text files with simulated encounter histories.  Each individual file does not contain identifying information about the parameters used to create it, so it will be good practice to create a folder structure to keep separate simulations separate.  Default is to create a folder for the runs (./rSPACE_X) and name all encounter histories ./rSPACE_X/rSPACEx??.txt where ?? indicates the run number.  If \code{Paramters} is not specified, \code{\link{enter.parameters}} is called by default.

\code{map} will define the landscape for the population simulation.  Values of \code{map} are used as weights for both selecting home range centers and scaling bivariate normal home ranges.  

There are additional arguments that can be included to specify other options, but some are more sketchy than others.  
}
\value{
returns list with directory location and output file names
}
\author{
  Martha Ellis
}
\references{ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).}

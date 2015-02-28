\name{createReplicates}
\alias{create.landscapes}
\alias{createReplicates}
\title{Wrapper function to create replicate simulated landscapes}
\description{
Creates a series of files with replicated complete encounter histories based on habitat and population parameters.
}
\usage{
createReplicates(n_runs, map, Parameters, ... )
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
optional options (see details)
}
}
\details{
This function will produce \code{n_run} text files with simulated complete encounter histories.  Each individual file does not contain identifying information about the parameters used to create it, so it will be good practice to create a folder structure to keep separate simulations separate.  Default is to create a folder for the runs (./rSPACE_X) and name all encounter histories ./rSPACE_X/rSPACEx??.txt where ?? indicates the run number.  If \code{Parameters} is not specified, \code{\link{enter.parameters}} will be called.

\code{map} raster layer defining the landscape for the population simulation.  Values of \code{map} are used as weights for both selecting individual activity centers and scaling bivariate normal movement distributions.  Raster layer should use either longlat or utm projections with meter units.

There are additional arguments that can be included to specify other options, but some are more sketchy than others. 

Optional argument \code{skipConfirm} will bypass the user-confirmation for creating a folder structure.

Optional argument \code{run.label} changes the folder name for the run scenario. Defaults to './rSPACE_X'

Optional argument \code{base.name} can be used to change the base file name for encounter history replicate files. Defaults to 'rSPACEx'

Optional argument \code{filter.map} can be used to specify an alternative sampling frame for the landscape.  The default sampling frame consists of a rectangular grid built from the habitat information in \code{map} and user-specifications for \code{Paramters$grid_size} and \code{Parameters$sample.cutoff}.  By providing a secondary raster layer via \code{filter.map}, users can restrict sampling to specific areas (e.g., to exclude inaccessible private lands or restrict sampling to a single jurisdiction such as USFS, NPS, BLM, etc) or to provide a user-specific grid for the entire landscape.  The raster provided in \code{filter.map} will be expanded to match the extent of \code{map}, with non-coded regions coded 0. Raster values of 0 in \code{filter.map} are treated as areas to exclude from sampling.  If \code{filter.map} consists of 0/1 values only, the default rectangular grid will be restricted to areas with values of 1 in \code{filter.map}.  If \code{filter.map} contains additional values, non-zero values will be treated as indentifiers for sampling cells in a user-specified grid.  

Optional argument \code{printN} defaults to TRUE.  Logical - prints population size by run by year to a text file.

Optional argument \code{saveParameters} defaults to TRUE.  Logical - save rData file with parameter list used in run scenario.

Optional argument \code{saveGrid} defaults to TRUE.  Logical - save raster layer with sampling grid used in run scenario.
 
}
\value{
returns list with directory location and output file names
}
\author{
  Martha Ellis
}
\references{ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2014).}
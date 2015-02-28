\name{rSPACE-package}
\alias{rSPACE-package}
\alias{rSPACE}
\docType{package}
\title{Spatially-Explicit Power Analysis for Conservation and Ecology}
\description{
Conducts a spatially-explicit, simulation-based power analysis for detecting trends in population abundance through occupancy-based modeling. Applicable for evaluating monitoring designs in conservation and ecological settings.
}
\details{
\tabular{ll}{
Package: \tab rSPACE\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-02-27\cr
License: \tab GPL (>=2)\cr
}

There are many possible options for how a simulation and power analysis could be conducted.  We set up this package to provide options for how to conduct both the simulations and the data analysis.  We would like to help make this framework accessible, but we also acknowledge that it may not always work in situations for which the software was not originally designed.  We will be working to make more of the functionality accessible as we go along and will be interested to hear how things work (or fail to work) for anyone else who is interested in using our program.

For now, we have only made the basics of the package available.  We originally designed our population simulation for territorial carnivores.  Given a raster habitat layer, we distribute individuals according to a number of spacing rules (e.g. territory size, percent overlap, etc).  We use bivariate normal movement distributions, adjusted based on an underlying raster habitat layer, to build a probability of use layer for each individual.  We then combine the individual probabilities to create a layer describing the probability of use for at least one individual for the entire landscape.  The landscape is gridded into cells and used to create an encounter history at each possible grid cell.  

There are a large number of parameters needed to set up this stage of the analysis.  The parameters are stored as a list, which can either be entered manually by the user or via a dialog box (for details see \code{\link{enter.parameters}}).  The function, \code{\link{encounter.history}}, takes the list of parameters and a raster map of habitat on the landscape to produce a single, complete encounter history file, with options to error check the building process.  \code{\link{createReplicates}} will create replicated encounter history files.  These data will be stored as text files (.txt) that can either be used as inputs for Program MARK, etc or read back in later in the power analysis side of the simulations.

Just as there are many options for simulating a population on a landscape, there are many ways to potentially analyze each encounter history.  The way our analysis is set up, we have one wrapper function \code{\link{testReplicates}} set up to subset the encounter history files to simulate varying sampling effort.  \code{\link{testReplicates}} takes as an argument a function name that will define the analysis that is run on each potential encounter history.  This function should take the encounter history and some information about that file as arguments and return a data frame with the simulation results.  We have provided the original test file that we used to analyze simulated wolverine occupancy in Ellis et al. (2014) in \code{\link{wolverine_analysis}}.  The final output of the simulations at this stage is a text file storing simulation results, based on the output of \code{\link{wolverine_analysis}}.   

There are many options at each of these steps, but the basic process is:
\enumerate{
\item Enter parameters: \code{\link{enter.parameters}}
\item Build encounter history files: \code{\link{createReplicates}}
\item Analyze encounter histories: \code{\link{testReplicates}}
}
The help files will focus on these three steps, and explain additional details needed at each step.  
}
\author{
Martha Ellis, Jake Ivan, Jody Tucker, Mike Schwartz

Maintainer: Martha Ellis <martha.ellis@gmail.com>
}
\references{
ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).
}

\examples{       
## Not run
# data(WolverineHabitat)           # Loads example habitat map
# Parameters<-enter.parameters()   # Open dialog box with default parameter values
#
# # Create encounter history replicates
# createReplicates(n_runs=10, map=WolverineHabitat, Parameters=Parameters)
#
# # Analyze encounter histories (
# # Default: RDOccup model via RMark (Program MARK must be available))      
# testReplicates(folder="./rSPACE.X", Parameters)
}

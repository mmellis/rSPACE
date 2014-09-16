\name{rSPACE-package}
\alias{rSPACE-package}
\alias{rSPACE}
\docType{package}
\title{Spatially-explicit Power Analysis for Conservation and Ecology}
\description{
This package provides subroutines necessary to run a simple, spatially-explicit population simulation and then conduct a power analysis for detecting population trends through occupancy-based modeling.
}
\details{
\tabular{ll}{
Package: \tab rSPACE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-15\cr
License: \tab GPL (>=2)\cr
}

There are many possible options for how a simulation and power analysis could be conducted.  We set up this package to provide options for how to conduct both the simulations and the data analysis.  We would like to help make this framework accessible, but we also acknowledge that it may not always work in situations for which the software was not originally designed.  We will be working to make more of the functionality accessible as we go along and will be interested to hear how things work (or fail to work) for anyone else who is interested in using our program.

For now, we have only made the basics of the package available.  We originally designed our population simulation for territorial carnivores.  Given an underlying landscape, we distribute individuals according to a number of spacing rules (e.g. home range size, percent overlap, etc).  We use bivariate normal homeranges, adjusted based on the underlying habitat/landscape layer, to build a probability of use layer for each individual.  We then combine the individual probabilities to create a layer describing the probability of use for at least one individual for the entire landscape.  The landscape is then gridded into cells and used to create an encounter history at each possible grid cell.  

There are a large number of parameters needed to set up this stage of the analysis.  The parameters are stored as a list, which can either be entered separately or we have set up a dialog box to create it (for details see \code{\link{enter.parameters}}).  The function, \code{\link{encounter.history}}, takes the list of parameters and a raster map of the landscape to produce a single encounter history file.  Or you can use \code{\link{create.landscapes}} to create replicated encounter histories.  The encounter histories will be created as text files (.txt) that can either be used as inputs for Program MARK, etc or read back in later in the power analysis side of our simulations.

Just as there are many options for simulating a population on a landscape, there are many ways to potentially analyze each encounter history.  The way our analysis is set up, we have one wrapper function \code{\link{test_samples}} set up to subset the encounter history files to simulate varying sampling effort.  \code{\link{test_samples}} takes as an argument a function name that will define the analysis that is run on each potential encounter history.  This function should take the encounter history and some information about that file as arguments and return a data frame with the simulation results.  We have provided the original test file that we used to analyze simulated wolverine occupancy in Ellis et al. (2013) in \code{\link{wolverine_analysis}}.  The final output of the simulations at this stage is a text file storing simulation results, based on the output of \code{\link{wolverine_analysis}}.   

There are many options at each of these steps, but the basic process is:
\enumerate{
\item Enter parameters \code{\link{enter.parameters}}
\item Build encounter history files \code{\link{create.landscapes}}
\item Analyze encounter histories \code{\link{test_samples}}
}
The help files will focus on these three steps, and explain additional details needed at each step.  
}
\author{
Martha Ellis, Jake Ivan, Mike Schwartz

Maintainer: Martha Ellis <martha.ellis@gmail.com>
}
\references{
ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).
}

\examples{       
#Parameters<-enter.parameters()
  #create.landscapes(10, map=map, Parameters=Parameters)      # Will need to have a map to do this.
  #test_samples(main_wd,folder,Parameters,function_name="test_file")   ##Impossible!
}

\name{wolverine_analysis}
\alias{wolverine_analysis}
\title{
Default function for analyzing an individual encounter history
}
\description{
  Original function designed for analyzing wolverine encounter histories in Ellis et al. 2013. 
}
\usage{
wolverine_analysis(n_yrs, ch = NULL, n_visit = NULL, gap_yr = 0, FPC = 1, ... )
}

\arguments{
  \item{n_yrs}{
Number of years in the encounter history. Run with all else NULL to produce an initial data.frame for setting up further results.
}
  \item{ch}{
Character vector with encounters
}
  \item{n_visit}{
Number of visits per year
}
  \item{gap_yr}{
Code for possible alternate model formulations
}
  \item{FPC}{
Finite population correction
}
  \item{...}{
additional arguments
}
}
\details{
This function provides an example of an analysis function, which will be called by \code{\link{test_samples}} on each individual encounter history.  Should include appropriate model options for every subset of encounter histories.  Returns a dataframe with results from current simulation, which will be compiled into a single output file by \code{\link{test_samples}}.  Can include multiple model runs on the same encounter history by returning a data.frame with multiple lines.

The current function runs a robust design occupancy model through RMark with options for continuous or alternate year sampling.  Trend is tested using variance.components procedure.
}
\value{
data.frame with model results
}
\author{
Jake Ivan, Martha Ellis
}
\note{
The analysis function you give \code{test_samples} via \code{function_name} will get called in two places in \code{test_samples}.  First, it will get called by \code{RunAnalysis(n_yrs)} to set up a data.frame header.  Then it will get called in each step of the loop over all the encounter history files and subsets using \code{RunAnalysis(n_yrs, ch, n_visit, gap_yr,fpc , ... )}.  This is clunky, but for now, if you want to change \code{wolverine_analysis}, you need to make sure the function you specify will work with those two function calls.
}

\references{ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).}
\seealso{
\code{\link{test_samples}}
}
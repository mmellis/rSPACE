\name{getResults}
\alias{getResults}

\title{
  Plot power curves and return analysis data
}
\description{
Load data and plot power curves for representing each subset in a specified rSPACE folder. 
}
\usage{
getResults(folder, CI = 0.95, returnData = 1, plot=T)
}

\arguments{
  \item{folder}{
  Full path for rSPACE folder containing both simulation and analysis output
}
  \item{CI}{
  Significance value for one-tailed test
}
  \item{returnData}{
  Option to return contents of results file as a data.frame.   
}
  \item{plot}{
  Option to plot power curves [Logical, default=T]
}
}

\value{
If returnData=0, \code{getResults} returns nothing \cr
   returnData=1, the complete output for each analysis is returned \cr
   returnData=2, output of the analyses summarized by sampling intensity is returned.


For returnData=1, column headers include:
\itemize{
\item \code{p_est} Detection probability estimated from occupancy model                           
\item \code{trend} Trend parameter estimate from GLM using occupancy ests and vcv
\item \code{trendSE} Standard error of trend parameter                                  
\item \code{singular} Number of singularities from estimated occupancy model
\item \code{X1}-\code{X10} Occupancy estimates by year of simulation
\item \code{n_grid} Number of cells used for observed encounter history (subset)
\item \code{n_visits} Number of sampling occasions for subset
\item \code{detP} Per visit detection probability for subset
\item \code{alt_model} Alternative model specification for subset
\item \code{rn} Base name for original complete encounter history file used for subset
\item \code{count} Indicator variable for whether statistical criteria for detecting trend was met
}


For returnData=2, column headers include:
\itemize{
\item \code{n_grid} Number of cells used for observed encounter history (subset)
\item \code{n_visits} Number of sampling occasions for subset
\item \code{detP} Per visit detection probability for subset
\item \code{alt_model} Alternative model specification for subset
\item \code{total} Total number of replicates in which models ran successfully
\item \code{count} Number of replicates in which statistical criteria for detecting trend was met
\item \code{n_runs} Original number of replicates used
}

}
\author{
Martha Ellis
}

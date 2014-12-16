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
   returnData=1, the complete output for each analyses are returned \cr
   returnData=2, output of the analyses summarized by sampling intensity is returned.
}
\author{
Martha Ellis
}

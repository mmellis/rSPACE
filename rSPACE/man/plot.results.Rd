\name{plot.results}
\alias{plot.results}

\title{
  Plot power curves
}
\description{
Plot power curves for representing each subset in a specified rSPACE folder
}
\usage{
plot.results(folder, CI = 0.95, returnData = T)
}

\arguments{
  \item{folder}{
  Full path for rSPACE folder containing both simulation and analysis output
}
  \item{CI}{
  Significance value for two-tailed test
}
  \item{returnData}{
  Logical: return contents of results file as a dataframe?
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If returnData is TRUE, return contents of results file as a data.frame.  Otherwise, just creates figure.
}
\author{
Martha Ellis
}

\examples{
}
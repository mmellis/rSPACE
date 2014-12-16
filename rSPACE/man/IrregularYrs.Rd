\name{IrregularYrs}
\alias{IrregularYrs}
\docType{data}
\title{
Possible sampling frame for irregular visits across years
}
\description{
For use in resampling, a data.frame with 223 observations of visiting schedules over a 10 year period.  Each row represents a sampling cell, with 1/0 in each year representing whether the cell was visited or not. For use in mimicking the implementation of a sampling design where logistics, funding, or site access (weather, permits, etc) may interfer with planned regular sampling.
}
\usage{data(IrregularYrs)}
\format{
  A data frame with 223 observations of 10 variables.
  \describe{
    \item{\code{Y0}}{Indicator of whether cell was visited at start of study}
    \item{\code{Y1}}{Indicator of whether cell was visited in year 1}
    \item{\code{Y2}}{Indicator of whether cell was visited in year 2}
    \item{\code{Y3}}{Indicator of whether cell was visited in year 3}
    \item{\code{Y4}}{Indicator of whether cell was visited in year 4}
    \item{\code{Y5}}{Indicator of whether cell was visited in year 5}
    \item{\code{Y6}}{Indicator of whether cell was visited in year 6}
    \item{\code{Y7}}{Indicator of whether cell was visited in year 7}
    \item{\code{Y8}}{Indicator of whether cell was visited in year 8}
    \item{\code{Y9}}{Indicator of whether cell was visited in year 9}
  }
}
\details{
The IrregularYrs sampling option was developed to simulate sampling that does not follow a regular schedule (i.e. every year, or alternate year sampling) which is a reality of many monitoring studies.  Long term monitoring programs to assess trend often face logistic, environmental, or financial challenges that may unexpectedly cause deviation  from the planned sampling schedule (i.e. a heavy snow year precludes access to a portion of the sampling grid).  The IrregularYrs option can be used to investigate how such irregular sampling may affect  statistical power to detect trend.  This example data is based on systematic sampling conducted by U.S Forest Service Sierra Nevada carnivore monitoring program designed to detect declines in occupancy for the southern Sierra Nevada fisher population (described in detail in Zielinski et al. 2013). 
}
\source{
Data expanded from dataset analyzed in Zielinski, W. J., Baldwin, J. A., Truex, R. L., Tucker, J. M., and Flebbe, P. A. (2013). Estimating trend in occupancy for the southern Sierra fisher (Martes pennanti) population. Journal of Fish and Wildlife Management, 4(1), 3-19 
}
\examples{
data(IrregularYrs)
IrregularYrs[sample(nrow(IrregularYrs), 10, replace=TRUE), ]
}
\keyword{datasets}
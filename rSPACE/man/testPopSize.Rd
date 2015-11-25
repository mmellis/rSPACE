\name{testPopSize}
\alias{testPopSize}
\title{
Produce estimates for use probabilities based on population size
}
\description{
Tests a range of population sizes on a landscape and outputs use probabilities based on a single replication of the population for one visit, the maximum specified number of visits, and the asymptotic value.
}
\usage{
testPopSize(Nst, Nen, length.out = 10, map, pList)
}

\arguments{
  \item{Nst}{
  Smallest population size to test.  Must be larger than the number of types of individuals.
}
  \item{Nen}{
  Largest population size to test
}
  \item{length.out}{
  How many different population sizes to test
}
  \item{map}{
  Habitat map
}
  \item{pList}{
  Parameter list
}
}
\details{
  This function will attempt to fit a population onto the landscape for each population size.  It attempt to fit \code{length.out} possible values for the population size, equally spaced from \code{Nst} to \code{Nen}.  These values may not be realized in the population, if individuals cannot be fit onto the landscape according to buffer and habitat rules (i.e., the population can become saturated on the landscape before \code{Nen} is met).
}
\value{
Data.frame with realized population size (N), expected occupancy rate after 1 visit, expected occupancy rate after \code{Parameters$n_visits}, and asymptotic occupancy.
}
\author{
Martha Ellis
}
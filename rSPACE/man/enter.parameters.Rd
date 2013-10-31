\name{enter.parameters}
\alias{enter.parameters}
\title{
Enter parameters
}
\description{
Clicky box for help setting up a parameter list
}
\usage{
enter.parameters(Parameters = NULL)
}
\arguments{
  \item{Parameters}{
  Optional parameters to use as defaults in clicky box.
}
}
\details{
Creates a minimal set of parameters necessay to run simulations and analysis.  Additional arguments can be added to the outputted list as needed (see Examples).

Parameter definitions
\tabular{ll}{
Population parameters \tab ----------------\cr
Initial population size \tab Population size at the start of the simulation\cr
Population growth rate \tab Annual population multiplication rate (Lambda=N_t/N_0)\cr
Number of groups \tab Types of individuals to treat separately in the population simulation\cr
Proportion of population \tab Proportion of population in each group (e.g. sex ratio)\cr
Movement parameters \tab ----------------\cr
Buffer distance \tab Minimum distance between movement centers for each type of individual\cr
Minimum habitat value \tab Cutoff value for where individual centers should be located\cr
Movement radius \tab With percent of time in radius, used to define bivariate normal use distribution\cr
Percent of time in radius \tab With movement radius, used to define bivariate normal use distribution\cr
Truncation \tab Limits to long-distance movements\cr
Sampling parameters \tab ----------------\cr
Number of years \tab Maximum number of years for the simulation\cr
Number of visits \tab Maximum number of visits to each grid cell per year\cr
Cell size \tab Area of grid cell\cr
Percent of habitat \tab Cutoff for home much habitat must be included in cell to sample\cr
}
}
\value{
A parameter list with the following items
\describe{
\item{N}{
Initial population size
}
\item{lmda}{
  Population growth rate
}
\item{n_yrs}{Maximum number of years in the simulation
}
\item{n_visits}{Maximum number of visits per year
}
\item{grid_size}{Area of each grid cell
}
\item{HRcenter.cutoff}{Cutoff for value of \code{map} in \code{create.landscapes} where individual movement centers should be located.
}
\item{sample.cutoff}{Cutoff sum of \code{map} values required in a grid cell for that cell to be included.
}
\item{buffer}{Distance between individual movement centers for each type of individual
}
\item{howmuch}{Percent of time within movement radius for each type of individual
}
\item{howfar}{Movement radius for each type of individual
}
\item{trunk}{Truncation value (SDs from movement center to truncate)
}
\item{MFratio}{Proportion of population by type of individual
}
\item{detP}{Initial detection probability (By default, encounter histories are created with perfect detection (detP=1), then sampled with imperfect detection via \code{test_samples})
}
}
}
\author{
Martha Ellis, Mike Schwartz
}
\references{ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).}
\note{
Please direct all complaints about the damn clicky box to Mike Schwartz
}
\name{enter.parameters}
\alias{enter.parameters}
\title{
Enter parameters
}
\description{
Dialogue box for help setting up a parameter list
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
Creates a minimal set of parameters necessay to run simulations and analysis.  Additional arguments can be added to the outputted list as needed.

Parameter explanations
\tabular{ll}{
Population parameters \tab ----------------\cr
Initial population size \tab Population size at the start of the simulation\cr
Population growth rate \tab Annual population multiplication rate (Lambda=N_t/N_0)\cr
Number of years \tab Maximum number of years for the simulation\cr
Number of types \tab Types of individuals to treat separately in the population simulation\cr
Proportion of population \tab Proportion of population in each individual-type (e.g. sex ratio)\cr
Movement parameters \tab ----------------\cr
Buffer distance \tab Minimum distance between movement centers for each type of individual\cr
Movement radius \tab With \code{moveDistQ}, used to define bivariate normal use distribution\cr
Proportions of movements within radius \tab With movement radius, used to define bivariate normal use distribution\cr
Truncation \tab Limits to long-distance movements (proportion of movement distribution to allow)\cr
Sampling parameters \tab ----------------\cr
Cell size \tab Area of grid cell\cr
Minimum habitat value \tab Cutoff value for where individual centers should be located\cr
Minimum sampling value \tab Cutoff for home much habitat must be included in cell to sample\cr
Number of visits \tab Maximum number of visits to each grid cell per year\cr
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
\item{MFratio}{Proportion of population by type of individual
}
\item{buffer}{Distance between individual movement centers for each type of individual
}
\item{moveDist}{Percent of time within movement radius for each type of individual
}
\item{moveDistQ}{Movement radius for each type of individual
}
\item{maxDistQ}{Proportion of movement distribution to allow
}
\item{grid_size}{Area of each grid cell
}
\item{habitat.cutoff}{Minimum value of \code{map} where individual movement centers can be located.
}
\item{sample.cutoff}{Cutoff of \code{map} values above \code{habitat.cutoff} required for a grid cell to be included in sampling
}
\item{n_visits}{Maximum number of visits per year
}
\item{detP}{Initial detection probability = 1 (perfect detection to established availability)
}
}
}
\author{
Martha Ellis, Jody Tucker, Mike Schwartz
}
\note{
Please direct all complaints about the damn clicky box to Mike Schwartz
}
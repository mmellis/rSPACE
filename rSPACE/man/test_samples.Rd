\name{testReplicates}
\alias{testReplicates}
\alias{test_samples}
\title{
  Analyze a set of replicate encounter histories
}
\description{
For a folder containing simulated encounter histories: reads in complete replicated encounter histories, subsets data if applicable, and runs an analysis function on each subsetted encounter history.  Results are accumulated into a data.frame and written into a results output file.
}
\usage{
testReplicates(folder, Parameters, ... )
}
\arguments{
  \item{folder}{
Folder containing encounter history files (full path).  
}
  \item{Parameters}{
list of parameter values
}
  \item{...}{
additional arguments
}
}
\details{   
The main role of \code{testReplicates} is to read in complete encounter histories, subset them based on varying sampling effort, and then run an analysis on each "observed" encounter history (e.g. subsetted file).  The actual analysis that gets applied will be determined by the function specified with optional argument, \code{function_name}.  The default function, \code{\link{wolverine_analysis}}, includes code from the original analysis used in Ellis et al. (2013).

Options for subsetting the encounter histories can be specified via \code{Parameters}.   The number of visits per year to test can be specified by \code{Parameters$n_visit_test}, the proprotion of the grid included in the sample \code{Parameters$grid_sample}, the detection probability per visit \code{Parameters$detP_test}, and the sampling scheme by year \code{Parameters$alt_model}.  If these objects are not in the \code{Parameters} list, the default values are:
\itemize{
\item \code{Parameters$n_visit_test = 2:Parameters$n_visits}                           
\item \code{Parameters$detP_test    = c(1,0.8,0.2)}                                     
\item \code{Parameters$grid_sample  = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.75, 0.95)}
\item \code{Parameters$alt_model    = c(0,1)}
}
For \code{Parameters$alt_model}, codes will need to match analyses in the function specified by \code{function_name}.  Using \code{\link{wolverine_analysis}}, \code{Parameters$alt_model = 0} indicates continuous sampling and \code{Parameters$alt_model = 1} indicates alternate year sampling.

Additional arguements include
Optional argument \code{skipConfirm} will bypass the user-confirmation for creating a folder structure.  Defaults to FALSE

Optional argument \code{results.file} specifies a file name for storing analysis results.  Defauls to "sim_results.txt"

Optional argument \code{add} can be used when restarting an analysis, adding lines to a previously existing \code{results.file}.  Defaults to FALSE

Optional argument \code{base.name} must be specified if an alternative base file name was used in \code{createReplicates()}.  This can be specified as a single base for all files (e.g.,\code{base.name=rSPACEx"}), a vector of specific files (\code{base.name=c("rSPACEx1.txt", "rSPACEx2.txt",rSPACEx3.txt")}) or as a regular expression (\code{base.name="rSPACEx|rSPACEy"}).

Optional argument \code{FPC} indicates whether to use a finite population correction.  \code{FPC} defaults to TRUE; however, this option requires that a finite population correction is available in the analysis function.
}

Optional argument \code{sample_matrix} allows for irregular gaps in sampling histories.  Input should be a matrix with one column for each year of sampling, each entry should be a 1/0 reflected that a cell is visited/not visited in a given year.  Analysis function must be able to handle missing values.  Visitation history for stations will be randomly selected from rows of \code{sample_matrix}.  See \code{data(IrregularYrs)} for example.

\value{
 Stores simulation results based on output specified in the analysis to \code{sim_results.txt}.  Returns time used for the function run.
 
 By default, variables outputted to \code{sim_results.txt} include:
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
}
}

\author{
Martha Ellis, Jake Ivan
}
\references{ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).}
\seealso{
\code{\link{wolverine_analysis}}
}
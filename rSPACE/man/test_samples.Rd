\name{test_samples}
\alias{test_samples}
\title{
  Analyze a set of encounter histories
}
\description{
For a folder containing simulated encounter histories: reads in complete encounter histories, subsets data if applicable, and runs an analysis function on each subsetted encounter history.  Results are accumulated into a data.frame and written into a results output file.
}
\usage{
test_samples(folder, Parameters, ... )
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
The main role of \code{test_samples} is to read in complete encounter histories, subset them based on varying sampling effort, and then run an analysis on each "observed" encounter history (e.g. subsetted file).  The actual analysis that gets applied will be determined by the function specified with optional argument, \code{function_name}.  The default function, \code{\link{wolverine_analysis}}, includes code from the original analysis used in Ellis et al. (2013).

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

Optional argument \code{base.name} must be specified if an alternative base file name was used in \code{create.landscapes()}

Optional argument \code{FPC} indicates whether to use a finite population correction.  \code{FPC} defaults to TRUE; however, this option requires that a finite population correction is available in the analysis function.
}

\value{
 Stores simulation results based on output specified in the analysis to \code{sim_results.txt}.  Returns time used for the function run.
}
\author{
Martha Ellis, Jake Ivan
}
\references{ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).}
\seealso{
\code{\link{wolverine_analysis}}
}
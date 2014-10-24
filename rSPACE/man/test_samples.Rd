\name{test_samples}
\alias{test_samples}
\title{
  Analyze a set of encounter histories
}
\description{
For a folder containing simulated encounter histories, reads in encounter histories, subsets data if applicable, and runs a function on each encounter history.  Results are accumulated into a data.frame and written into a results output file.
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
The main role of \code{test_samples} is to read in encounter histories, subset them based on varying sampling effort, and then run an analysis on each "observed" encounter history (e.g. subsetted file).  The actual analysis that gets applied will be determined by the function specified with \code{function_name}.  The default function, \code{\link{wolverine_analysis}}, includes code from the original analysis used in Ellis et al. (2013).

Options for subsetting the encounter histories can be specified via \code{Parameters}.   For example, the number of visits per year to test can be specified by \code{Parameters$n_visit_test}, the percent of the grid included in the sample \code{Parameters$grid_sample}, the detection probability per visit \code{Parameters$detP_test}, and the sampling scheme by year \code{Parameters$sample_yrs}.  If these objects are not in the \code{Parameters} list, the default values are:
\itemize{
\item \code{Parameters$n_visit_test = 2:Parameters$n_visits}                           
\item \code{Parameters$detP_test    = c(1,0.8,0.2)}                                     
\item \code{Parameters$grid_sample  = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.75, 0.95)}
\item \code{Parameters$sample_yrs   = c(0,1)}
}
For \code{Parameters$sample_yrs}, codes will need to match analyses in the function specified by \code{function_name}.  Using \code{\link{wolverine_analysis}}, \code{Parameters$sample_yrs = 0} indicates continuous sampling and \code{Parameters$sample_yrs = 1} indicates alternate year sampling.
}

\value{
 Stores simulation results based on output specified in the analysis to \code{sim_results.txt}.  Any messages returned from the analysis function will be displayed on the console.  Returns time used for the function run.
}
\author{
Martha Ellis, Jake Ivan
}
\references{ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).}
\seealso{
\code{\link{wolverine_analysis}}
}
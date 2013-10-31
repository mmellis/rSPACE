\name{test_samples}
\alias{test_samples}
\title{
  Analyze a set of encounter histories
}
\description{
For a folder containing simulated encounter histories, reads in encounter histories, subsets data if applicable, and runs a function on each encounter history.  Results are accumulated into a data.frame and written into a results output file.
}
\usage{
test_samples(main_wd, folder, Parameters, function_name = "test_file", SubPop = NULL, sample_matrix = NULL, base.name = "tmp", n_runs = NULL)
}
\arguments{
  \item{main_wd}{
Base working directory, for simplicity if running multiple tests for a single species
}
  \item{folder}{
Folder containing encounter history files.  Should be be set as current working directory
}
  \item{Parameters}{
list of parameter values
}
  \item{function_name}{
\code{character()} name of function that takes an individual encounter history as input, runs analysis, and outputs results as data.frame
}
  \item{SubPop}{
  Optional raster layer including a subset of the population from which to sample  
}
  \item{sample_matrix}{
  Optional matrix matching encounter histories matrix with cells to include in sample
}
  \item{base.name}{
base name for encounter history files
}
  \item{n_runs}{
number of encounter histories to include in analysis.  By default all files with \code{base.name} in \code{folder} will be included.
}
}
\details{   
The main role of \code{test_samples} is to read in encounter histories, subset them based on varying sampling effort, and then run an analysis on each subset of file.  The actual analysis that gets applied will be determined by whatever function you name under \code{function_name}.  We've provided the original analysis we used in Ellis et al. (2013), \code{\link{test_file}}.  The analysis function should be able to handle every subset you create in \code{test_samples}.

The options for subsetting the encounter histories include changing the number of visits per year \code{Parameters$n_visit_test}, the percent of the grid included in the sample \code{Parameters$grid_sample}, the detection probability per visit \code{Parameters$detP_test}, and the sampling scheme by year \code{Parameters$sample_yrs}.  If these objects are not in the \code{Parameters} list, the default values are:
\itemize{
\item \code{Parameters$n_visit_test = 2:Parameters$n_visits}                           
\item \code{Parameters$detP_test    = c(1,0.8,0.2)}                                     
\item \code{Parameters$grid_sample  = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.75, 0.95)}
\item \code{Parameters$sample_yrs   = c(0,1)}
}
For \code{Parameters$sample_yrs}, the different codes allow for different analyses to be run in \code{\link{test_file}}.  For the code we provide, \code{Parameters$sample_yrs = 0} indicates continuous sampling and \code{Parameters$sample_yrs = 1} indicates alternate year sampling.

Changing the value of \code{function_name} allows you to run whatever analysis you choose on the subsetted encounter histories.
}

\value{
 Stores results data.frame to file named \code{sim_results1.txt}.  The only output to the console is the time used for the function run.
}
\author{
Martha Ellis, Jake Ivan
}
\references{ELLIS, MARTHA M., JACOB S. IVAN, and MICHAEL K. SCHWARTZ. "Spatially Explicit Power Analyses for Occupancy-Based Monitoring of Wolverine in the US Rocky Mountains." Conservation Biology (2013).}
\seealso{
\code{\link{test_file}}
}
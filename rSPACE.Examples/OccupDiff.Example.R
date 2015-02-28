################################################################################
# Script to run rSPACE as a comparison of two occupancy estimates              #
#   M. Ellis, 2015/02/20, rSPACE v1.0.7                                        #
################################################################################
getwd()

### Part 1: Load habitat layer, enter parameters, and create landscapes as usual
library(raster)
library(rSPACE)

WolverineHabitat<-raster::raster(system.file('external/WolvHabitat_Bitterroot.tif', package='rSPACE'))

BaseParameters<-list(  # Base wolverine parameters
  N               = 50,
  MFratio         = c(0.6, 0.4),
  buffer          = c( 16, 25),
  moveDist        = c(8.5, 12.5),
  moveDistQ       = c(0.9, 0.7),
  maxDistQ        = c(.95, .95),
  grid_size       = 100,
  habitat.cutoff  = 1,
  sample.cutoff   = 0.5,
  n_visits        = 6)

  # Option 1: Effect size is still in terms of lmbda, especially if interested
  #  in the number of years required to detect a population change.
  #  Alternatively, you'd need to rig it so that the combination of n_yrs and
  #  lmda produce the difference in abundance that you would like to detect.
  BaseParameters<-c(BaseParameters,
    lmda  = 0.933,
    n_yrs = 10)

#  # Option 2: Effect size is in terms of change in abundance.  Produce change
#  #  in the first year, only need two years of simulation, but less interesting
#  #  in terms of time to effect, potentially.
#  #  Need to calculate lmda require to produce change in N in one time step.
#  #  E.g. to detect a 10% decline between samples, Nt/N0 = lmda = 0.9
#  BaseParameters<-c(BaseParameters,
#    lmda  = 0.9,
#    n_yrs = 2)

Example<-encounter.history(map=WolverineHabitat, Parameters=BaseParameters, showSteps=T)

createReplicates(n_runs = 10, ## Usually good to try running on a small example first
   map = WolverineHabitat, Parameters = BaseParameters,
   run.label = 'rSPACE.Example', skipConfirm=T)


### Part 2: Replace analysis function in testReplicates

source('OccupDiff.R') # Code for analysis function, based on Guillera-Arroita & Lahoz-Monfort 2012 Appendix S3
                      # (Need to pick up additional file from AnalysisFunctions for this script to work)

# Add subsetting options
BaseParameters<-c(BaseParameters,
  list(
    n_visit_test=c(3:6),                  # Number of sampling occasions per year
    detP_test   =c(0.8),                  # Per visit detection probability (if spp present)
    grid_sample =c(0.05,0.15,0.25,0.35,   # Percent of grid to sample
                  0.45,0.55,0.75,0.95),
    alt_model   =c(2,5,10)                # For  OccupDiff, alt_model specifies
  ))                                      #  which year to use for the second sample
                                          #  Here, I'm picking years 2, 5, and 10

testReplicates('rSPACE.Example', BaseParameters, function_name='OccupDiff')

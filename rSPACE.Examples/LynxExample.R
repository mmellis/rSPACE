################################################################################
#     Power analysis for Lynx in Colorado - rSPACE Example                     #
#             Prepared 11/21/2014, M. Ellis, J. Ivan                           #
#                                                                              #
# Additional options demonstrated:                                             #
#  - Varying population growth rate in time                                    #
#  - Editing analysis function                                                 #
################################################################################

library(rSPACE)
library(raster)
library(ggplot2)

####
# Check working directory (Should be the base directory encompassing all scenarios)
getwd()

#### Build landscapes

# Habitat suitability layer (borrowed from Wolverine example)
  Habitat<-raster::raster(system.file('external/WolvHabitat_Bitterroot.tif', package='rSPACE'))

# Input parameter list
  BaseParameters<-list(
    N               = 250,                # Initial population size
    lmda            = 0.933,              # Population growth rate
    n_yrs           = 10,                 # Maximum number of years in simulation
    n_visits        = 6,                  # Maximum number of visits per year
    grid_size       = 75,                 # Cell size in grid
    MFratio         = c(0.6, 0.4),        # Ratio of types of individuals
    buffer          = 2*c(4.9, 5.7),      # Distance between individual center locations
    howfar          = c(4.9, 5.7),        # Movement radius
    howmuch         = c(0.95, 0.95),      # Proportion of time in radius
    trunk           = c(  2, 2),          # Truncate movements above 1 SD
    HRcenter.cutoff = 0.5,                # Minimum habitat value required for individual center locations
    sample.cutoff   = 0.75                # Proportion of pixels above minimum habitat value to include cell
    )

# Check steps for building encounter histories
  Example<-encounter.history(map=Habitat, Parameters=BaseParameters, showSteps=T)

# Create replicate landscapes
  create.landscapes(n_runs=100, map=Habitat, Parameters=BaseParameters, run.label='BaseScenario')

#### Analyze landscapes

## Analysis function to be applied to each observed encounter history ##################################
lynx_analysis<-function(n_yrs, ch=NULL, n_visit, M, FPC, ...){                                         #
  TrueModel=list(...)$TrueModel  #Passes argument from test_samples additional argument list           #
  DF<-data.frame(count=NA, Model=1:4, p_est=NA, AICc=NA, lnl=NA, PSI=matrix(NA,4,n_yrs))               #
  if(is.null(ch)) return(DF)                                                                           #
                                                                                                       #
  tryN<-suppressMessages(tryCatch(x, error=function(e) return(NA)))                                    #
                                                                                                       #
  mark_data<-data.frame(ch=ch,freq=rep(1,length(ch)),stringsAsFactors=F)                               #
  test_processed=process.data(mark_data,model="RDOccupPE",time.intervals=time_int(n_visit,n_yrs))      #
  test_ddl=make.design.data(test_processed)                                                            #
    test_ddl$Psi$Trend=c(1:n_yrs)                                                                      #
    test_ddl$Psi$Trend2=test_ddl$Psi$Trend^2                                                           #
    test_ddl$Psi$Trend3=test_ddl$Psi$Trend^3                                                           #
                                                                                                       #
  p.session=list(formula=~session)                                                                     #
  epsilon.t=list(formula=~time)                                                                        #
  psi.options= list(                                                                                   #
    psi.dot=list(formula=~1),                                                                          #
    psi.T  =list(formula=~Trend),                                                                      #
    psi.T2 =list(formula=~Trend+Trend2),                                                               #
    psi.T3 =list(formula=~Trend+Trend2+Trend3))                                                        #
                                                                                                       #
  Models<-lapply(1:length(psi.options), function(opt){                                                 #
              tryN(mark(test_processed, test_ddl,                                                      #
                model.parameters=list(Psi=psi.options[[opt]], Epsilon=epsilon.t, p=p.session),         #
                delete=T,output=F,silent=T))})                                                         #
                                                                                                       #
   DF$AICc <-sapply(Models, function(x) ifelse(is.null(x$results),NA,x$results$AICc)                   #
   DF$lnl  <-sapply(Models, function(x) ifelse(is.null(x$results),NA,x$results$lnl)                    #
   DF$p_est<-sapply(Models, function(x) ifelse(is.null(x$results),NA,                                  #
                x$results$real$estimate[which(row.names(x$results$real)=="p g1 s1 t1")])               #
   DF$PSI  <-sapply(Models, function(x) ifelse(is.null(x$results),NA,                                  #
                x$results$real$estimate[grepl('Psi',row.names(x$results$real))]                        #
                                                                                                       #
   DF$count[TrueModel]<-as.numeric(which.min(DF$AICc) == TrueModel)                                    #
        # DF$count returns 1/0 for the true model indicating if it                                     #
        #    would be the best model according to the data                                             #
        #   NAs returned for all other models                                                          #
                                                                                                       #
  return(DF)                                                                                           #
  }                                                                                                    #
                                                                                                       #
               ### End of code for specifying analysis function ########################################

 SubsettingParameters<-list(
      n_visit_test=c(3:6),                  # Number of sampling occasions per year
      detP_test   =c(0.8),                  # Per visit detection probability (if spp present)
      grid_sample =c(0.05,0.15,0.25,0.35,   # Percent of grid to sample
                    0.45,0.55,0.75,0.95),
      sample_yrs  =c(0)                     # Alternative models specifications:
    )                                       #  0 = continuous sampling only


 test_samples('./BaseScenario', Parameters=c(BaseParameters,SubsettingParameters), 
                function_name="lynx_analysis", TrueModel=1)   
                #TrueModel is an argument specific to lynx_analysis()






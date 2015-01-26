################################################################################
#     Power analysis for Sierra Neveda Fisher - rSPACE Example                 #
#             Prepared 11/13/2014, M. Ellis, J. Tucker                         #
#                                                                              #
# Additional options demonstrated:                                             #
#  - Secondary filter map                                                      #
#  - Effective sample area                                                     #
#  - Alternative model specification (alternate or irregular year sampling)    #
################################################################################

library(rSPACE)
library(raster)
library(ggplot2)



#### Build landscapes

# Habitat suitability layer (borrowed from Wolverine example)
  Habitat<-raster::raster(system.file('external/WolvHabitat_Bitterroot.tif', package='rSPACE'))

# Secondary filter map to restrict sampling area
  SampleFrame<-raster::raster(system.file('external/ExampleSampleFrame.tif', package='rSPACE'))
    table(getValues(SampleFrame)) # 1/0 values indicate areas to include/exclude in sampling

# Input parameter list
  BaseParameters<-list(
    N               = 250,                # Initial population size
    lmda            = 0.933,              # Population growth rate
    n_yrs           = 10,                 # Maximum number of years in simulation
    n_visits        = 5,                  # Maximum number of visits per year
    grid_size       = 25,                 # Cell size in grid
    MFratio         = c(0.6, 0.4),        # Ratio of types of individuals
    buffer          = c(3.75, 6.30),      # Distance between individual center locations
    moveDist        = c(2.5, 4.2),        # Movement radius
    moveDistQ       = c(0.9, 0.9),        # Proportion of time in radius
    maxDistQ        = c(0.68,0.68),       # Truncate movements above 1 SD
    habitat.cutoff  = 0.5,                # Minimum habitat value required for individual center locations
    sample.cutoff   = 0.5,                # Proportion of pixels above minimum habitat value to include cell
    filter.cutoff   = 0.95                # Proportion of pixels within SampleFrame to include cell
    )                

# Comparison parameters - What if sampling only effectively covers 10km2 of the 25km2 sample cell?
  TestParameters<-c(BaseParameters, list(Effective.sample.area=10))
 
# Check steps for building encounter histories
  Example1<-encounter.history(map=Habitat, Parameters=BaseParameters, filter.map=SampleFrame, showSteps=T)
  Example2<-encounter.history(map=Habitat, Parameters=TestParameters, filter.map=SampleFrame, showSteps=T)

# Create replicate landscapes
  create.landscapes(n_runs=100, map=Habitat, Parameters=BaseParameters, filter.map=SampleFrame, run.label='Base')
  create.landscapes(n_runs=100, map=Habitat, Parameters=TestParameters, filter.map=SampleFrame, run.label='Test')

 
 
 
  
#### Analyze landscapes

 SubsettingParameters<-list(
      n_visit_test=c(4),                    # Number of sampling occasions per year
      detP_test   =c(0.8),                  # Per visit detection probability (if spp present)
      grid_sample =c(0.05,0.15,0.25,0.35,   # Percent of grid to sample
                    0.45,0.55,0.75,0.95),
      alt_model   =c(0,1,2)                 # Alternative models specifications:
    )                                       #  0 = continuous sampling 
                                            #  1 = alternate year sampling
                                            #  2 = irregular year sampling


 data(IrregularYrs) # Possible set of irregular sampling histories over time  ##!!!!!
 
 test_samples('./Base', Parameters=c(BaseParameters,SubsettingParameters), sample_matrix=IrregularYrs, skipConfirm=T)
 test_samples('./Test', Parameters=c(TestParameters,SubsettingParameters), sample_matrix=IrregularYrs, skipConfirm=T)





#### Compare results

  BaseResults<-getResults('./Base/output', CI=0.95, returnData=2, plot=F)
  TestResults<-getResults('./Test/output', CI=0.95, returnData=2, plot=F)
  
  Results<-rbind(data.frame(Run='Base',BaseResults), data.frame(Run='Test',TestResults))
  
  ggplot(subset(Results, detP<1),  
            aes(x=n_grid, y=(count/n_runs), colour=factor(gap_yr), linetype=factor(detP)))+
          geom_line(aes(group=interaction(n_visits, detP, gap_yr)), size=1.25)+
          scale_colour_discrete(name = "Model Specification",guide="legend",
            labels=c('Continuous sampling', 'Alternate-year sampling','Irregular-year sampling')) +
          scale_linetype_discrete(name=expression(p["sim"]))+
          scale_y_continuous(limits=c(0,1))+
          labs(x="Number of cells sampled", 
               y="Detected trend/Number of runs")+
          facet_grid(.~Run, labeller=label_both) 
    
      
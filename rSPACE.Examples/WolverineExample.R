################################################################################
#     Script to demonstrate power analysis for Wolverine - MEE manuscript      #
#             Started 9/16/2014, M. Ellis                                      #
################################################################################

## 1.) Check working directory (Should be the base directory encompassing all scenarios)
getwd()


## 2.) Load package
  library(rSPACE)

  #-- Package overview and guidance
  help(rSPACE, package='rSPACE')



## 3.) Load habitat layer  
  #data(WolverineHabitat) ### This doesn't work anymore?! Stupid! Awful!
  WolverineHabitat<-raster::raster(system.file('external/WolvHabitat_Bitterroot.tif', package='rSPACE'))



## 4.) Input parameter list 
  
  #-- Manual input version
  BaseParameters<-list(
    N               = 50,                 # Initial population size
    lmda            = 0.933,              # Population growth rate
    n_yrs           = 10,                 # Maximum number of years in simulation
    MFratio         = c(0.6, 0.4),        # Ratio of types of individuals                 
    buffer          = c( 16, 25),         # Distance between individual center locations  
    moveDist        = c(8.5, 12.5),       # Movement radius                               
    moveDistQ       = c(0.9, 0.7),        # Proportion of movements in radius                     
    maxDistQ        = c(.95, .95),        # Truncate movements beyond 95% quantile                 
    grid_size       = 100,                # Cell size in grid                                                   
    habitat.cutoff  = 1,                  # Minimum habitat value required for individual center locations
    sample.cutoff   = 0.5,                # % pixels in cell above habitat value to include cell in sample frame  
    n_visits        = 6)                  # Maximum number of visits per year  
                                          
  #-- Alternative dialogue input for parameters
  # BaseParameters <- enter.parameters()   
  
  #-- Optional: specify subsetting levels for later use.  For this example, we 
  #    opted change the defaults to make it run faster (reduced number of 
  #    subsets to test).
  BaseParameters<-c(BaseParameters,
    list(
      n_visit_test=c(3:6),                # Number of sampling occasions per year
      detP_test   =c(0.8),                  # Per visit detection probability (if spp present)
      grid_sample =c(0.05,0.15,0.25,0.35,   # Percent of grid to sample
                    0.45,0.55,0.75,0.95),
      alt_model   =0                        # Possible alternative models specifications
    ))


## 5.) Optional: Check encounter history steps
Example1<-encounter.history(map=WolverineHabitat, Parameters=BaseParameters, showSteps=T)


## 6.) Run each scenario
  
  ## Scenario 1: N = 50, lmda = 0.933
  #- 6a.) Set parameters for scenario
    Plist_Scenario1<-BaseParameters
      Plist_Scenario1$N    <- 50
      Plist_Scenario1$lmda <- 0.933
  
  #- 6b.) Create simulated encounter history replicates 
    create.landscapes(n_runs=10, map=WolverineHabitat, Parameters=Plist_Scenario1,
                    run.label='rSPACE.Scenario1')
         # Note: 100 runs will take a long time to run both create.landscapes and test_samples.
         #        Reducing the number of runs is good for examples and testing 
         #        because it will run more quickly but will cause figures to 
         #        jump around more.
                    
  #- 6c.) Analyze subsetted encounter history replicates
    test_samples('./rSPACE.Scenario1', Parameters=Plist_Scenario1, skipConfirm=T)
   
  #- 6d.) Check out the results!
    #-- Power curves    
    getResults('./rSPACE.Scenario1', CI=0.9, returnData=0)

    #-- Power of monitoring when 20 cells are sampled
    findPower('./rSPACE.Scenario1', CI=0.9, n_grid=20)
    
    #-- Number of sampled cells needed to have 80% power
    findPower('./rSPACE.Scenario1', CI=0.9, pwr=0.8)
  
  #- 6e.) Rinse, lather, repeat for the other scenarios  
                    
  ## Scenario 2: N = 25, lmda = 0.933
  Plist_Scenario2<-BaseParameters
    Plist_Scenario2$N    <- 25
    Plist_Scenario2$lmda <- 0.933
    
  create.landscapes(n_runs=100, map=WolverineHabitat, Parameters=Plist_Scenario2,
                    run.label='rSPACE.Scenario2', skipConfirm=T) 
                                       
  test_samples('./rSPACE.Scenario2', Parameters=Plist_Scenario2, skipConfirm=T)
  
  
  ## Scenario 3: N = 50, lmda = 1.041
  Plist_Scenario3<-BaseParameters
    Plist_Scenario3$N    <- 50
    Plist_Scenario3$lmda <- 1.041
  
  create.landscapes(n_runs=100, map=WolverineHabitat, Parameters=Plist_Scenario3,
                    run.label='rSPACE.Scenario3', skipConfirm=T)
                    
  test_samples('./rSPACE.Scenario3', Parameters=Plist_Scenario3, skipConfirm=T)
  
  
  ## Scenario 4: N = 25, lmda = 1.041
  Plist_Scenario4<-BaseParameters
    Plist_Scenario4$N    <- 25
    Plist_Scenario4$lmda <- 1.041
   
  create.landscapes(n_runs=100, map=WolverineHabitat, Parameters=Plist_Scenario4,
                    run.label='rSPACE.Scenario4', skipConfirm=T)  
                     
  test_samples('./rSPACE.Scenario4', Parameters=Plist_Scenario4, skipConfirm=T)   
 


                    
## 7.) Combine results to create a fancy figure

  #-- Set up folder locations and parameters  
  ResultsList<-list()
  N   <- c(50, 25, 50, 25)             # Population size for each scenario
  lmda<- c(0.933, 0.933, 1.041, 1.041) # Population growth rate for each scenario
  
  folders<-c('./rSPACE.Scenario1',     # Folder location for each scenario
             './rSPACE.Scenario2',
             './rSPACE.Scenario3',
             './rSPACE.Scenario4')
  
  #-- Loop through folders
  for(i in 1:length(folders))          
  {
    # Read in results file from each folder
    ResultsList[[i]]<-getResults(folders[i], CI=0.9, returnData=2, plot=F) 
    
    # Add identifying info for each scenario
    ResultsList[[i]]<-data.frame(ResultsList[[i]], N=N[i], lmda=lmda[i])
  }
  
  #-- Combine results into one data.frame
  ResultsList<-rbind.fill(ResultsList) 
  
  #-- Facetted plot of scenarios for comparison
  facet_labels<-function(var, value){
    if(var=='lmda') value<-lapply(value, function(v) bquote(lambda==.(v)))
    if(var=='N')    value<-lapply(value, function(v) bquote(N==.(v)))
    return(value)}
  
  ggplot(subset(ResultsList, detP<1),  
            aes(x=n_grid, y=(count/n_runs), colour=n_visits, linetype=factor(detP)))+
          geom_line(aes(group=interaction(n_visits, detP)), size=1.25)+
          scale_colour_gradient(name = "# visits",guide="legend", breaks=c(3,4,5,6)) +
          scale_linetype_discrete(name=expression(p["sim"]))+
          scale_y_continuous(limits=c(0,1))+
          labs(x="Number of cells sampled", 
               y="Detected trend/Number of runs")+
          facet_grid(N~lmda, labeller=facet_labels) 
                   
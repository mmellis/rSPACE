################################################################################
#     Covariate power analysis using rSPACE                                    #
#       M. Ellis, 3/5/2015 (in progress)                                       #
#       Will require development version, rSPACE v1.1.1                        #
################################################################################

# Premise: suppose habitat suitability layer is built on other underlying
#   covariate values.  We want to test the ability of occupancy models to model
#   covariate relationships.

library(rSPACE)
library(raster)

## 1. Build fake data set to simulate from
# First, we need a fake map...
FakeMap<-local({ #Code suggested by reviewer from MEE manuscript...THANKS!
  temp<-list()
  temp$x<-seq(0, 100000, by=1000)
  temp$y<-seq(0, 100000, by=1000)
  temp$z<-matrix(1, length(temp$x), length(temp$y))
  raster(temp, crs=CRS("+proj=utm +zone=12 +ellps=GRS80 +datum=NAD83 +units=m"))
  })
  
# Add fake covariates...       
# To make things easier, grid the landscape and create covariates on grid instead
#  of trying to make a spatial layer with nicely varying covariates.  In the end,
#  for analysis, covariates need to be summarized by cell anyway.  This simplifies.
#  No autocorrelation.

pList<-list(grid_size=100, sample.cutoff=0, habitat.cutoff=0) #100km2 grid cells
grd<-rSPACE:::createGrid(FakeMap, pList)

CovariateData<-data.frame(grdID=unique(grd[grd!=0]))
  CovariateData$A<-rnorm(nrow(CovariateData))
  CovariateData$B<-rnorm(nrow(CovariateData))
  CovariateData$Y<- 0.8 + 0.2*CovariateData$A - 0.5*CovariateData$B
  CovariateData$HSI<-(CovariateData$Y-min(CovariateData$Y))/(max(CovariateData$Y)-min(CovariateData$Y))
  
HabitatMap<-setValues(FakeMap, values=c(0,CovariateData$HSI)[grd+1])

## 2. Run population simulation as normal
pList<-c(pList, list(
    N               = 50,                 
    #lmda            = exp(rnorm(9,sd=0.1)), # Avg no trend over study period             
    n_yrs           = 10,                    
    MFratio         = c(0.6, 0.4),        
    buffer          = c( 16, 25),         
    moveDist        = c(8.5, 12.5),       
    moveDistQ       = c(0.9, 0.7),        
    maxDistQ        = c(.95, .95),        
    n_visits        = 6)) 


# Replacing lmda=exp(rnorm(pList$n_yrs-1, sd=0.1)) for each replicate adds random 
#  yearly fluctation around a constant trend that will be unique for each
#  replicate.  Important since rSPACE uses pretty extreme site fidelity each year.
Example<-encounter.history(map=HabitatMap, 
          Parameters=c(pList,list(lmda=exp(rnorm(pList$n_yrs-1, sd=0.1)))), 
          showSteps=T)

createReplicates(n_runs=10, 
  map=HabitatMap, 
  Parameters=c(pList,list(lmda=exp(rnorm(pList$n_yrs-1, sd=0.1)))),   
  run.label='rSPACE.Example')



## 3. Change analysis function to run covariate analysis
source('CovariateAnalysis.R') # Must get new function code from source

pList<-c(pList, list(
      lmda=1,
      n_visit_test=c(3:6),                  
      detP_test   =c(0.8),                  
      grid_sample =c(0.05,0.15,0.25,0.35,   
                    0.45,0.55,0.75,0.95),
      alt_model   =0                        
    )) 

testReplicates('rSPACE.Example', pList, function_name='CovariateAnalysis',
                  CovDF = CovariateData[,c(1:3)], # Passed to CovariateAnalysis() via ...
                                                  #   Ids in col 1, covariates ONLY in additional cols
                  skipConfirm=T)
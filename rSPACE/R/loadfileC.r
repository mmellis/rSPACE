## links the functions from C++
## UPDATED 3/1/13 to match loadfileC_fisher & loadfileC_cyclic
#dyn.load("SPACE.dll")

##check function names from objdump  (In CMD: objdump -x fisher.dll >myfile.txt)
#	[   0] _Z10reduce_popPiS_PdS_S_S0_S0_S0_S_S0_S0_S0_S0_S0_
#	[   1] _Z10sample_indPdS_PiS_S0_S0_S0_
#	[   2] _Z11filter_gridPiPdS0_S_S0_
#	[   3] _Z11rdist_earthdddd
#	[   4] _Z11use_surfacePdS_PiS_S_S_S0_S_S_S_
#	[   5] _Z7wrapperPdS_S_PiS0_S_S_S_S0_S_S_S0_S0_S0_S_S_S_S_
#	[   6] _Z9calc_probPdPiS0_S_S0_S0_S_
#	[   7] _Z9make_gridPdS_S_PiS0_

sampling_fxn = 'sample_ind'
use_fxn = 'use_surface'
wrapper_fxn =  'wrapper'
calc_prob = 'calc_prob'
make_grid = 'make_grid'
filter_grid = 'filter_grid'
rd_pop = 'reduce_pop'

is.loaded(c(sampling_fxn, use_fxn, wrapper_fxn, calc_prob, make_grid))

smple <- function(map, use, N, buffer, wght=F){      #functions for individual calls to C++
   n = length(use)
   if(wght==T)
    {  smp =  sample(n,n, prob=getValues(map)[use])-1  # Changed 8/8/2012: needs to be -1 to avoid bad reference in C++
    }else{ smp=  sample(n,n)-1 }
   x = coordinates(map)[use,1]
   y = coordinates(map)[use,2]
   use10=rep(0,n)
    use10[smp[1]+1]= as.numeric(grepl('+proj=longlat',proj4string(map)))  
   
   USE = .C(sampling_fxn,
                  as.double(x),               #Longitude
                  as.double(y),               #Latitude
                  as.integer(N),              #Desired number of wolverines to place
                  as.double(buffer),          #Minimum distance between points
                  use = as.integer(use10),    #Indicator of whether to include or not
                  as.integer(n),              #Number of snow_points to check (maxid)
                  as.integer(smp)            #Random order to use
                )$use
   return((1:n)[USE==1]) }
   
use_surface<-function(Wolv, howmuch, howfar, map, trunk = trunk){
  sd_xy<-var_solver(howmuch, howfar, map)
  USE = .C(use_fxn,
                  as.double(coordinates(map)[Wolv,1]),               #x_wolv
                  as.double(coordinates(map)[Wolv,2]),               #y_wolv
                  as.integer(length(Wolv)),                       #N_wolv

                  as.double(coordinates(map)[,1]),                #x_snow
                  as.double(coordinates(map)[,2]),                #y_snow
                  use = as.double(getValues(map)),                     #snow
                  as.integer(length(getValues(map))),                        #pixels
                  
              as.double(c(sd_xy[,1])),                # double sd_x[]
              as.double(c(sd_xy[,2])),                # double sd_y[]
              as.double(c(trunk))                   # double trunc_cutoff[]
              )$use
  return(USE)
  }
#calc_prob(double use[], int grid[], bool detection[], double *detectionP, int *pixels, int *max_grid)  
probPRES<-function(surface, gridvec, detP=1){
  #tmp = as.double(runif(n=(max(gridvec)+1)))
  USE = .C(calc_prob,
                  as.double(surface),                                   #use  =  use surface with prob of at least one wolverine
                  as.integer(gridvec),                                  #grid = vector of grid indices
                  use = as.integer(rep(0,max(gridvec)+1)),              #detection

                  as.double(detP),                                      #detection probability
                  as.integer(length(surface)),                          #pixels
                  as.integer(max(gridvec)),                             #max_grid
                  occP = as.double(runif(n=(max(gridvec)+1)))          #test values, output detection probabilities
                 )$occP[(unique(gridvec)+1)[-1]]
  #tmp$use[(unique(gridvec)+1)[-1]]   ##OUTPUT: R index number = grid number (0 = out of grid dropped)
  #tmp = cbind(unique(gridvec)[-1],USE$use[unique(gridlist[[4]])+1][-1],USE$occP[unique(gridlist[[4]])+1][-1], tmp[unique(gridlist[[4]])+1][-1])
  #tmp = cbind(USE$use,USE$occP, tmp)
  return(USE)
  }
    

#make_grid(double x[], double y[], double *grid_size, int *pixels, int grid[])
#filter_grid(int grid[], double snow[], double *cutoff, int *pixels, double *snow_cutoff)
make.grid <- function(map, gridsize, cutoff=0, snow_cutoff=NULL, filtered=T){      #function for individual calls to C++
   n = length(getValues(map))
   x = coordinates(map)[,1]
   y = coordinates(map)[,2]

   if(grepl('+proj=longlat',proj4string(map))){
     USE = .C(make_grid,
                    as.double(x),                     #Longitude
                    as.double(y),                     #Latitude
                    as.double(gridsize),              #Desired grid size (100km2)
                    as.integer(n),                    # #pixels
                    gridvec = as.integer(rep(0,n))   # grid vector
                  )$gridvec
   } else USE = make.grid.utm(map, gridsize)   
       
   if(filtered){
   USE = .C(filter_grid,
                  gridvec=as.integer(USE),          # grid vector
                  as.double(getValues(map)),           # snow
                  as.double(cutoff),                # % snow pixels in cell required
                  as.integer(n),                    # #pixels
                  as.double(snow_cutoff)           # habitat value for pixel to be considered habitat (if(snow[i] >= *snow_cutoff))
                )$gridvec }
   return(USE) }

   

##void reduce_pop(int IN[], int N[], double useTotal[], int *pixels, int *snowpoints,
##				double x[], double y[], double snow[], int *n_grps, double *lmda,
##				double sd_x[], double sd_y[])   
#reduce_pop <- function(Wolverines, pts, lmda, howmuch, howfar, trunk, snow_cutoff){      #function for individual calls to C++
#   n = length(getValues(pts))
#   N_wolv = length(Wolverines)
#   n_grps = length(levels(Wolverines$sex))
#   snowpoints = which(getValues(pts) > snow_cutoff)
#
#   IN = rep(0, n_grps*length(snowpoints))
#   b = paste(coordinates(pts)[snowpoints,1],coordinates(pts)[snowpoints,2])
#   for(st in 1:n_grps){
#      use = Wolverines$sex == levels(Wolverines$sex)[st]
#      a = paste(coordinates(Wolverines)[use,1],coordinates(Wolverines)[use,2])
#      IN[which(b %in% a)+(st-1)*length(snowpoints)] = 1  }
#
#   sd_xy<-var_solver(howmuch, howfar, map)
#
#   USE = .C(rd_pop,
#              IN  = as.integer(IN),                       #int IN[]
#              N   = as.integer(table(Wolverines$sex)),    #int N[]
#              use = as.double(getValues(pts)),               #double useTotal[]
#              as.integer(n),                              #int *pixels,
#              as.integer(length(snowpoints)),             #int *snowpoints
#              as.double(coordinates(pts)[,1]),            #double x[]
#              as.double(coordinates(pts)[,2]),            #double y[]
#              as.double(pts$band1),                       #double snow[]
#              as.integer(n_grps),                         # int *n_grps
#              as.double(lmda),                            # double *lmda
#              as.double(c(sd_xy[,1])),                # double sd_x[]
#              as.double(c(sd_xy[,2])),                # double sd_y[]
#              as.double(trunk),                       # double trunc_cutoff[]
#              as.double(snow_cutoff),                 # double *snow_cutoff
#                )
#   return(USE) }
#
### R Subfunctions #### 
second.filter<-function(grid_layer, map, condition=1, cutoff=1){
   V1<-table(grid_layer[getValues(map)>=condition])[-1]
   tmp<-as.numeric(names(V1)[V1>=cutoff*max(V1)])
   grid_layer[!(grid_layer %in% tmp)]=0
   return(grid_layer)
   }
   
third.filter<-function(grid_layer, map, sample.area=1.5){
  small.grid<-make.grid(map, gridsize=sample.area, cutoff=0, filtered=F)
  for(x in unique(grid_layer[grid_layer>0])){
    V1<-table(small.grid[grid_layer==x])
    V1<-V1[names(V1)!='0']
    cutoff<-as.numeric(names(which.max(table(V1))))
    use<-as.numeric(names(V1)[V1>=cutoff])
    use<-sample(use,1)
    grid_layer[grid_layer==x & small.grid!=use]=0
    }
  return(grid_layer)
    }
    
make.grid.utm<-function(map, gridsize){
  map_xy<-dim(map)[2:1] # x = columns, y=rows; values by row
  nxy=ceiling(sqrt(gridsize)*1000/res(map))
  skip_xy<-map_xy %% nxy
  n_cells<-(map_xy-skip_xy) %/% nxy
  
  grid_layer<-expand.grid(x=1:map_xy[1], y=1:map_xy[2], value=0)
  keeppixels<-(grid_layer$x > floor(skip_xy[1]/2)) &
              (grid_layer$x <= map_xy[1]-ceiling(skip_xy[1]/2)) &
              (grid_layer$y > floor(skip_xy[2]/2)) &
              (grid_layer$y <= map_xy[2]-ceiling(skip_xy[2]/2))
  grid_layer$value[keeppixels]<-
    rep(rep(rep(1:n_cells[1], each=nxy[1]),nxy[2]), n_cells[2])+
    rep(seq(from=0,by=n_cells[1], length.out=n_cells[2]), each=(n_cells[1]*nxy[1]*nxy[2]))
  
  return(grid_layer$value)
  }
  
  
  
      
  

var_solver<-function(howmuch, howfar, map){ 
  if(grepl('+proj=utm', proj4string(map))){  # km per map unit conversion
    km_per_x = 1/1000       
    km_per_y = 1/1000
  } else if(grepl('+proj=longlat', proj4string(map))){
    center=c(mean(c(slot(extent(map), 'xmin'), slot(extent(map),'xmax'))),
             mean(c(slot(extent(map), 'ymin'), slot(extent(map),'ymax'))))
    km_per_x = pointDistance(center,center+c(1,0), lonlat=T)/1000
    km_per_y = pointDistance(center,center+c(0,1), lonlat=T)/1000   
  }
    
  cutoff=qnorm((1+howmuch)/2)  # Cutoff for a standard normal distribution       
  howfar_x<-howfar/km_per_x    # Convert km to map units 
  howfar_y<-howfar/km_per_y   
  sd_x <- howfar_x/cutoff      # Convert cutoff to map scale to get sd 
  sd_y <- howfar_y/cutoff
  
   return(cbind(sd_x,sd_y))   
}
        
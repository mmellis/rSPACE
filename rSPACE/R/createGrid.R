# Create grid layer
createGrid<-function(map, pList, filter.map=NULL){
  reNumber=T
  if(is.null(filter.map)){
      gridLayer<-makeGrid(map, gridsize=pList$grid_size)
      gridLayer<-filterByHabitat(gridLayer, map, pList$sample.cutoff, pList$habitat.cutoff)
  } else {
    if(any(is.nan(getValues(filter.map))))
      stop('NaNs in filter map')

    filter.map <- reclassify(filter.map, cbind(NA, 0))

    if(all(getValues(filter.map) %in% c(0,1))){
      gridLayer <- makeGrid(map, gridsize=pList$grid_size)
      gridLayer <- filterByHabitat(gridLayer, map, pList$sample.cutoff, pList$habitat.cutoff)
      gridLayer <- filterByMap(gridLayer, filter.map, pList$filter.cutoff)
    } else {
      gridLayer <- getValues(filter.map)
      reNumber <- F
    }
  }

  if(!is.null(pList$Effective.sample.area))
      gridLayer<-reduceArea(gridLayer, pList$grid_size, pList$Effective.sample.area)

  if(reNumber)
    gridLayer<-(match(gridLayer,unique(c(0, gridLayer)))-1)

  return(gridLayer)
}



# Make rectangular grid -----------------------------------------
makeGrid<-function(map, gridsize){
   n <- length(getValues(map))
   x <- coordinates(map)[, 1]
   y <- coordinates(map)[, 2]
   
   isLongLat <- grepl('+proj=longlat', proj4string(map))
   isUTM     <- grepl('+proj=utm', proj4string(map))

   if(isLongLat){
     use <- .C('make_grid',
                    as.double(x),                     #Longitude
                    as.double(y),                     #Latitude
                    as.double(gridsize),              #Desired grid size (100km2)
                    as.integer(n),                    # #pixels
                    gridvec = as.integer(rep(0,n))    # grid vector
                  )$gridvec
   } else if(isUTM){
       map_xy<-dim(map)[2:1] # x = columns, y=rows; values by row
       nxy=ceiling(sqrt(gridsize)*1000/res(map))
       skip_xy<-map_xy %% nxy
       n_cells<-(map_xy-skip_xy) %/% nxy
  
      gridDF<-expand.grid(x=1:map_xy[1], y=1:map_xy[2], value=0)
      keeppixels<-(gridDF$x > floor(skip_xy[1]/2)) &
                  (gridDF$x <= map_xy[1]-ceiling(skip_xy[1]/2)) &
                  (gridDF$y > floor(skip_xy[2]/2)) &
                  (gridDF$y <= map_xy[2]-ceiling(skip_xy[2]/2))
      gridDF$value[keeppixels]<-
        rep(rep(rep(1:n_cells[1], each=nxy[1]),nxy[2]), n_cells[2])+
        rep(seq(from=0,by=n_cells[1], length.out=n_cells[2]), each=(n_cells[1]*nxy[1]*nxy[2]))
  
      use <- gridDF$value
  } else stop('Map must be either longlat or utm')
  
  return(use) 
}

# 1st filter: meet minimum habitat threshold --------------------
filterByHabitat<-function(gridLayer, map, sample.cutoff, habitat.cutoff){  
  habitatOK <- getValues(map) >= habitat.cutoff  
  V1 <- table(gridLayer[habitatOK])[-1]
  V2 <- table(gridLayer)[-1]
  matchV1<-match(names(V1), names(V2)) 
  
  keep    <- (V1 >= sample.cutoff * V2[matchV1])
  keepIds <- as.numeric(names(V1))[keep]
  
  gridLayer[!(gridLayer %in% keepIds)] <- 0
  return(gridLayer)
}

# 2nd filter: apply filter map ----------------------------------
filterByMap<-function(gridLayer, filter.map, cutoff=0.95){
   mapOK <- getValues(filter.map)==1
   V1 <- table(gridLayer[mapOK])[-1]
   V2 <- table(gridLayer)[-1]
   matchV1 <- match(names(V1),names(V2))

   keep    <- (V1 >= cutoff*V2[matchV1])
   keepIds <- as.numeric(names(V1))[keep]

   gridLayer[!(gridLayer %in% keepIds)] <- 0
   return(gridLayer)
}

# 3rd filter: reduce effective sampling area --------------------
reduceArea<-function(gridLayer, grid_size, sample.area){
  if(sample.area < grid_size){
    for(x in unique(gridLayer[gridLayer>0])){
      IDs<-which(gridLayer==x)
        big.nr<-min(which(diff(IDs)>1))
        big.nc<-round(length(IDs)/big.nr)
      Cell<-matrix(IDs[1:c(big.nr*big.nc)], nrow=big.nr)

      keep <- sample.area/grid_size
        little.nr<- max(c(1,round(keep * big.nr)))
        little.nc<- max(c(1,round(keep * big.nc)))

      Start<-sample((big.nr-little.nr)*(big.nc-little.nc),1)
        start.nr<- Start %% (big.nr-little.nr)
          start.nr<-ifelse(start.nr==0,big.nr-little.nr,start.nr)
        start.nc<- (Start-1)%/%(big.nr-little.nr) + 1

      KeepID<-c(Cell[seq(start.nr, length.out=little.nr), seq(start.nc, length.out=little.nc)])
      DropIDs<-IDs[!(IDs %in% KeepID)]

      gridLayer[DropIDs]<-0
      }}

  return(gridLayer)
}
# Create grid layer
createGrid<-function(map, pList, filter.map=NULL){
  reNumber=T
  if(is.null(filter.map)){
      gridLayer<-makeGrid(map, gridsize=pList$grid_size)
      gridLayer<-filterByHabitat(gridLayer, map, pList$sample.cutoff, pList$habitat.cutoff)
  } else {
    if(any(is.nan(getValues(filter.map))))
      stop('NaNs in filter map')

    if(any(is.na(getValues(filter.map))))
      stop('NAs in filter map. Replace with 0s')

    if(all(getValues(filter.map) %in% c(0,1))){
      gridLayer <- makeGrid(map, gridsize=pList$grid_size)
      gridLayer <- filterByHabitat(gridLayer, map, pList$sample.cutoff, pList$habitat.cutoff)
      gridLayer <- filterByMap(gridLayer, filter.map, pList$filter.cutoff)
    } else {
      gridLayer <- getValues(filter.map)
      reNumber <- F
    }
  }

  if(length(unique(gridLayer))==1)
    stop('No grid cells in layer. Check filtering steps')

  if(!is.null(pList$Effective.sample.area))
      gridLayer<-reduceArea(gridLayer, pList$grid_size, pList$Effective.sample.area)

  if(reNumber)
    gridLayer<-(match(gridLayer,unique(c(0, gridLayer)))-1)

  return(gridLayer)
}



# Make rectangular grid -----------------------------------------
makeGrid<-function(map, gridsize, type='square'){
  if (!type %in% c("square", "hexagonal")) {
    stop("Type must be either 'square' or 'hexagonal'")}
 
 x<-as(extent(map), 'SpatialPolygons')

 if(grepl('longlat', projection(map))){
  stop("Please switch to planar projection") }
    
 if (type == "square") {
    cell_width <- sqrt(gridsize)    # Must match units of raster
  } else if (type == "hexagonal") {
    cell_width <- sqrt(2 * gridsize / sqrt(3))
  }
  
  projection(x) <- projection(map)

 if (type == "square") {
    g <- raster(x, resolution = cell_width )
    g <- as(g, "SpatialPolygons")
  } else if (type == "hexagonal") {
    g <- spsample(x, type = "hexagonal", cellsize = cell_width, offset = c(0, 0))
    g <- HexPoints2SpatialPolygons(g, dx = cell_width)
  }

  g <- rgeos::gIntersection(g, x, byid = TRUE)
  g <- g[round(area(g)) == max(round(area(g))), ]

  return(g) 
}

# 1st filter: meet minimum habitat threshold --------------------
filterByHabitat<-function(gridLayer, map, sample.cutoff, habitat.cutoff){  
  if(sample.cutoff == 0 ) return(gridLayer)

  habitatOK <- getValues(map) >= habitat.cutoff

  IDs <- sort(unique(gridLayer))
  if(IDs[1] == 0)
    IDs <- IDs[-1]
   
  V1 <- tabulate(gridLayer[habitatOK], nbins=max(IDs))[IDs]
  V2 <- tabulate(gridLayer, nbins=max(IDs))[IDs]
  
  keep    <- IDs[V1 >= sample.cutoff * V2]
  
  gridLayer[!(gridLayer %in% keep)] <- 0
  return(gridLayer)
}

# 2nd filter: apply filter map ----------------------------------
filterByMap<-function(gridLayer, filter.map, cutoff=0.95){
   if(cutoff == 0) return(gridLayer)

   mapOK <- getValues(filter.map)==1

   Cells<-sort(unique(gridLayer))
   if(Cells[1] == 0)
    Cells<-Cells[-1]

   V1 <- tabulate(gridLayer[mapOK], nbins=max(Cells))[Cells]
   V2 <- tabulate(gridLayer, nbins=max(Cells))[Cells]

   keep    <- Cells[V1 >= cutoff*V2]

   gridLayer[!(gridLayer %in% keep)] <- 0
   return(gridLayer)
}

# 3rd filter: reduce effective sampling area --------------------
reduceArea<-function(gridLayer, grid_size, sample.area){
  if(sample.area > grid_size)
    stop('Effective sample area is larger than grid_size')

  if(sample.area < grid_size){
    for(x in unique(gridLayer[gridLayer>0])){
      IDs<-which(gridLayer==x)
        big.nr<-min(which(diff(IDs)>1))
        big.nc<-round(length(IDs)/big.nr)
      Cell<-matrix(IDs[1:c(big.nr*big.nc)], nrow=big.nr)

      keep <- sqrt(sample.area/grid_size)
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
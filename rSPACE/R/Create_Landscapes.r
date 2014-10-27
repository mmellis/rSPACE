### Functions for creating encounter history files
default.value<-function(x,val) ifelse(is.null(x),val,x) 

# Create grid_layer
create.grid<-function(map, Parameters, filter.map=NULL){
  reNumber=T 
  if(is.null(filter.map)) {
      grid_layer<-make.grid(map, gridsize=Parameters$grid_size, 
                    cutoff=Parameters$sample.cutoff, 
                    snow_cutoff=Parameters$HRcenter.cutoff)
    } else {
      filter.map <- reclassify(filter.map, cbind(NA, 0))
      filter.map <- extend(filter.map, extent(map), value=0)
      if(any(is.nan(getValues(filter.map)))) stop('NaNs in filter map')
  
      if(all(getValues(filter.map) %in% c(0,1))){
        grid_layer<-make.grid(map, Parameters$grid_size, filtered=F)
        grid_layer<-second.filter(grid_layer, map, filter.map)
      }else{ 
        grid_layer<-getValues(filter.map) 
        reNumber=F
      }
    }
    if(reNumber) 
      grid_layer<-(match(grid_layer,unique(c(0, grid_layer)))-1)
    
    return(grid_layer)
  }
    

# Add wolverines to simulation
add.wolv<-function(dN, map, Parameters, wolv.df = NULL){
  if(is.null(Parameters$wghts)) Parameters$wghts=F
  if(is.null(wolv.df)){
   use<-which(getValues(map)>=Parameters$HRcenter.cutoff)
   new.wolv<-lapply(1:length(Parameters$MFratio), function(x) {use[smple(map,use,dN*Parameters$MFratio[x],Parameters$buffer[x],Parameters$wghts)]})
  } else {
   use<-lapply(1:length(Parameters$MFratio), function(x) {which(getValues(map)>=Parameters$HRcenter.cutoff & getValues(distanceFromPoints(map, coordinates(map)[wolv.df[wolv.df$sex==x,]$locID,]))>Parameters$buffer[x])}) 
   new.wolv<-lapply(1:length(Parameters$MFratio), function(x) use[[x]][smple(map, use[[x]], round(dN*Parameters$MFratio[x]),Parameters$buffer[x])])      
  }
  return(new.wolv)
}

# Remove wolverines from simulation
drop.wolv<-function(dN, map, wolv.df, Parameters){
  if(is.null(Parameters$wghts)) Parameters$wghts=F
  if(Parameters$wghts==T){
    drop.rows<-sample(nrow(wolv.df),dN, prob=getValues(map)[wolv.df$locID])
  } else { drop.rows<-sample(nrow(wolv.df),dN) }
    
  lost.wolv<-wolv.df[drop.rows,]
  lost.wolv<-lapply(1:length(Parameters$MFratio), function(x) lost.wolv$locID[lost.wolv$sex==x])
  return(list(drop.rows, lost.wolv))
}

# Convert wolv list to data.frame
wolv.dataframe<-function(wolv.list, map=NULL){
  wolv.df<-data.frame(sex=rep(1:length(wolv.list),sapply(wolv.list,length)),locID=unlist(wolv.list))
  if(!is.null(map)){
    xy<-coordinates(map)[wolv.df$locID,]
    wolv.df$x<-xy[,1]
    wolv.df$y<-xy[,2]
    }
  return(wolv.df)
}
  
 
# Create a use surface from a list of wolverine locations 
build.useLayer<-function(map, wolv, Parameters, Example=F){
  useLayer<-1-use_surface(wolv[[1]],Parameters$howmuch[1], Parameters$howfar[1], map, trunk = Parameters$trunk[1])
  if(!Example & length(Parameters$MFratio)>1){
    for(ii in 2:length(Parameters$MFratio)){ 
      useLayer<-useLayer*(1-use_surface(wolv[[ii]],Parameters$howmuch[ii], Parameters$howfar[ii], map, trunk = Parameters$trunk[ii]))}
  }
  useLayer<-1-useLayer 
  return(useLayer)
}  

# Main simulation function...simulate landscape and create encounter history
encounter.history<-function(map, Parameters, ...){
  add.args<-list(...)
    n_cells<-add.args$n_cells
    grid_layer<-add.args$grid_layer
    filter.map<-add.args$filter.map
    showSteps<-default.value(add.args$showSteps, F)
    printN<-default.value(add.args$printN, 0)

  n_visits<-Parameters$n_visits
  n_yrs<-Parameters$n_yrs
  
  lmda<-rep(Parameters$lmda, length.out=(n_yrs-1))

  if(is.null(grid_layer)){
    grid_layer<-create.grid(map, Parameters, filter.map)
    n_cells<-length(unique(grid_layer)[unique(grid_layer)>0])
  }

  encounter_history<-matrix(0, nrow=n_cells, ncol=n_yrs)

  # 1. Place individuals
  wolv<-add.wolv(Parameters$N, map, Parameters)
  wolv.df<-wolv.dataframe(wolv)

  # 2. Make use surface
  useLayer<-build.useLayer(map, wolv, Parameters)
  
  # 3. Calculate probability by grid
  P.pres<-matrix(0,nrow=n_cells, ncol=n_yrs)
  P.pres[,1]<-probPRES(useLayer, grid_layer)
 
  # 3b. (Optional) Output plots with first year data
  if(showSteps){
    printN<-""  # Cause printN to print to console
    
    prev.ask<-devAskNewPage(ask=T)
    # Habitat map
    plot(map); title('User supplied habitat map')
    
    # Initial locations
    plot(map)
      wolv.points<-wolv.dataframe(wolv, map)
    points(wolv.points$x, wolv.points$y, pch=c(15:18)[wolv.points$sex], col='black')
    title('Activity center locations - Year 1')
    rm('wolv.points')
    
    # Grid
    plot(setValues(map, grid_layer), col=c("white", sample(26:137, n_cells, replace=T)))
    title('Grid (colored for contrast only)')
    
    # Individual use layer example
    ExampleLayer<-build.useLayer(map, list(wolv[[1]][sample(length(wolv[[1]]),1)]), Parameters, Example=T)
    plot(setValues(map,ExampleLayer))
    title('Example individual use layer')
    rm('ExampleLayer')
    
    # Full use surface
    plot(setValues(map,useLayer))
    title('Probability of at least one individual - Year 1') 
    
    # Gridded use
    plot(setValues(map, c(0,P.pres[,1])[match(grid_layer,unique(c(0, grid_layer)))]))
    title('Probability of at least one individual by cell - Year 1')
    
    prev.ask<-devAskNewPage(ask=prev.ask)
    rm('prev.ask')
    
    cat('\nTotal number of individuals by year\n')
  } 
 
  # 4. Sample detections in the first year
  encounter_history[,1]<-sapply(1:n_cells, 
    function(i) paste(rbinom(n=n_visits, size=1, prob=P.pres[i,1]), collapse=''))
  
  # 5. Loop over years to fill in encounter_history
  if(n_yrs > 1){  
   for(tt in 2:n_yrs){
      if(printN!=0)
       cat(nrow(wolv.df),' ',sep='',file=printN,append=T)  # Store population sizes by year
      
      # 6. Calculate population change between t and t+1 
      dN<-round(nrow(wolv.df)*(lmda[tt-1] - 1))
        
      # 7. Implement population change  
       if(dN>0){
           new.wolv<-add.wolv(dN, map, Parameters, wolv.df)
           useLayer<-1-(1-useLayer)*(1-build.useLayer(map, new.wolv, Parameters))
           wolv.df<-rbind(wolv.df, wolv.dataframe(new.wolv))
       } else if(dN<0) {
          lost.wolv<-drop.wolv(abs(dN), map, wolv.df, Parameters)
          useLayer<-1-(1-useLayer)/(1-build.useLayer(map, lost.wolv[[2]], Parameters))
          wolv.df<-wolv.df[-lost.wolv[[1]],]
       }
    
      # 8. Sample detections and update encounter_history
      P.pres[,tt]<-probPRES(useLayer, grid_layer)

      encounter_history[,tt]<-sapply(1:n_cells, 
          function(i) paste(rbinom(n=n_visits, size=1, prob=P.pres[i,tt]), collapse=''))
    }} # End year loop (and if statement)


  if(printN!=0) 
    cat(nrow(wolv.df),'\n',sep='',file=printN,append=T)  # Initial population size

  if(showSteps){
    cat('\nOutputting P.pres and encounter histories by year...\n')
    encounter_history<-data.frame(round(P.pres,3), encounter_history)
      names(encounter_history)<-paste(rep(paste0('Yr', 1:n_yrs),2), 
                                      rep(c('p','ch'),each=n_yrs), sep='.')
      encounter_history<-encounter_history[,order(rep(1:n_yrs,2))]
    return(encounter_history)
  } else return(apply(encounter_history,1,paste, collapse=''))
}  

# Landscape wrapper
create.landscapes<-function(n_runs, map, Parameters, ... ){
  # 0. Match argument list 
  additional.args<-list(...)
    folder.dir<-default.value(additional.args$folder.dir, getwd())    
    run.label<- default.value(additional.args$run.label, 'rSPACE_X')
    base.name<- default.value(additional.args$base.name, 'rSPACEx')
    filter.map<-additional.args$filter.map
    printN<-default.value(additional.args$printN, 1)
    saveParameters<-default.value(additional.args$saveParameters, 1)
    skipConfirm<-default.value(additional.args$skipConfirm, F)

  # 0. Set up files
   if(!skipConfirm){
      askConfirm<-("" ==readline(prompt="\n rSPACE creates text files.
        If you're ok with this, press ENTER to continue.
        Typing anything else will exit.\n"))
      if(!askConfirm){
        message('Exiting function')
        return(0)
      }
   }

   folder.dir <- paste(folder.dir, run.label, sep='/') 
   if(!file.exists(folder.dir)) dir.create(folder.dir)


   output.dir <- paste(folder.dir, 'output', sep='/')
   if(!file.exists(output.dir)) dir.create(output.dir)
   if(printN) printN<-paste0(output.dir,'/N_final.txt')

  
  # 1. Enter parameters
  if(missing(Parameters)) {Parameters<-enter.parameters()}
  if(is.null(Parameters$detP)) Parameters$detP=1
  if(is.null(Parameters$trunk)) Parameters$trunk=rep(0,length(Parameters$MFratio))

  # 2. Set up map + grid layer
  if(missing(map)) stop("Missing habitat layer")
  if(!grepl('proj=utm|proj=longlat',proj4string(map))) stop("Projection needs to be in utm or longlat")
  if(grepl('+proj=utm.*',proj4string(map))) 
    if(!grepl('+units=m',proj4string(map))) 
      message('Assuming UTM +units=m')   
  map <- reclassify(map, cbind(c(NA), c(0)))
  if(any(is.nan(getValues(map)))) stop('NaNs in habitat map')

  grid_layer<-create.grid(map, Parameters, filter.map)
  gridIDs<-unique(grid_layer)[unique(grid_layer)>0]

  # 3. Simulate encounter histories loop ##
  for(rn in 1:n_runs){
    cat(rn,'\n');flush.console()
    ch<-encounter.history(map, Parameters, grid_layer=grid_layer, n_cells=length(gridIDs), printN=printN)
  
    # 4. Output encounter history
    output_file<-paste(folder.dir,'/',base.name,rn,".txt",sep='')
    cat(paste("/*", gridIDs, "*/", ch, "1;"), sep="\n",file=output_file) 
  } # End runs loop
  
  if(saveParameters==1)
    save(Parameters, file=paste0(output.dir,'/Parameters.Rdata'))
    
  return(list(DIR = folder.dir, 
          filenames=paste0(base.name,1:n_runs,'.txt')))
} # End function

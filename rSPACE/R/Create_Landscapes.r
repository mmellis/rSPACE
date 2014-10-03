### Functions for creating encounter history files
default.value<-function(x,val) ifelse(is.null(x),val,x) 

# Add wolverines to simulation
add.wolv<-function(dN, map, Parameters, wolv.df = NULL){
  if(is.null(Parameters$wghts)) Parameters$wghts=F
 if(is.null(wolv.df)){
   use<-which(getValues(map)>=Parameters$HRcenter.cutoff)
   new.wolv<-lapply(1:length(Parameters$MFratio), function(x) {use[smple(map,use,dN*Parameters$MFratio[x],Parameters$buffer[x],Parameters$wghts)]})
 } else {  
   use<-sapply(1:length(Parameters$MFratio), function(x) {which(getValues(map)>Parameters$HRcenter.cutoff & getValues(distanceFromPoints(map, coordinates(map)[subset(wolv.df,sex==x)$locID,]))>Parameters$buffer[x])}) 
   new.wolv<-sapply(1:length(Parameters$MFratio), function(x) use[[x]][smple(map, use[[x]], round(dN*Parameters$MFratio[x]),Parameters$buffer[x])])      
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
build.useLayer<-function(map, wolv, Parameters){
  useLayer<-1-use_surface(wolv[[1]],Parameters$howmuch[1], Parameters$howfar[1], map, trunk = Parameters$trunk[1])
  if(length(Parameters$MFratio)>1){
    for(ii in 2:length(Parameters$MFratio)){ 
      useLayer<-useLayer*(1-use_surface(wolv[[ii]],Parameters$howmuch[ii], Parameters$howfar[ii], map, trunk = Parameters$trunk[ii]))}
  }
  useLayer<-1-useLayer 
  return(useLayer)
}  

# Main simulation function...simulate landscape and create encounter history
encounter.history<-function(Parameters, map, grid_layer, n.cells, printN=NULL){
  
  encounter_history<-matrix(0,nrow=n.cells, ncol=Parameters$n_yrs*Parameters$n_visits)
  ch.rplc<-lapply(seq(1,(Parameters$n_yrs)*Parameters$n_visits, by=Parameters$n_visits), function(x) x:(x+Parameters$n_visits-1)) 

  # 1. Place individuals
 wolv<-add.wolv(Parameters$N, map, Parameters)
 wolv.df<-wolv.dataframe(wolv)

  # 2. Make use surface
 useLayer<-build.useLayer(map, wolv, Parameters)
  
  # 3. Calculate probability by grid
 P.pres<-probPRES(useLayer, grid_layer)
 P.pres[P.pres<0]=0
 
  # 4. Sample detections in the first year
 encounter_history[,ch.rplc[[1]]]<-replicate(Parameters$n_visits, rbinom(n=length(P.pres),size=1, prob=P.pres))
 
if(Parameters$n_yrs>1){  # 5. Loop over years to fill in encounter_history
   for(tt in 2:Parameters$n_yrs){
      if(printN!=0)
       cat(nrow(wolv.df),' ',sep='',file=printN,append=T)  # Store population sizes by year
      # 6. Calculate population change between t and t+1 
      if(length(Parameters$lmda)>1) {lmda.tt<-Parameters$lmda[tt-1]} else {lmda.tt<-Parameters$lmda}
      dN<-lmda.tt*nrow(wolv.df)-nrow(wolv.df)
        
      # 7. Implement population change  
       if(round(dN)>0){
           new.wolv<-add.wolv(round(dN), map, Parameters, wolv.df)
           useLayer<-1-(1-useLayer)*(1-build.useLayer(map, new.wolv, Parameters))
           wolv.df<-rbind(wolv.df, wolv.dataframe(new.wolv))
       } else if(round(dN)<0) {
          lost.wolv<-drop.wolv(round(-dN), map, wolv.df, Parameters)
          useLayer<-1-(1-useLayer)/(1-build.useLayer(map, lost.wolv[[2]], Parameters))
          wolv.df<-wolv.df[-lost.wolv[[1]],]
       }
    
      # 8. Sample detections and update encounter_history
      P.pres<-probPRES(useLayer, grid_layer)
      P.pres[P.pres<0]=0
      encounter_history[,ch.rplc[[tt]]]<-replicate(Parameters$n_visits, rbinom(n=length(P.pres),size=1, prob=P.pres))
    
    }} # End year loop (and if statement)

  if(printN!=0) 
    cat(nrow(wolv.df),'\n',sep='',file=printN,append=T)  # Initial population size

  return(encounter_history)
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

  # 0. Set up files
   folder.dir <- paste(folder.dir, run.label, sep='/') 
   if(!file.exists(folder.dir)) dir.create(folder.dir)
   if(printN) printN<-paste0(folder.dir,'/N_final.txt')
  
  # 1. Enter parameters
  if(missing(Parameters)) {Parameters<-enter.parameters()}
  if(is.null(Parameters$detP)) Parameters$detP=1
  if(is.null(Parameters$trunk)) Parameters$trunk=rep(0,length(Parameters$MFratio))

  # 2. Set up map + grid layer
  if(missing(map)) stop("Missing habitat layer")
  if(grepl('+proj=utm.*',proj4string(map))) 
    if(!grepl('+units=m',proj4string(map))) 
      message('Assuming UTM +units=m')   
  map <- reclassify(map, cbind(c(NA), c(0)))
  if(any(is.nan(getValues(map)))) stop('NaNs in habitat map')

  if(is.null(filter.map)) {
    grid_layer<-make.grid(map, gridsize=Parameters$grid_size, 
                  cutoff=Parameters$sample.cutoff, 
                  snow_cutoff=Parameters$HRcenter.cutoff)
  } else {
    filter.map<-reclassify(filter.map, cbind(NA,0))
    filter.map<-setValues(filter.map, ifelse(!(getValues(filter.map)==0),1,0))
    
    grid_layer<-make.grid(map, Parameters$grid_size, filtered=F)
    grid_layer<-second.filter(grid_layer, map, filter.map)
  }
  
  gridIDs<-unique(grid_layer)[unique(grid_layer)>0]

  # 3. Simulate encounter histories loop ##
  for(rn in 1:n_runs){
    cat(rn,'\n');flush.console()
    encounter_history<-encounter.history(Parameters, map, grid_layer, length(gridIDs), printN)
  
    # 4. Output encounter history
    output_file<-paste(folder.dir,'/',base.name,rn,".txt",sep='')
    ch<-apply(encounter_history, 1, function(x) paste(x,collapse=""))
    cat(paste("/*", gridIDs, "*/", ch, "1;"), sep="\n",file=output_file) 
  } # End runs loop
  
  if(saveParameters==1)
    save(Parameters, file=paste0(folder.dir,'/Parameters.Rdata'))
    
  return(list(DIR = folder.dir, 
          filenames=paste0(base.name,1:n_runs,'.txt')))
} # End function

### Functions for creating encounter history files
# Landscape wrapper
create.landscapes<-function(n_runs, map, Parameters, ... ){
  # 0. Match argument list 
  additional.args<-list(...)
    folder.dir<-setDefault(additional.args$folder.dir, getwd())    
    run.label<- setDefault(additional.args$run.label, 'rSPACE_X')
    base.name<- setDefault(additional.args$base.name, 'rSPACEx')
    filter.map<-additional.args$filter.map
    printN<-setDefault(additional.args$printN, 1)
    saveParameters<-setDefault(additional.args$saveParameters, 1)
    skipConfirm<-setDefault(additional.args$skipConfirm, F)

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
  Parameters<-checkParameters(Parameters, additional.args)

  # 2. Set up map + grid layer
  if(missing(map)) stop("Missing habitat layer")
  if(!grepl('proj=utm|proj=longlat',proj4string(map))) stop("Projection needs to be in utm or longlat")
  if(grepl('+proj=utm.*',proj4string(map))) 
    if(!grepl('+units=m',proj4string(map))) 
      message('Assuming UTM +units=m')   
  map <- reclassify(map, cbind(c(NA), c(0)))
  if(any(is.nan(getValues(map)))) stop('NaNs in habitat map')

  grid_layer<-createGrid(map, Parameters, filter.map)
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

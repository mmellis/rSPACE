### Subfunctions for analyzing output files with SPACE ###
########## Subfunctions ################
getData<-function(folder, CI=0.95){
    filename<-dir(folder, pattern='sim_results.*txt', recursive=T, full.names=T)
    key <- dir(folder, pattern='Parameters.Rdata', recursive=T, full.names=T)
    if(length(filename)==0) stop("Can't find the file")
    if(length(key)==0) stop("Can't find the parameters key")
    
    if((length(filename)>1)|(length(key)>1))
      warning('More than 1 file. Only the first will be used.')
    
    lmda<-local({
        load(key[1])
        get(ls())$lmda 
      })
    
    CI<-qnorm(CI) 
    dta<-read.table(filename[1], header=T)
    dta$count<-as.numeric(sign(lmda-1)*dta$trend > CI*dta$trendSE) 
 return(dta)
}


## Summing %detected
sumData<-function(dta){
  variables<-c("n_grid","n_visits","detP","gap_yr")
  count <- NULL
  dta<-dta[, names(dta) %in% c(variables, 'count')]
  dta<-ddply(dta, names(dta)[!grepl("count",names(dta))], summarise,
    total=sum(!is.na(count)),
    n_runs=length(count),
    count=sum(count, na.rm=T))
    
    variables<-variables[!(variables %in% names(dta))]
    if(length(variables)>0)
       eval(parse(text=paste0('dta$',variables,'=0')))
       
  return(dta)
}

## Single file plot function
getResults<-function(folder, CI=0.95, returnData=1, plot=T) {  
    dta<-getData(folder, CI) 
    dtaS<-sumData(dta)
    
    count<-n_runs<-n_grid<-n_visits<-detP<-gap_yr<-NULL
          
    if(plot){
      print(ggplot(subset(dtaS, detP<1), 
          aes(x=n_grid, y=(count/n_runs), group=interaction(n_visits,detP,gap_yr)))+
        geom_line(aes(colour=n_visits, linetype=factor(detP)),size=1.25)+
        scale_colour_gradient(name = "# visits",guide="legend") +
        scale_linetype_discrete(name=expression(p["sim"]))+
        scale_y_continuous(limits=c(0,1))+
        labs(x="Number of cells sampled", 
             y="Detected trend/Number of runs", 
             title=basename(folder))+
        facet_wrap(~gap_yr) )
    }
     
    if(returnData==1) return(dta)
    if(returnData==2) return(dtaS)
 }
 
 
# Uses loess smoother to interpolate number of cells required for power
findPower<-function(folder, data, CI=0.95, pwr=.8, n_grid=NULL){
  tryNA<-function(xx) tryCatch(xx, error=function(x) return(NA))
  if(missing(folder) & missing(data)) stop("Must supply either a folder or data")
  if(missing(data)){
    dta<-getData(folder, CI)
    dtaS<-sumData(dta)
  } else { dtaS<-data }
  dtaS$percent<-dtaS$count/dtaS$total

  extra.var<-which(!(names(dtaS) %in% c('n_grid','total','n_runs','count','percent'))&
              unname(apply(dtaS, 2, function(x) length(unique(x))>1)))

  if(is.null(n_grid)){
    fit_loess_pwr<-function(DF, percent){
      fit<-tryNA(loess(n_grid~percent, data=DF))
      fit<-tryNA(as.numeric(predict(fit, data.frame(percent=percent))))
      return(data.frame(n_grid=fit))}
    fit<-ddply(dtaS, names(dtaS)[extra.var], fit_loess_pwr, percent=pwr)
  } else {
    fit_loess_grid<-function(DF, n_grid){
      fit<-tryNA(loess(percent~n_grid, data=DF))
      fit<-tryNA(as.numeric(predict(fit, data.frame(n_grid=n_grid))))
      return(data.frame(pwr=fit))}
    fit<-ddply(dtaS, names(dtaS)[extra.var], fit_loess_grid, n_grid=n_grid)
  }
  return(fit)
}
#########################
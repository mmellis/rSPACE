### Subfunctions for analyzing output files with SPACE ###
########## Subfunctions ################
get_data<-function(folder, CI=0.95){
    filename<-dir(folder, pattern='sim_results.*txt', recursive=T, full.names=T)
    key <- dir(folder, pattern='Parameters.Rdata', recursive=T, full.names=T)
    if(length(filename)==0) stop("Can't find the file")
    if(length(key)==0) stop("Can't find the parameters key")
    
    if((length(filename)>1)|(length(key)>1))
      warning('More than 1 file. Only the first will be used.')
    
    load(key)
    
    CI<-qnorm(CI) 
    dta<-read.table(filename[1], header=T)
    lmda<-Parameters$lmda
    dta$count<-as.numeric(sign(lmda-1)*dta$trend > CI*dta$trendSE) 
 return(dta)
}


## Summing %detected
sum_data<-function(dta){
  dta<-dta[, !(grepl('trend',names(dta)) | grepl('X',names(dta)) | names(dta) %in% c('p_est','rn'))]
  dta<-ddply(dta, names(dta)[!grepl("count",names(dta))], summarise,
    total=sum(!is.na(count)),
    n_runs=length(count),
    count=sum(count, na.rm=T))
  return(dta)
}

## Single file plot function
plot.results<-function(folder, CI=0.95, returnData=1, plotResults=T) {  
    dta<-get_data(folder, CI) 
    dta_sum<-sum_data(dta)
    
    if(plotResults){
      print(ggplot(subset(dta_sum, detP<1), 
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
    if(returnData==2) return(dta_sum)
 }
 
 
# Uses loess smoother to interpolate number of cells required for power
find.power<-function(folder, data, CI=0.95, pwr=.8, n_grid=NULL){
  tryNA<-function(xx) tryCatch(xx, error=function(x) return(NA))
  if(missing(folder) & missing(data)) stop("Must supply either a folder or data")
  if(missing(data)){
    dta<-get_data(folder, CI)
    dta_sum<-sum_data(dta)
  } else { dta_sum<-data }
  dta_sum$percent<-dta_sum$count/dta_sum$total

  extra.var<-which(!(names(dta_sum) %in% c('n_grid','total','n_runs','count','percent'))&
              unname(apply(dta_sum, 2, function(x) length(unique(x))>1)))

  if(is.null(n_grid)){
    fit.loess<-function(DF, percent){
      fit<-tryNA(loess(n_grid~percent, data=DF))
      fit<-tryNA(as.numeric(predict(fit, data.frame(percent=percent))))
      return(data.frame(n_grid=fit))}
    fit<-ddply(dta_sum, names(dta_sum)[extra.var], fit.loess, percent=pwr)
  } else {
    fit.loess<-function(DF, n_grid){
      fit<-tryNA(loess(percent~n_grid, data=DF))
      fit<-tryNA(as.numeric(predict(fit, data.frame(n_grid=n_grid))))
      return(data.frame(pwr=fit))}
    fit<-ddply(dta_sum, names(dta_sum)[extra.var], fit.loess, n_grid=n_grid)
  }
  return(fit)
}
#########################
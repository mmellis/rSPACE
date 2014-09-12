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
    
    CI<-qnorm(CI+0.5*(1-CI))
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
plot.results<-function(folder, CI=0.95, returnData=T) {  
    dta<-get_data(folder, CI) 
    dta_sum<-sum_data(dta)
    
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
     
    if(returnData==T) return(dta)
 }
 
 
# Uses loess smoother to interpolate number of cells required for power
find.power<-function(folder, pwr=.8, CI=0.95, n_grid=NULL){
  dta<-get_data(folder, CI)
  dta_sum<-sum_data(dta)
    
  dta_sum$percent<-dta_sum$count/dta_sum$total
  if(is.null(n_grid)){ 
    fit<-tryCatch(loess(n_grid~percent, data=dta_sum), error=function(x) return(NA))
    fit<-tryCatch(as.numeric(predict(fit, data.frame(percent=pwr))), error=function(x) return(NA))
  } else {
    fit<-tryCatch(loess(percent~n_grid, data=dta_sum), error=function(x) return(NA))
    fit<-tryCatch(as.numeric(predict(fit, data.frame(n_grid=n_grid))), error=function(x) return(NA))
  } 
  return(fit)
}
#########################
### Subfunctions for analyzing output files with SPACE ###
########## Subfunctions ################
## Summing %detected
sum_data<-function(dta){
  dta<-dta[, !(grepl('trend',names(dta)) | grepl('X',names(dta)) | names(dta) %in% c('p_est','rn'))]
  dta<-plyr::ddply(dta, names(dta)[!grepl("count",names(dta))], summarise,
    total=sum(!is.na(count)),
    n_runs=length(count),
    count=sum(count, na.rm=T))
  return(dta)
}

## Single file plot function
plot.results<-function(key, CI=0.95, REML=T) {  # here 'key' should just be the one row you want to plot
    CI<-qnorm(CI+0.5*(1-CI))
    dta<-read.table(key$filename, header=T)
    lmda<-key$lmda
    if(REML==T & !is.null(dta$trend.REML)){dta$count<-as.numeric(sign(lmda-1)*dta$trend.REML > CI*dta$trendSE.REML)
    }else{dta$count<-as.numeric(sign(lmda-1)*dta$trend > CI*dta$trendSE)} 
    dta_sum<-sum_data(dta)
    
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(1,1)))
    use<-droplevels(subset(dta_sum,detP<1 & gap_yr==F))
    d<-ggplot(use, aes(x=n_grid, y=(count/n_runs), group=interaction(n_visits,detP,gap_yr)))
    print(d+geom_line(aes(colour=n_visits, linetype=factor(detP)),size=1.25)+
      scale_colour_gradient(name = "# visits",guide="legend") +
      scale_linetype_discrete(name=expression(p["sim"]))+
      scale_y_continuous(limits=c(0,1))+
      labs(x="Number of cells sampled", y="Detected trend/Number of runs", title=filename),
     vp=viewport(layout.pos.row=1,layout.pos.col=1))
     
    return(dta)
 }
 
## Combine a subset of files and plot facetted graph in gpplot
combine.subset<-function(key, facet.form, CI=0.95, REML=T, subset.fun=NULL, do.plot=T){
  facet.form<-gsub(' ','',facet.form)
  CI<-qnorm(CI+0.5*(1-CI))
  for(i in 1:nrow(key)){  
    dta<-read.table(key$filename[i], header=T)
    if(!is.null(subset.fun)) dta<-subset.fun(dta)
    lmda<-key$lmda[i]
     if(REML & !is.null(dta$trend.REML)) {dta$count<-as.numeric(sign(lmda-1)*dta$trend.REML > CI*dta$trendSE.REML)
     }else{    dta$count<-as.numeric(sign(lmda-1)*dta$trend > CI*dta$trendSE)}
     dta$key<-i  
    if(i==1) {dta_sum<-sum_data(dta)} else {dta_sum<-rbind(dta_sum, sum_data(dta))}}
    
    dta_sum<-data.frame(dta_sum,key[dta_sum$key, -1])
    
    if(do.plot==T){
    facet.sides<-unlist(strsplit(facet.form,'~'))
      facet.sides<-unlist(lapply(facet.sides, function(x) tryCatch(length(unique(dta_sum[,x])),error=function(e) return(1))))
    dev.new(width=1.5+4*facet.sides[2], height=1+3*facet.sides[1])
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(1,1)))
    d<-ggplot(dta_sum, aes(x=n_grid, y=(count/n_runs), group=interaction(n_visits,detP,gap_yr)))
    print(d+geom_line(aes(colour=n_visits, linetype=factor(detP)),size=1.25)+
      scale_colour_gradient(name = "# visits",guide="legend") +
      scale_linetype_discrete(name=expression(p["sim"]))+
      scale_y_continuous(limits=c(0,1))+
      labs(x="Number of cells sampled", y="Detected trend/Number of runs")+
      facet_grid(facet.form, labeller="label_both"),
     vp=viewport(layout.pos.row=1,layout.pos.col=1)) }
    
    return(dta_sum)
}     

# Uses loess smoother to interpolate number of cells required for power
find.power<-function(dta, pwr=.8, n_grid=NULL){
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

look.up<-function(df){
  use<-subset(key,N==df$N[1] & lmda==df$lmda[1] & grid_size==df$grid_size[1])
  if(nrow(use)>1){ use<-use[1,] }
  if(nrow(use)==1 & !is.na(use$file)){
    dta<-load_data(use)
    dta<-subset(dta, n_visits==5 & gap_yr==F)
    tmp<-ddply(dta, .(detP), find.power)
    tmp<-data.frame(tmp, Test=use$test)
    return(tmp)}} 
#########################
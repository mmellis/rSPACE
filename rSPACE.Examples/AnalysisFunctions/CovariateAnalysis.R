CovariateAnalysis<-function(n_yrs, ch=NULL, n_visit=NULL, sample_yr=0, FPC=1, grdID, ...){

## Set up output when ch=NULL
  sim_results<-data.frame(p_est=0, trend=0, trendSE=0, singular=0, matrix(0,1,n_yrs))
  if(is.null(ch)) return(sim_results)

## If ch has a value, do the rest of it...


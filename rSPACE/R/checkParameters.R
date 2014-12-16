# Check parameter list for missing values, update names, etc

checkParameters<-function(pList,argList){

  if(is.null(pList$detP))
    pList$detP <- 1

  if(is.null(pList$maxDistQ))
    pList$maxDistQ <- rep(1, length(pList$MFratio))

  if(is.null(pList$wghts))
    pList$wghts <- T

  if(is.null(pList$filter.cutoff)){
    if(!is.null(argList$filter.map)){
      if(all(getValues(argList$filter.map)<=1)){
        pList$filter.cutoff <- 0.95
      }
    }
  }

  if('trunk' %in% names(pList))
    stop("Truncation parameter has been reworked to use quantiles instead of SDs.
      Use 'maxDistQ' instead of 'trunk'")

  oldnames <- c('howfar', 'howmuch', 'trunk', 'HRcenter.cutoff')
  newnames <- c('moveDist', 'moveDistQ', 'maxDistQ', 'habitat.cutoff')

  if(any(oldnames %in% names(pList))){
    oldnames<-match(oldnames, names(pList))
    names(pList)[oldnames] <- newnames
    message('Parameter names have been updated')
  }

  return(pList)
}
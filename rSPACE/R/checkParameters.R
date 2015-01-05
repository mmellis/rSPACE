# Check map for NAs, NaNs, etc
checkMap<-function(map){
  if(!grepl('proj=utm|proj=longlat',proj4string(map))) stop("Projection needs to be in utm or longlat")

  if(grepl('+proj=utm.*',proj4string(map)))
    if(!grepl('+units=m',proj4string(map)))
      message('Assuming UTM +units=m')

  map <- reclassify(map, cbind(c(NA), c(0)))

  if(any(is.nan(getValues(map)))) stop('NaNs in habitat map')

  return(map)
}


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
    stop("Truncation parameter has been reworked using probabilities instead of SDs.
      Use 'maxDistQ' instead of 'trunk'")

  oldnames <- c('howfar', 'howmuch', 'trunk', 'HRcenter.cutoff')
  newnames <- c('moveDist', 'moveDistQ', 'maxDistQ', 'habitat.cutoff')

  if(any(oldnames %in% names(pList))){
    old<-match(oldnames, names(pList))
    names(pList)[old] <- newnames[!is.na(old)]
    warning(paste('The following parameter names are deprecated:', 
                    oldnames[!is.na(old)],
                  '\n  Replace with: ', newnames[!is.na(old)]))
  }

  return(pList)
}
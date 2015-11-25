testPopSize<-function(Nst, Nen, length.out=10, map, pList){
   if(Nst < length(pList$MFratio)) stop('Cannot fit population with fewer than nTypes of individuals')
   
   pList$trendtype='abundance-linear'
   pList$N = Nst
   pList$n_yrs= length.out
   pList$lmda = (Nen/Nst)^(1/(length.out-1))

   Pop <- encounter.history(map, pList, printN=0)
   return(Pop[,-c(1:2)])
}
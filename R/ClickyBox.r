## tclck2 package - testing for wolverine
enter.parameters<-function(Parameters=NULL){
require(tcltk2)
tt <- tktoplevel()            # Create window
tkwm.title(tt,"Parameters")   # Title

done <- tclVar(0)
if(is.null(Parameters)){
  Parameters=list(N=250,
         lmda=0.933,
         n_yrs=10,
         n_visits=5,
         grid_size=25,
         sample.cutoff=0.5,
         HRcenter.cutoff=0.5,
         buffer = 1.5*c(2.5,4.2),
         howmuch = c(0.9,0.9),
         howfar = c(2.5,4.2),
         trunk = c(1,1),
         wghts=F,
         MFratio = c(0.6,0.4),
         detP=1)}


tkgrid(tklabel(tt,text=" "))

tkgrid(tklabel(tt,text="-------Population parameters-------"))
N.Val <- tclVar(Parameters$N)
enterN <-tk2entry(tt,textvariable=N.Val, width=5)
tkgrid(tklabel(tt,text="Initial population size? "),enterN)

lmda.val<-tclVar(Parameters$lmda)
enterL <-tk2entry(tt,textvariable=lmda.val, width=5)
tkgrid(tklabel(tt,text="Population growth rate? "),enterL)

grps.val<-tclVar(length(Parameters$MFratio))
enterG <-tk2entry(tt,textvariable=grps.val, width=2)
tkgrid(tklabel(tt,text="Number of separate groups? "),enterG)

MF.val<-tclVar(Parameters$MFratio)
enterMF<-tk2entry(tt,textvariable=MF.val, width=15)
tkgrid(tklabel(tt,text="Proportion of population by group? "),enterMF)

tkgrid(tklabel(tt,text=" "))
tkgrid(tklabel(tt,text="-----Territoriality/movement parameters-----"))

buff.val<-tclVar(Parameters$buffer)
enterBuff<-tk2entry(tt,textvariable=buff.val, width=15)
tkgrid(tklabel(tt,text="Buffer distance between homerange centers (km)? "),enterBuff)

HRcutoff.val<-tclVar(Parameters$HRcenter.cutoff)
enterHRcut<-tk2entry(tt,textvariable=HRcutoff.val, width=15)
tkgrid(tklabel(tt,text="Minimum habitat value for home range centers? "),enterHRcut)

#cb <- tkcheckbutton(tt)
#cbValue <- tclVar("0")
#tkconfigure(cb,variable=cbValue)
#tkgrid(tklabel(tt,text="Map provides use probabilities/weights? "),cb)
#
how.far.val<-tclVar(Parameters$howfar)
enterHF<-tk2entry(tt,textvariable=how.far.val, width=15)
tkgrid(tklabel(tt,text="Movement radius (km)? "),enterHF)

how.much.val<-tclVar(Parameters$howmuch)
enterHM<-tk2entry(tt,textvariable=how.much.val, width=15)
tkgrid(tklabel(tt,text="Percent of time spent in movement radius? "),enterHM)

trunk.val<-tclVar(Parameters$trunk)
enterTrunk<-tk2entry(tt,textvariable=trunk.val, width=15)
tkgrid(tklabel(tt,text="Truncation to movement (# SDs, 0 = no truncation) "),enterTrunk)


tkgrid(tklabel(tt,text=" "))
tkgrid(tklabel(tt,text="---------------Sampling Effort---------------"))

yrs.val<-tclVar(Parameters$n_yrs)
enter.yrs <-tk2entry(tt,textvariable=yrs.val, width=5)
tkgrid(tklabel(tt,text="Maximum number of years? "),enter.yrs)

visit.val <- tclVar(Parameters$n_visits)
enter.visit <-tk2entry(tt,textvariable=visit.val, width=5)
tkgrid(tklabel(tt,text="Maximum visits per year? "),enter.visit)

grid.val<-tclVar(Parameters$grid_size)
enter.grid <-tk2entry(tt,textvariable=grid.val, width=5)
tkgrid(tklabel(tt,text="Cell size (in km\u00b2)? "),enter.grid)

gridcut.val<-tclVar(Parameters$sample.cutoff)
enter.grid2 <-tk2entry(tt,textvariable=gridcut.val, width=5)
tkgrid(tklabel(tt,text="Percent of cell in habitat? "),enter.grid2)

tkgrid(tklabel(tt,text=" "))

# Subfunction #
check.values<-function(Parameters,n.grps){
  ok.check<-1
    
  if(length(Parameters$buffer) < n.grps) {ok.check <- 0}
  if(length(Parameters$trunk) < n.grps) {ok.check <- 0}
  if(length(Parameters$howmuch) < n.grps) {ok.check <- 0}   
  if(length(Parameters$howfar) < n.grps) {ok.check <- 0}  
  
  
  if(ok.check == 0) {cat("Something's wrong!\n");flush.console()}
  return(ok.check)
  } 
##

# Create two buttons to set the value of done
OK.but <- tkbutton(tt,text="  OK  ",    command=function() tclvalue(done)<-1)
Cancel.but <- tkbutton(tt,text="Cancel",command=function() tclvalue(done)<-2)
tkgrid(OK.but)
tkgrid(Cancel.but)

tkbind(tt,"<Destroy>",function() tclvalue(done)<-2)
tkwait.variable(done)

doneVal <- as.integer(tclvalue(done))


if(doneVal==1){
	Parameters=list(N=as.numeric(tclvalue(N.Val)),
         lmda=as.numeric(tclvalue(lmda.val)),
         n_yrs=as.numeric(tclvalue(yrs.val)),
         n_visits=as.numeric(tclvalue(visit.val)),
         grid_size=as.numeric(tclvalue(grid.val)),
         sample.cutoff=as.numeric(tclvalue(gridcut.val)),
         HRcenter.cutoff=as.numeric(tclvalue(HRcutoff.val)),
         buffer = as.numeric(unlist(strsplit(tclvalue(buff.val),split=' '))),
         howmuch = as.numeric(unlist(strsplit(tclvalue(how.much.val),split=' '))),
         howfar = as.numeric(unlist(strsplit(tclvalue(how.far.val),split=' '))),
         trunk = as.numeric(unlist(strsplit(tclvalue(trunk.val),split=' '))),
         MFratio = as.numeric(unlist(strsplit(tclvalue(MF.val),split=' '))),
         detP=1)
         
  n.grps<-as.numeric(tclvalue(grps.val))
  check.values(Parameters,n.grps)
}   

tkdestroy(tt)
detach(package:tcltk2)


return(Parameters)
}

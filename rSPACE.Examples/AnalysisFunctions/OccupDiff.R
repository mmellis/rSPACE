################################################################################
## Analysis function for rSPACE to run a comparison of single season occupancy models
##  Attempt to replicate power analysis from Guillera-Arroita and Lahoz-Monfort 2012,
##  Methods in Ecology and Evolution
##
## This approach treats the two occupancy estimates as independent samples.  In
##  rSPACE simulations sampling will reflect Markovian dependence, with occupancy
##  in cells positively correlated with previous occupancy status due to site fidelity.
##  Using separate single season occupancy models to test for a difference thus
##  represents a conservative test of the difference in estimates, since the
##  positive covariation between occupancy at t1 and t2 will reduce variance in
##  the estimate of the difference.
##
## This implementation of the analysis function for rSPACE uses the alternative
##  model formulations to looking how many years between sampling are necessary
##  to detect a constant decline.  Thus, alt_model in the subsetting options
##  should reflect the number of years between the first and second occupancy
##  estimates.  However, it could also be used to test for a fixed difference
##  in abundance between two samples by setting n_yrs=2 and lmda = 1-DIFF in the
##  parameters for createReplicates().  In this case, use M=0 only.
##
##  M. Ellis, 2/20/2015, rSPACE v1.0.7
################################################################################

OccupDiff<-function(n_yrs, ch=NULL, n_visit, M, FPC, ...){

  DF<-data.frame(trend = NA, trendSE=NA, p1=NA, X1=NA, X1se=NA, p2=NA, X2=NA, X2se=NA)
  if(is.null(ch)) return(DF)

  ConvertToMatrix<-function(ch){
    ch<-strsplit(ch, split='')
    ch<-lapply(ch, function(x) as.numeric(x))
    ch<-do.call(rbind, ch)
    return(ch)}

  stM<-n_visit * (M-1) + 1
  enM<-n_visit * M
  # Use alt_model specifications to select which year to use
  h1<-ConvertToMatrix(substr(ch, 1, n_visit))
  h2<-ConvertToMatrix(substr(ch, stM, enM))


 # From Guillera-Arroita and Lahoz-Monfort 2012, Appendix S3b (MEE3_225_sm_appendixS3b.R)
 #  Assuming independence between samples as the more conservative option

    loglik <- function(param, h)
    {
      s   <- dim(h)[1] # nr sites
      k   <- dim(h)[2] # nr sampling occasions
      psi <- 1/(1+exp(-param[1]))  # to probabilities
      p   <- 1/(1+exp(-param[2]))
      d  	<- sum(sum(h)) # summary statistics
      Sd  <- sum(rowSums(h)>0)
      loglik <- -(Sd*log(psi)+d*log(p)+(k*Sd-d)*log(1-p)+(s-Sd)*log((1-psi)+psi*(1-p)^k))
      return(loglik)
    }

    ## Analyze two datasets with different psi and different p
    fitmA <-function(h1,h2)
    {
      fm1 <- optim(par=runif(2), fn=loglik, h=h1, hessian=T)
      fm2 <- optim(par=runif(2), fn=loglik, h=h2, hessian=T)
      VC1 <- try(solve(fm1$hessian),silent = TRUE)
      VC2 <- try(solve(fm2$hessian),silent = TRUE)
      pars <- 1/(1+exp(-c(fm1$par,fm2$par)))    # to probabilities
      if ((class(VC1)=="try-error")||(class(VC2)=="try-error"))
      {
        SEs <- rep(NA,4)
      }else{
        SEs <- c(sqrt(diag(VC1))*pars[1:2]*(1-pars[1:2]),sqrt(diag(VC2))*pars[3:4]*(1-pars[3:4]))
      }
      psi1  <- list(est = pars[1], se = SEs[1]) # arrange output
      p1    <- list(est = pars[2], se = SEs[2])
      psi2  <- list(est = pars[3], se = SEs[3])
      p2    <- list(est = pars[4], se = SEs[4])
      myres <- list(psi1 = psi1, p1 = p1, psi2 = psi2, p2 = p2,L = fm1$value + fm2$value)
      return(myres)
    }
  #
  # Back from G-A/L-M appendix to rSPACE implementation...

   fits<-fitmA(h1,h2)

   DF$p1<-fits$p1
   DF$p2<-fits$p2
   DF$X1<-fits$psi1$est
   DF$X1se<-fits$psi1$se
   DF$X2<-fits$psi2$est
   DF$X2se<-fits$psi2$se

   DF$trend  <- DF$X2 - DF$X1
   DF$trendSE<- sqrt(DF$X2se^2 + DF$X1se^2)

   return(DF)
   }


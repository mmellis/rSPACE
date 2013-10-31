### Main function determining what MARK model is run, etc.
       
test_file<-function(n_yrs, ch=NULL, n_visit=NULL, gap_yr=0, FPC=1, sample_matrix=NULL, xxx=1){
  
  sim_results<-data.frame(p_est=0, trend=0, trendSE=0, matrix(0,1,n_yrs))
  if(is.null(ch)) return(sim_results)

  if(gap_yr == 0){
          mark_data<-data.frame(ch=ch,freq=rep(1,length(ch)),stringsAsFactors=F)
          test_processed=process.data(mark_data,model="RDOccupEG",time.intervals=time_int(n_visit,n_yrs))
          test_ddl=make.design.data(test_processed)
            test_ddl$Epsilon$eps=-1
            test_ddl$Gamma$eps=1
  
            Psi.t=list(formula=~time)
            p.session=list(formula=~session)
            Epsilon.random.shared=list(formula=~-1+eps:time, share=TRUE)
          model.parameters=list(Epsilon=Epsilon.random.shared,p=p.session)

          RDoccupancy=tryN(mark(test_processed,test_ddl,model.parameters=model.parameters,delete=T,output=F,prefix=paste("mark",LETTERS[xxx],sep='')))
                             
          derived_psi <- tryN(RDoccupancy$results$derived[,1])
          derived_psi_vcv <-tryN(RDoccupancy$results$derived.vcv)
          P_est <-tryN(RDoccupancy$results$real$estimate[which(row.names(RDoccupancy$results$real)=="p g1 s1 t1")])
   
          Trend_DM=cbind(1,1:n_yrs)
          Random.effects.model=tryN(variance.components(derived_psi,Trend_DM,derived_psi_vcv, REML=T))
                
  } else if(gap_yr == 1) {    #runs gap year analysis

          ch_gap<-drop_years(ch, n_visit)
          mark_data<-data.frame(ch=ch_gap,freq=rep(1,length(ch_gap)),stringsAsFactors=F)
          test_processed=process.data(mark_data,model="RDOccupEG",time.intervals=time_int(n_visit,n_yrs-1))
          test_ddl=make.design.data(test_processed)
            gap_matrix = matrix(rep(c(1,1,rep(0,n_yrs-2)),length.out=4*(n_yrs-2)),n_yrs-2,4)
            test_ddl$Epsilon$gap = -1*gap_matrix
            test_ddl$Gamma$gap = gap_matrix
  
            p.even.indices = which(as.numeric(test_ddl$p$session) %% 2 == 0)
            p.fixed.values = rep(0,length(p.even.indices))

            p.session.fixed=list(formula=~session,fixed=list(index=p.even.indices,value=p.fixed.values))
            Epsilon.gapsampling=list(formula=~-1+gap, share=TRUE)
          model.parameters=list(Epsilon=Epsilon.gapsampling,p=p.session.fixed)

          RDoccupancy=tryN(mark(test_processed,test_ddl,model.parameters=model.parameters,delete=T,output=F,prefix=paste("mark",LETTERS[xxx],sep='')))
          derived_psi <- tryN(RDoccupancy$results$derived[,1])
          derived_psi_vcv <-tryN(RDoccupancy$results$derived.vcv)
          P_est <-tryN(RDoccupancy$results$real$estimate[which(row.names(RDoccupancy$results$real)=="p g1 s1 t1")])
          
          if(!is.null(RDoccupancy)){
            derived_psi <-c(derived_psi,NA)
            new_derived_psi<-derived_psi[seq(1,n_yrs-1,by=2)]
            new_derived_psi_vcv<-derived_psi_vcv[seq(1,n_yrs-1,by=2),seq(1,n_yrs-1,by=2)]
          }
          
          Trend_DM = cbind(1,seq(1,n_yrs-1,by=2))
          Random.effects.model=tryN(variance.components(new_derived_psi,Trend_DM,new_derived_psi_vcv, REML=T))            
      
  } else if(gap_yr == 2) {    #runs skip year analysis
          ch<-drop_years(ch, n_visit, samples=sample_matrix)
          mark_data<-data.frame(ch=ch,freq=rep(1,length(ch)),stringsAsFactors=F)
          test_processed=process.data(mark_data,model="RDOccupEG",time.intervals=time_int(n_visit,n_yrs))
          test_ddl=make.design.data(test_processed)
            test_ddl$Epsilon$eps=-1
            test_ddl$Gamma$eps=1
  
            Psi.t=list(formula=~time)
            p.session=list(formula=~session)
            Epsilon.random.shared=list(formula=~-1+eps:time, share=TRUE)
          model.parameters=list(Epsilon=Epsilon.random.shared,p=p.session)

          RDoccupancy=tryN(mark(test_processed,test_ddl,model.parameters=model.parameters,delete=T,output=F,prefix=paste("mark",LETTERS[xxx],sep='')))
                             
          derived_psi <- tryN(RDoccupancy$results$derived[,1])
          derived_psi_vcv <-tryN(RDoccupancy$results$derived.vcv)
          P_est <-tryN(RDoccupancy$results$real$estimate[which(row.names(RDoccupancy$results$real)=="p g1 s1 t1")])
   
          Trend_DM=cbind(1,1:n_yrs)
          Random.effects.model=tryN(variance.components(derived_psi,Trend_DM,derived_psi_vcv, REML=T)) 
   }
   
   if(!is.null(Random.effects.model)){
     sim_results$trend          <- Random.effects.model$beta[2,1]
     sim_results$trendSE        <- FPC_trendSE(Random.effects.model, k=nrow(Trend_DM), FPC)
   }
   if(!is.null(derived_psi)){
     sim_results[1,4:(n_yrs+3)] <- matrix(derived_psi,nrow=1)
     sim_results$p_est          <- P_est
   }
   
   return(sim_results)
}
   

            
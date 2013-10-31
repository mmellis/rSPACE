build.key<-function(species.path){
  tests<-list.dirs(path=species.path)
  tests<-tests[grepl('test',tests)]
  parameters_hold<-list()
  for(ii in 1:length(tests)){
    if('Parameters.txt' %in% dir(tests[ii]))
      {
        source(file=dir(tests[ii], pattern="Parameters.txt", full.names=T))
        Parameters<-list(n_yrs=n_yrs, 
            n_visits=n_visits, 
            N=N, 
            lmda=lmda, 
            MFratio=MFratio, 
            buffer=buffer, 
            howmuch=howmuch, 
            howfar=howfar, 
            trunk=trunk, 
            grid_size=grid_size, 
            HRcenter.cutoff=cutoff)
        parameters_hold[[ii]]<-Parameters
      } else if('Parameters.rdata' %in% dir(tests[ii]))
      {
        load(paste(tests[ii], 'Parameters.rdata', sep=''))
        parameters_hold[[ii]]<-Parameters
      }}
      return(parameters_hold)
      }
      
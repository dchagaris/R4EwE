#plotting----
fn.ecospace_plot_ts <- function(predB=predB, predC=predC, timestep='annual',obs.ts=obs.ts, scale2run=1, pltB.dims=c(1,1), pltC.dims=c(1,1), scaleCatch=FALSE,
                                dir.plts = dir.pred, plot2pdf=TRUE){
  xtime = as.numeric(dimnames(predB)[[1]])
  if(plot2pdf) pdf(file.path(dir.plts,"biomass timeseries fits.pdf"), onefile=T)
  par(mfrow=pltB.dims, mar=c(2,4,2,1))
  for(s in 1:dim(predB)[2]){
    #s=6
    has.obs = ifelse(s%in%obs.ts$obsB.head$Poolcode,TRUE,FALSE)
    pred.s = predB[,s,]
    pred.base.s = predB[,s,scale2run]
    obs.s=NULL
    plt.ylims = c(0,max(pred.s,na.rm=T)*1.2)
    if(has.obs){
      obs.s = as.matrix(obs.ts$obsB[,which(obs.ts$obsB.head$Poolcode==s)])
      colnames(obs.s) = obs.ts$obsB.head$Title[which(obs.ts$obsB.head$Poolcode==s)]
      #rescale obs to EwE units
      obs.q = colMeans(obs.s, na.rm=T)/mean(pred.base.s, na.rm=T)  
      obs.scaled.s = sweep(obs.s,2,obs.q,"/")
      plt.ylims = c(0,max(pred.s,obs.scaled.s,na.rm=T)*1.2)
    }
    matplot(xtime,pred.s, type='l', lty=1, lwd=2, ylim=plt.ylims, main=dimnames(predB)[[2]][s])
    if(has.obs){
      legend('topleft',legend=colnames(obs.scaled.s), pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=.8)    
      matpoints(xtime,obs.scaled.s, type='p', pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=0.75)
    } 
  }
  if(plot2pdf) dev.off()
  
  if(plot2pdf) pdf(file.path(dir.plts,"catch timeseries fits.pdf"), onefile=T)
  par(mfrow=pltC.dims, mar=c(2,4,2,1))
  for(s in 1:dim(predC)[2]){
    #s=12
    pred.name = paste(unique(strsplit(dimnames(predC)[[2]][s],split="\\.")[[1]][-c(1,2)],fromLast=T),collapse="_")
    ss = match(pred.name, gsub("\\+","",gsub("-","_",df.names$group.names)))
    has.obs = ifelse(ss%in%obs.ts$obsC.head$Poolcode,TRUE,FALSE)
    pred.s = predC[,s,]
    pred.base.s = predC[,s,scale2run]
    obs.s=NULL
    plt.ylims = c(0,max(pred.s,na.rm=T)*1.2)
    if(has.obs){
      obs.s = as.matrix(obs.ts$obsC[,which(obs.ts$obsC.head$Poolcode==ss)])
      obs.type.s = obs.ts$obsC.head$Type[which(obs.ts$obsC.head$Poolcode==ss)]
      colnames(obs.s) = obs.ts$obsC.head$Title[which(obs.ts$obsC.head$Poolcode==ss)]
      #rescale obs to EwE units
      if(scaleCatch){ obs.q = colMeans(obs.s, na.rm=T)/mean(pred.base.s, na.rm=T) 
      obs.q[obs.type.s==6] = 1
      }
      if(!scaleCatch) obs.q = rep(1,ncol(obs.s))
      obs.scaled.s = sweep(obs.s,2,obs.q,"/")
      plt.ylims = c(0,max(pred.s,obs.scaled.s,na.rm=T)*1.2)
    }
    matplot(xtime,pred.s, type='l', lty=1, lwd=2, ylim=plt.ylims, main=pred.name)
    if(has.obs){
      legend('topleft',legend=colnames(obs.scaled.s), pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=.8)    
      matpoints(xtime,obs.scaled.s, type='p', pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=0.75)
    } 
  }
  if(plot2pdf) dev.off()
  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title Plot Ecospace predicted timeseries.
#' @description Makes a multipanel plot of 1 or more ecopace runs with observed data where available.
#' @param predB An array of biomass predictions, as produced by fn.ecospace_ts2array()
#' @param predC An array of catch predictions, as produced by fn.ecospace_ts2array()
#' @param obs.ts A list object containing reference timeseries information, as that created by 
#' fn.read_ecosim_timeseries().
#' @param timestep Read 'annual' or 'monthly' data.
#' @param scale2run If multiple models, which to scale the observed values to.
#' @param pltB.dims Numeric vector of length 2 specifying the number of rows and columns for multipanel plotting.
#' @param pltC.dims Numeric vector of length 2 specifying the number of rows and columns for multipanel plotting.
#' @param scaleCatch Logical. TRUE will scale the predictions so the mean is equal to that of the observed.  Default is FALSE.
#' @param plt.cols Vector of line colors.
#' @param plot2pdf Logical.  TRUE (default) will save files to the specified directory.
#' @param dir.plts  Location to save pdf plots. Defaults to dir.pred.
#' @return Generate a figure in the active plotting device or saves as pdf file.
#' @examples
#' # example code:
#' \dontrun{result <- fn.runEwE.parallel(runlist=myrunlist, obj.fxn=1, cl.export=list('myrunlist','obs.ts'))}
#' @export
fn.ecospace_plot_ts <- function(predB=predB, predC=predC, timestep='annual',obs.ts=obs.ts, scale2run=1, pltB.dims=c(3,3), pltC.dims=c(1,1), 
                                scaleCatch=FALSE,plt.cols=1:dim(predB)[3],
                                dir.plts = dir.pred, plot2pdf=TRUE, run.label=Sys.Date()){
  
  # timestep='annual'
  # scale2run=1
  # pltB.dims=c(4,3)
  # pltC.dims=c(4,3)
  # scaleCatch=TRUE
  # plt.cols=1:dim(predB)[3]
  # dir.plts = dir.pred
  # plot2pdf=FALSE
  
  xtime = as.numeric(dimnames(predB)[[1]])
  if(plot2pdf) pdf(file.path(dir.plts,paste0("biomass timeseries fits_",run.label,".pdf")), onefile=T)
  par(mfrow=pltB.dims, mar=c(2,4,2,1), oma=c(4,0,0,1),xpd=F)
  
  for(s in 1:dim(predB)[2]){
    #s=2
    par(xpd=F)
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
    matplot(xtime,pred.s, type='l', lty=1, lwd=2, ylim=plt.ylims, main=dimnames(predB)[[2]][s], col=plt.cols, xlab='',ylab='biomass')
    if(has.obs){
      legend('topleft',legend=colnames(obs.scaled.s), pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=.8)    
      matpoints(xtime,obs.scaled.s, type='p', pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=0.75)
    }
    if(s%in%c(prod(pltB.dims)*(1:ceiling(dim(predB)[2]/prod(pltB.dims))),dim(predB)[2])){
      #add legend
    }
  }
  if(plot2pdf) dev.off()
  
  
  if(plot2pdf) pdf(file.path(dir.plts,paste0("catch timeseries fits_",run.label,".pdf")), onefile=T)
  par(mfrow=pltC.dims, mar=c(2,4,2,1))
  predC.agg <- fn.agg_catch_by_group(predC)
  dimnames(predC.agg)
  for(s in 1:dim(predC.agg)[2]){
    #s=26
    #pred.name = paste(unique(strsplit(dimnames(predC)[[2]][s],split="\\.")[[1]][-c(1,2)],fromLast=T),collapse="_")
    #ss = match(pred.name, gsub("\\+","",gsub("-","_",df.names$group.names)))
    pred.name = gsub(" ","_",dimnames(predC.agg)[[2]][s])
    ss = match(pred.name, df.names$group.names)
    has.obs = ifelse(ss%in%obs.ts$obsC.head$Poolcode,TRUE,FALSE)
    pred.s = predC.agg[,s,]
    pred.base.s = predC.agg[,s,scale2run]
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
    matplot(xtime,pred.s, type='l', lty=1, lwd=2, ylim=plt.ylims, main=pred.name, col=plt.cols, xlab='',ylab='catch')
    if(has.obs){
      legend('topleft',legend=colnames(obs.scaled.s), pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=.8)    
      matpoints(xtime,obs.scaled.s, type='p', pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=0.75)
    } 
  }
  if(plot2pdf) dev.off()
  
}#eof


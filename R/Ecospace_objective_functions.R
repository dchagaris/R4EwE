#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#OBJECTIVE FUNCTIONS---------------------------------------------------------------------------------  
#Loop to run objfxn for sensitivity and fill in runlist ////////////////////////
#print(paste("Sensitivity run obj. fxn.", start, "to", stop))

#...............................................................................
##Objective fxn 1-----
#fits to timeseries data only, using the ecosim timeseries csv file
fn.objfxn1 <- function(dir.pred, obs.ts=obs.ts, run.idx=1, fit.abs.catch=TRUE){  
  #get annual timeseries predictions
  #dir.pred = "C:/NWACS MICE/GA output 2026-01-06/GA_Run_20260106_114454/run_00abeb6e03a8897940ee4aa45b593526"
  predB = fn.ecospace_predB_ts2array(dir.out = dir.pred, timestep='annual', n.reg=0)[,,run.idx]
  predC = fn.ecospace_predC_ts2array(dir.out = dir.pred, timestep='annual', n.reg=0)[,,run.idx] #unable to read region catch so far
  dim(predB); dim(predC)
  ###biomass timeseries annual----
  obs.ts.head = obs.ts$obsB.head
  obs.ts.biomass = obs.ts$obsB
  names(obs.ts.head) = c('title','weight','poolcode','type')
  lk.ts.biomass <- cbind(obs.ts.head, dattype='biomass timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=5
    #grp = group.names[j]
    grpnum = obs.ts.head$poolcode[j]
    grp = group.names[grpnum]
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    lk.dat = data.frame(pred=predB[,grpnum], obs=obs.ts.biomass[,j])
    lk.dat$obs = ifelse(lk.dat$obs==0,NA,lk.dat$obs)
    rownames(lk.dat) <- rownames(predB)
    if(obs.ts.head$type[j]==0){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$pred = lk.dat$pred/q
    lk.dat$ll = log(lk.dat$pred/lk.dat$obs)^2/(2*obs.sd^2)
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.biomass$obs.cv[j] = obs.cv
    lk.ts.biomass$obs.sd[j] = obs.sd
    lk.ts.biomass$loglik[j] = lk.sum
  }
  
  ###catch timeseries annual----
  obs.ts.head = obs.ts$obsC.head
  obs.ts.catch = obs.ts$obsC
  names(obs.ts.head) = c('title','weight','poolcode','type')
  lk.ts.catch <- cbind(obs.ts.head, dattype='catch timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=2
    #grp = group.names[j]
    grpnum = obs.ts.head$poolcode[j]
    grp = gsub("/+",".",gsub("-",".",gsub("_", ".", group.names[grpnum])))
    idx_group = which(grepl(grp,colnames(predC)))
    #idx_group = match(grp,colnames(predC))
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    if(length(idx_group)==1) lk.dat = data.frame(pred=predC[,idx_group], obs=obs.ts.catch[,j])
    if(length(idx_group)>1) lk.dat = data.frame(pred=rowSums(predC[,idx_group]), obs=obs.ts.catch[,j])
    lk.dat$obs[lk.dat$obs <=0] <- NA
    rownames(lk.dat) <- rownames(predC)
    
    if(obs.ts.head$type[j]==61 | !fit.abs.catch){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$pred = lk.dat$pred/q
    lk.dat$ll = log(lk.dat$pred/lk.dat$obs)^2/(2*obs.sd^2)
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.catch$obs.cv[j] = obs.cv
    lk.ts.catch$obs.sd[j] = obs.sd
    lk.ts.catch$loglik[j] = lk.sum
  } 
  
  lk.ts <<- rbind(lk.ts.biomass,lk.ts.catch)
  
  ## Make output vector................................
  ##
  lk.vec <- rep(0,length(group.names))
  lk.agg <- aggregate(loglik~poolcode, data=lk.ts, sum,na.rm=T)
  lk.vec[lk.agg$poolcode] <- lk.agg$loglik
  
  outvec = round(c(sum(lk.vec), lk.vec),2)
  names(outvec) = c('LL.total',paste0('LL.',group.names))
  
  #write.csv(outvec,paste0(runlist$dir.out,"/objvals.csv"))
  return(outvec)
}

fn.objfxn2 <- function(dir.pred, obs.ts=obs.ts, obs.maps=obs.maps, obs.maps.meta=obs.maps.meta, autoweight.LL=T){
  #get annual timeseries predictions
  predB = fn.ecospace_predB_ts2array(dir.out = dir.pred, timestep='annual', n.reg=0)[,,1]
  predC = fn.ecospace_predC_ts2array(dir.out = dir.pred, timestep='annual', n.reg=0)[,,1] #unable to read region catch so far
  
  ###biomass timeseries annual----
  obs.ts.head = obs.ts$obsB.head
  obs.ts.biomass = obs.ts$obsB
  names(obs.ts.head) = c('title','weight','poolcode','type')
  lk.ts.biomass <- cbind(obs.ts.head, dattype='biomass timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=1
    #grp = group.names[j]
    grpnum = obs.ts.head$poolcode[j]
    grp = group.names[grpnum]
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    lk.dat = data.frame(pred=predB[,grpnum], obs=obs.ts.biomass[,j])
    rownames(lk.dat) <- rownames(predB)
    if(obs.ts.head$type[j]==0){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$pred = lk.dat$pred/q
    lk.dat$ll = log(lk.dat$pred/lk.dat$obs)^2/(2*obs.sd^2)
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.biomass$obs.cv[j] = obs.cv
    lk.ts.biomass$obs.sd[j] = obs.sd
    lk.ts.biomass$loglik[j] = lk.sum
  }
  
  ###catch timeseries annual----
  obs.ts.head = obs.ts$obsC.head
  obs.ts.catch = obs.ts$obsC
  names(obs.ts.head) = c('title','weight','poolcode','type')
  lk.ts.catch <- cbind(obs.ts.head, dattype='catch timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=2
    #grp = group.names[j]
    grpnum = obs.ts.head$poolcode[j]
    grp = gsub("/+",".",gsub("-",".",gsub("_", ".", group.names[grpnum])))
    idx_group = which(grepl(grp,colnames(predC)))
    #idx_group = match(grp,colnames(predC))
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    lk.dat = data.frame(pred=predC[,idx_group], obs=obs.ts.catch[,j])
    lk.dat$obs[lk.dat$obs <=0] <- NA
    rownames(lk.dat) <- rownames(predC)
    
    if(obs.ts.head$type[j]==61){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$pred = lk.dat$pred/q
    lk.dat$ll = log(lk.dat$pred/lk.dat$obs)^2/(2*obs.sd^2)
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.catch$obs.cv[j] = obs.cv
    lk.ts.catch$obs.sd[j] = obs.sd
    lk.ts.catch$loglik[j] = lk.sum
  } 
  
  lk.ts = rbind(lk.ts.biomass,lk.ts.catch)
  
  #fit to raster grids----
  #probability distribution maps (maps should sum to 1)
  lk.maps = cbind(obs.maps.meta, dattype='prob dist map',obs.cv=NA, obs.sd=NA, loglik=NA)
  
  for(i in 1:nlayers(obs.maps)){
    #i=1
    obs.i = rast(obs.maps[[i]])
    spp.i = strsplit(names(obs.i),"_")[[1]][1]
    grpnum.i = as.numeric(strsplit(names(obs.i),"_")[[1]][3])
    grpname.i = gsub("_"," ",df.names$group.names[grpnum.i])
    source.i = obs.maps.meta$source[i]
    files.pred.asc = list.files(path=file.path(dir.pred,'asc'), pattern=paste0("EcospaceMapBiomass-",grpname.i), full.names=T)
    
    if(i==1){
      tsteps <- as.numeric(substr(basename(files.pred.asc),nchar(basename(files.pred.asc))-8,nchar(basename(files.pred.asc))-4)  )
      tsteps.y <- rep(styear:(styear-1+length(tsteps)/12),each=12)
      tsteps.m <- rep(1:12,length(tsteps)/12)
    }
    
    if(source.i=='NEFSC sdm rasters'){
      get.files.i = which(tsteps.y>=obs.maps.meta$yr.start[i] & tsteps.y<=obs.maps.meta$yr.end[i] &
                            tsteps.m>=obs.maps.meta$mo.start[i] & tsteps.m<=obs.maps.meta$mo.end[i] )
      stack.i = terra::rast(files.pred.asc[get.files.i])
      mean.i = terra::mean(stack.i,na.rm=T) #terra::app(stack.i, fun=mean, na.rm=T)
      prop.i <- mean.i / global(mean.i, sum, na.rm=TRUE)[1,1]
    }
    obs_vals <- values(obs.i)
    pred_vals <- values(prop.i)
    valid_idx <- !is.na(obs_vals) & !is.na(pred_vals)
    obs_probs <- obs_vals[valid_idx] / sum(obs_vals[valid_idx])
    pred_probs <- pred_vals[valid_idx] / sum(pred_vals[valid_idx])
    
    # Compute log-likelihood (cross-entropy)
    lk.maps$loglik[i] <- -sum(obs_probs * log(pred_probs))
  }
  
  
  ## combine likelihoods
  if(autoweight.LL){
    lk.ts$loglik.w = lk.ts$loglik*(1/mean(lk.ts$loglik))
    lk.maps$loglik.w = lk.maps$loglik*(1/mean(lk.maps$loglik))
  } else{
    lk.ts$loglik.w = lk.ts$loglik
    lk.maps$loglik.w = lk.maps$loglik
  }
  lk.comb = data.frame(obs.name=c(lk.ts$title, lk.maps$layername),
                       poolcode = as.numeric(c(lk.ts$poolcode, lk.maps$grp.num)),
                       dattype = c(lk.ts$dattype, lk.maps$dattype),
                       loglik = c(lk.ts$loglik, lk.maps$loglik),
                       loglik.w = c(lk.ts$loglik.w, lk.maps$loglik.w))
  
  lk.agg <- aggregate(loglik.w~poolcode, data=lk.comb, sum,na.rm=T)
  
  ## Make output vector................................
  lk.vec <- rep(0,length(group.names))
  lk.vec[lk.agg$poolcode] <- lk.agg$loglik.w
  
  outvec = round(c(sum(lk.vec), lk.vec),2)
  names(outvec) = c('LL.total',paste0('LL.',group.names))
  
  #write.csv(outvec,paste0(runlist$dir.out,"/objvals.csv"))
  return(outvec)
}
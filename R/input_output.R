library('raster')
library('terra')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#fn.read ecosim timeseries---------------------------------------------------------------------------------------
fn.read_ecosim_timeseries = function(filename){
  fnm.obs_ts = filename
  obs.ts.head = as.data.frame(t(read.csv(fnm.obs_ts,header=F,nrows=4)))
  names(obs.ts.head) = obs.ts.head[1,]; obs.ts.head = obs.ts.head[-1,]
  obs.ts.head[,2:4] = as.numeric(as.matrix(obs.ts.head[,2:4]))
  obs.ts = read.csv(fnm.obs_ts,header=F,skip=4)
  rownames(obs.ts) = obs.ts[,1]; obs.ts[,1] = NULL
  if(nrow(obs.ts.head) != ncol(obs.ts)) print('HEADER AND TIMESERIES DIMENSION DO NOT MATCH!!!')
  
  obsB.head = obs.ts.head[obs.ts.head[,4] %in% c(0,1),]
  obsB = obs.ts[,which(obs.ts.head[,4] %in% c(0,1))]
  obsC.head = obs.ts.head[obs.ts.head[,4] %in% c(6,61,-6),]
  obsC = obs.ts[,which(obs.ts.head[,4] %in% c(6,61,-6))]
  names(obsB.head) = names(obsC.head) = gsub(" ","_",names(obsB.head))
  ret = list(obsB.head,obsB,obsC.head,obsC)
  names(ret) = c('obsB.head','obsB','obsC.head','obsC')
  return(ret)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecospace output to arrays---------------------------------------------------------------------------------  
fn.ecospace_predB_ts2array = function(dir.out=dir.pred, timestep='annual',n.reg=0){
  
  #dir.out = "C:/NWACS MICE/GA output 2026-01-06/GA_Run_20260106_114454/run_00abeb6e03a8897940ee4aa45b593526"
  if(n.reg==0){
    if(timestep=='annual') files.bio = list.files(dir.out,pattern="Ecospace_Annual_Average_Biomass.csv",recursive=T,full.names = T)
    if(timestep=='monthly') files.bio = list.files(dir.out,pattern="Ecospace_Average_Biomass.csv",recursive=T,full.names = T)
    
    bio = lapply(files.bio,FUN=function(x){
      if(timestep=='annual') nskip = which(substr(readLines(x),1,4)=='Year')-1
      if(timestep=='monthly') nskip = which(substr(readLines(x),1,8)=='TimeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1)
    })
    
    bio.array = array(dim=c(dim(bio[[1]])[1],dim(bio[[1]])[2],length(files.bio)),
                      dimnames=list(rownames(bio[[1]]),names(bio[[1]]),basename(dirname(files.bio))))
    for(r in 1:length(bio)){
      tmp = as.matrix(bio[[r]])
      bio.array[,,r] <- tmp
    }
  }
  
  if(n.reg>0){
    if(timestep=='annual') files.bio = list.files(dir.out,pattern="^Ecospace_Annual_Average_Region.*\\Biomass.csv$",recursive=T,full.names = T)
    if(timestep=='monthly')  files.bio = list.files(dir.out,pattern="^Ecospace_Average_Region.*\\Biomass.csv$",recursive=T,full.names = T)
    
    bio = lapply(files.bio,FUN=function(x){
      if(timestep=='annual') nskip = which(substr(readLines(x),1,4)=='Year')-1
      if(timestep=='monthly') nskip = which(substr(readLines(x),1,8)=='TimeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1)
    })
    
    bio.array = array(dim=c(dim(bio[[1]])[1],dim(bio[[1]])[2],n.reg+1,length(unique(dirname(files.bio)))),
                      dimnames=list(rownames(bio[[1]]),names(bio[[1]]),paste0('reg',0:n.reg),unique(basename(dirname(files.bio)))))
    for(r in 1:length(bio)){
      if(timestep=='annual') reg.idx=as.numeric(substr(basename(files.bio[r]),32,32))+1
      if(timestep=='monthly') reg.idx=as.numeric(substr(basename(files.bio[r]),25,25))+1
      run.idx=match(basename(dirname(files.bio[r])),dimnames(bio.array)[[4]])
      tmp = as.matrix(bio[[r]])
      bio.array[,,reg.idx,run.idx] <- tmp
    }
  }
  if(timestep=='annual') dimnames(bio.array)[[1]] <- styear:(styear+dim(bio.array)[1]-1)
  if(timestep=='monthly') dimnames(bio.array)[[1]] <- length(seq(styear,(styear+dim(bio.array)[1]-1/12),1/12))
  return(bio.array)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecospace output to arrays---------------------------------------------------------------------------------  
fn.ecospace_predC_ts2array = function(dir.out=dir.pred, timestep='annual',n.reg=0){
  if(n.reg==0){
    if(timestep=='annual') files.cat = list.files(dir.out,pattern="Ecospace_Annual_Average_Catch.csv",recursive=T,full.names = T)
    if(timestep=='monthly') files.cat = list.files(dir.out,pattern="Ecospace_Average_Catch.csv",recursive=T,full.names = T)
    
    cat = lapply(files.cat,FUN=function(x){
      if(timestep=='annual') nskip = which(substr(readLines(x),1,4)=='Year')-1
      if(timestep=='monthly') nskip = which(substr(readLines(x),1,8)=='TimeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1, check.names=F)
    })
    sppnames = sapply(strsplit(names(cat[[1]]),"\\|"),tail,1)
    cat.array = array(dim=c(dim(cat[[1]])[1],dim(cat[[1]])[2],length(files.cat)),
                      dimnames=list(rownames(cat[[1]]),names(cat[[1]]),basename(dirname(files.cat))))
    for(r in 1:length(cat)){
      tmp = as.matrix(cat[[r]])
      cat.array[,,r] <- tmp
    }
  }
  
  if(n.reg>0){
    if(timestep=='annual') files.cat = list.files(dir.out,pattern="^Ecospace_Annual_Average_Region.*\\Catch.csv$",recursive=T,full.names = T)
    if(timestep=='monthly')  files.cat = list.files(dir.out,pattern="^Ecospace_Average_Region.*\\Catch.csv$",recursive=T,full.names = T)
    
    cat = lapply(files.cat,FUN=function(x){
      if(timestep=='annual') nskip = which(substr(readLines(x),1,4)=='Year')-1
      if(timestep=='monthly') nskip = which(substr(readLines(x),1,8)=='TimeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1, check.names=F)
    })
    
    cat.array = array(dim=c(dim(cat[[1]])[1],dim(cat[[1]])[2],n.reg+1,length(unique(dirname(files.cat)))),
                      dimnames=list(rownames(cat[[1]]),names(cat[[1]]),paste0('reg',0:n.reg),unique(basename(dirname(files.cat)))))
    for(r in 1:length(cat)){
      if(timestep=='annual') reg.idx=as.numeric(substr(basename(files.cat[r]),32,32))+1
      if(timestep=='monthly') reg.idx=as.numeric(substr(basename(files.cat[r]),25,25))+1
      run.idx=match(basename(dirname(files.cat[r])),dimnames(cat.array)[[4]])
      tmp = as.matrix(cat[[r]])
      cat.array[,,reg.idx,run.idx] <- tmp
    }
  }
  if(timestep=='annual') dimnames(cat.array)[[1]] <- styear:(styear+dim(cat.array)[1]-1)
  if(timestep=='monthly') dimnames(cat.array)[[1]] <- length(seq(styear,(styear+dim(cat.array)[1]-1/12),1/12))
  return(cat.array)
}#eof

#aggregate catch by group
fn.agg_catch_by_group <- function(predC){
  print(dimnames(predC)[[2]])
  grps = sapply(strsplit(dimnames(predC)[[2]],"\\|"),tail,1)
  print(grps)
  length(unique(grps))
  
  #now aggregate columns in predC.agg to grps
  predC.agg <- apply(predC, c(1,3), function(v) {tapply(v, grps, sum)})
  predC.agg <- aperm(predC.agg, c(2,1,3))
  dimnames(predC.agg)
  return(predC.agg)
}#eof

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
#fn.vul_testval_tags-----
#' @title Read Ecospace output ascii files into a raster stack.
#' @description This function reads a subset of asc grids produced by an Ecospace run into a raster stack.  Requires the terra package for faster processing.
#' @param dir.out The run output folder containing.  It should contain a folder named 'asc', which holds all the output maps.
#' @param do.vars The model output variables to read. Must match the naming convention of the output files, all  lowercase.  Current options are 'biomass','catch', and 'effort'.
#' @param do.grps The model groups for which to read data.
#' @param do.fleets The fleet numbers for which to read effort data.
#' @param do.months The months to read output data, only used if monthly maps are output.
#' @return A SpatRaster stack from the terra package.
#' @examples
#' # example code:
#' predB.maps <- fn.ecospace_ascii2stack(dir.out=dir.pred, do.vars=c('biomass'), do.grps=c(3,5,6,8,10,12), do.fleets=1:6, do.months=1, do.years=2010:2019, annualmaps=TRUE)
#' @export
fn.ecospace_ascii2stack <- function(dir.out=dir.pred,   
                                    do.vars=c('biomass'), 
                                    do.grps=1:12, 
                                    do.fleets=1:6, 
                                    do.months=1, 
                                    do.years=2010:2019,
                                    annualmaps=TRUE){
  # dir.out=dir.pred
  # do.vars=c('biomass')
  # do.grps=1:12
  # do.fleets=1:6
  # do.months=1
  # do.years=2010:2019
  # annualmaps=TRUE
  library('terra')
  files.ascii <- list.files(file.path(dir.out,'asc'),pattern=".asc$", recursive=T, full.names=T)
  grpsplit = strsplit(basename(files.ascii),"-")
  grpnames = character()
  for(g in 1:length(grpsplit)){
    #g=1326
    name.g = grpsplit[[g]]
    name.g = name.g[-c(1,length(name.g))]
    name.g = paste(name.g,collapse="-")
    grpnames = c(grpnames,name.g)
  }
  unique(grpnames)
  timestep = substr(basename(files.ascii),nchar(basename(files.ascii))-8,nchar(basename(files.ascii))-4)
  time = styear + (as.numeric(timestep)-1)/12
  #output map metadata----
  ascii.df = data.frame(type=ifelse(grepl('Biomass',basename(files.ascii)),'biomass',
                                    ifelse(grepl('Catch',basename(files.ascii)),'catch',
                                           ifelse(grepl('Effort',basename(files.ascii)),'effort','other'))),
                        group=grpnames,
                        group.num = match(gsub(" ","_",grpnames),df.names$group.names),
                        fleet.num = match(grpnames,df.names$fleet.names),
                        timestep =timestep,
                        time = time,
                        year = floor(time),
                        month = round(1+((time-styear)-floor(time-styear))*12,0),
                        file=files.ascii)
  #subset----
  if(annualmaps) ascii.df.sub = subset(ascii.df, type%in%do.vars & (group.num%in%do.grps | fleet.num%in%do.fleets) & year%in%do.years)
  if(!annualmaps) ascii.df.sub = subset(ascii.df, type%in%do.vars & (group.num%in%do.grps | fleet.num%in%do.fleets) & year%in%do.years & month%in%do.months)
  ascii.df.sub$label = ifelse(ascii.df.sub$type=='effort',
                              paste(ascii.df.sub$type,formatC(ascii.df.sub$fleet.num,width=2,flag="0"),ascii.df.sub$year,formatC(ascii.df.sub$month,width=2,flag="0"),sep="_"),
                              paste(ascii.df.sub$type,formatC(ascii.df.sub$group.num,width=2,flag="0"),ascii.df.sub$year,formatC(ascii.df.sub$month,width=2,flag="0"),sep="_"))
  
  #stack and return----
  #out.stack <- stack(ascii.df.sub$file)
  out.stack <- rast(ascii.df.sub$file)
  names(out.stack) <- ascii.df.sub$label
  return(out.stack)
}

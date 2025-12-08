library('fields')
library('colorRamps')
library('ncdf4')
library('jsonlite')
library('terra')
library('raster')

fn.pull_glorys <- function(dir.out=dir.glorys, 
                           bbox=bbox,
                           startdate = startdate,
                           enddate = enddate,
                           vars.phy = c('thetao','so','bottomT'),
                           vars.bgc = c('chl','nppv','o2','phyc')){
  
  varstring.phy <- paste(paste0("--variable ", vars.phy), collapse = " ")
  varstring.bgc <- paste(paste0("--variable ", vars.bgc), collapse = " ")
  
  #PHYSICAL REANALYSIS DATA-------------------------------------------------------
  ##get times----
  command <- paste(shQuote(path.copernicusmarine),
                   "describe --dataset-id cmems_mod_glo_phy_my_0.083deg_P1M-m --contains time --return-fields coordinates,maximum_value,minimum_value")
  raw.output = system(command, intern=T)
  json.start = which(startsWith(raw.output,"{"))
  json.output = raw.output[json.start:length(raw.output)]
  json.str = paste(json.output, collapse = "\n")
  json.data = fromJSON(json.str)
  basetime = json.data$products$datasets[[1]]$versions[[1]]$parts[[1]]$services[[1]]$variables[[1]]$coordinates[[1]]$coordinate_unit
  basetime = as.Date(gsub("[A-Za-z()]", "", basetime))
  time.millisec = range(json.data$products$datasets[[1]]$versions[[1]]$parts[[1]]$services[[1]]$variables[[1]]$coordinates[[1]]$values[[1]])
  time.range = as.Date(as.POSIXct(time.millisec/1000, origin=basetime, tz="UTC"))
  
  ##query the data----
  if(as.Date(startdate)<time.range[2]){
    startdate1 = max(as.Date(startdate),time.range[1])
    enddate1 = min(as.Date(enddate),time.range[2])
    command <- paste(
      shQuote(path.copernicusmarine),
      "subset",
      "--dataset-id cmems_mod_glo_phy_my_0.083deg_P1M-m",
      "--dataset-version 202311",
      varstring.phy,
      "--start-datetime ", paste0(startdate1,"T00:00:00"),
      "--end-datetime ",  paste0(enddate1,"T00:00:00"),
      "--minimum-longitude",bbox[1],
      "--maximum-longitude",bbox[2],
      "--minimum-latitude",bbox[3],
      "--maximum-latitude",bbox[4],
      "--coordinates-selection-method strict-inside",
      "--netcdf-compression-level 1",
      "-o",shQuote(dir.out)
    )
    system(command)
  }
  
  #BIOGEOCHEMICAL HINDCAST----------------------------------------------------------------------
  ##get times----
  command <- paste(shQuote(path.copernicusmarine),
                   "describe --dataset-id cmems_mod_glo_bgc_my_0.25deg_P1M-m --contains time --return-fields coordinates,maximum_value,minimum_value")
  raw.output = system(command, intern=T)
  json.start = which(startsWith(raw.output,"{"))
  json.output = raw.output[json.start:length(raw.output)]
  json.str = paste(json.output, collapse = "\n")
  json.data = fromJSON(json.str)
  basetime = json.data$products$datasets[[1]]$versions[[1]]$parts[[1]]$services[[1]]$variables[[1]]$coordinates[[1]]$coordinate_unit
  basetime = as.Date(gsub("[A-Za-z()]", "", basetime))
  time.millisec = range(json.data$products$datasets[[1]]$versions[[1]]$parts[[1]]$services[[1]]$variables[[1]]$coordinates[[1]]$values[[1]])
  time.range = as.Date(as.POSIXct(time.millisec/1000, origin=basetime, tz="UTC"))
  time.range
  ##query the data----
  if(as.Date(startdate)<time.range[2]){
    startdate1 = max(as.Date(startdate),time.range[1])
    enddate1 = min(as.Date(enddate),time.range[2])
    
    command <- paste(
      shQuote(path.copernicusmarine),
      "subset",
      "--dataset-id cmems_mod_glo_bgc_my_0.25deg_P1M-m",
      "--dataset-version 202406",
      varstring.bgc,
      "--start-datetime ", paste0(startdate1,"T00:00:00"),
      "--end-datetime ",  paste0(enddate1,"T00:00:00"),
      "--minimum-longitude",bbox[1],
      "--maximum-longitude",bbox[2],
      "--minimum-latitude",bbox[3],
      "--maximum-latitude",bbox[4],
      "--coordinates-selection-method strict-inside",
      "--netcdf-compression-level 1",
      "-o",shQuote(dir.out)
    )
    system(command)
  }
  
  
  #GLORYS bathymetry--------------------------------------------------------------
  ##5 min grid for phy----
  command <- paste(
    shQuote(path.copernicusmarine),
    "subset",
    "--dataset-id cmems_mod_glo_phy_my_0.083deg_static",
    "--dataset-part bathy",
    "--dataset-version 202311",
    "--variable deptho",
    "--variable deptho_lev",
    "--variable mask",
    "--minimum-longitude",bbox[1],
    "--maximum-longitude",bbox[2],
    "--minimum-latitude",bbox[3],
    "--maximum-latitude",bbox[4],
    "--coordinates-selection-method strict-inside",
    "--netcdf-compression-level 1",
    "--disable-progress-bar",
    "--log-level ERROR",
    "-o",shQuote(dir.out)
  )
  system(command)
  
  ##15 min grid for bgc----
  command <- paste(
    shQuote(path.copernicusmarine),
    "subset",
    "--dataset-id cmems_mod_glo_bgc_my_0.25deg_static",
    "--dataset-part mask",
    "--dataset-version 202406",
    "--variable deptho",
    "--variable deptho_lev",
    "--variable mask",
    "--minimum-longitude",bbox[1],
    "--maximum-longitude",bbox[2],
    "--minimum-latitude",bbox[3],
    "--maximum-latitude",bbox[4],
    "--coordinates-selection-method strict-inside",
    "--netcdf-compression-level 1",
    "--disable-progress-bar",
    "--log-level ERROR",
    "-o",shQuote(dir.out)
  )
  system(command)
}

#Process Netcdf files-----------------------------------------------------------

fill_na_iter <- function(r, w = 3, max_iter = 50) {
  #filled = ras.v[[t]]
  filled <- r
  kernel <- matrix(1, w, w)
  for (i in 1:max_iter) {
    focal_r <- focal(filled, w = kernel, fun = mean, na.rm = TRUE, fillvalue = NA)
    filled[is.na(filled)] <- focal_r[is.na(filled)]
    #if (global(is.na(filled), "sum") == 0) break
    if(length(which(!is.na(depth[]) & is.na(filled[])))==0) break
  }
  return(filled)
}


fn.write_ascii_glorys <- function(dir.in=dir.glorys, dir.stdriver=dir.stdriver, depth=depth){
  # dir.in=dir.glorys
  # dir.stdriver=file.path(dirname(getwd()),'ST drivers','15min')
  # depth = depth.15min
  assign("depth", depth, envir = .GlobalEnv)
  #list netcdf files
  nc.files <- list.files(dir.in,pattern=".nc$",full.names=T)
  nc.files <- nc.files[which(!grepl("static",basename(nc.files)))]
  
  for(i in 1:length(nc.files)){  
    #i=2
    nc = nc_open(nc.files[i])
    nc.vars = names(nc$var)
    print(nc.vars)
    nc.times = as.POSIXct('1950-01-01 00:00')+as.difftime(nc$dim$time$vals,units='hours')
    
    for(v in 1:length(nc.vars)){
      #v=1
      nc.v = ncvar_get(nc,nc.vars[v])
      ndims = length(dim(nc.v))
      
      if(ndims==3){  #bottom temp is the only variable with just 3 dimensions
        ntimes = dim(nc.v)[3]
        dimnames(nc.v)[[3]] <- gsub("-","_",substr(nc.times,1,10))
        out.v = stack()
        nc.tmp = aperm(nc.v,c(2,1,3))
        nc.tmp = nc.tmp[dim(nc.tmp)[1]:1,,]
        brick1 = brick(nc.tmp, crs=crs(depth))
        extent(brick1) = extent(depth)
        ras.v <- resample(rast(brick1),rast(depth))
        
        for(t in 1:nlyr(ras.v)){
          #t=1
          print(paste0(nc.vars[v],"--month ",t," of ",ntimes));flush.console()
          ras.t.filled <- fill_na_iter(ras.v[[t]])
          ras.t <- mask(ras.t.filled,rast(depth))
          out.v <- addLayer(out.v, raster(ras.t))
          rm(ras.t, ras.t.filled); gc()
        }
        names(out.v) <- paste0("X",format(nc.times,"%Y%m"))
        
        #write ascii
        dir.out = file.path(dir.stdriver,paste0(nc.vars[v],'_glorys'))
        if(!dir.exists(dir.out)) dir.create(dir.out)
        writeRaster(out.v,file.path(dir.out,paste0(nc.vars[v])),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
        rm(out.v);gc()
      }
      
      if(ndims==4){
        ntimes = dim(nc.v)[4]
        dimnames(nc.v)[[4]] <- gsub("-","_",substr(nc.times,1,10))
        out.surf = out.mean = out.bott = out.sum = stack() 
        
        #the next 2 lines transpose the array
        nc.tmp = aperm(nc.v,c(2,1,3,4))
        nc.tmp = nc.tmp[dim(nc.tmp)[1]:1,,,]
        
        #loop over timesteps----
        for(t in 1:ntimes){  #big monthly loop takes time, need to parallel
          #t=1
          print(paste0(nc.vars[v],"--month ",t," of ",ntimes));flush.console()
          brick1 = brick(nc.tmp[,,,t],crs=crs(depth))
          extent(brick1) = extent(depth)
          
          #surface layer----
          if(!nc.vars[v] %in% c('chl','no3','nppv','phyc')){
            ras.s = rast(brick1[[1]])
            ras.s = resample(ras.s,rast(depth)) #resample to output resolution
            ras.s.filled <- fill_na_iter(ras.s) #fill missing cells iteratively
            ras.s <- mask(ras.s.filled, rast(depth))  #mask the land
            out.surf = addLayer(out.surf,raster(ras.s)) #add to stack
            rm(ras.s, ras.s.filled);gc()
          }
          
          #bottom layer----
          if(!nc.vars[v] %in% c('chl','no3','nppv','phyc')){
            #fill down to last level, so easier to extract bottom values
            brick2 = approxNA(brick1, method='constant', rule=2)
            ras.b = rast(brick2[[nlayers(brick2)]])
            ras.b = resample(ras.b,rast(depth))
            ras.b.filled <- fill_na_iter(ras.b)
            ras.b <- mask(ras.b.filled, rast(depth))
            out.bott = addLayer(out.bott,raster(ras.b))
            rm(brick2, ras.b, ras.b.filled); gc()
          }
          
          #mean over water column----
          if(!nc.vars[v] %in% c('chl','no3','nppv','phyc')){
            ras.m = mean(rast(brick1),na.rm=T)
            ras.m = resample(ras.m,rast(depth))
            ras.m.filled <- fill_na_iter(ras.m)
            ras.m <- mask(ras.m.filled, rast(depth))
            out.mean = addLayer(out.mean,raster(ras.m))
            rm(ras.m, ras.m.filled);gc()
          }
          
          #sum over water column----
          if(nc.vars[v] %in% c('chl','no3','nppv','phyc')){
            ras.sum = sum(rast(brick1),na.rm=T)
            ras.sum = resample(ras.sum, rast(depth))
            ras.sum.filled <- fill_na_iter(ras.sum)
            ras.sum <- mask(ras.sum.filled, rast(depth))
            out.sum <- addLayer(out.sum, raster(ras.sum))
            rm(ras.sum, ras.sum.filled); gc()
          }
        }
        rm(brick1);gc()
        
        if(nlayers(out.surf)>0) names(out.surf) <- paste0("X",format(nc.times,"%Y%m"))
        if(nlayers(out.bott)>0) names(out.bott) <- paste0("X",format(nc.times,"%Y%m"))
        if(nlayers(out.mean)>0) names(out.mean) <- paste0("X",format(nc.times,"%Y%m"))
        if(nlayers(out.sum)>0) names(out.sum) <- paste0("X",format(nc.times,"%Y%m"))
        
        #write ascii------------------------------------------------------------
        if(nlayers(out.surf)>0){
          dir.out = file.path(dir.stdriver,paste0(nc.vars[v],'_surf_glorys'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(out.surf,file.path(dir.out,paste0(nc.vars[v],"_surf")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
          rm(out.surf);gc()
        }
        
        if(nlayers(out.bott)>0){
          dir.out = file.path(dir.stdriver,paste0(nc.vars[v],'_bott_glorys'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(out.bott,file.path(dir.out,paste0(nc.vars[v],"_bott")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
          rm(out.bott);gc()
        }
        
        if(nlayers(out.mean)>0){
          dir.out = file.path(dir.stdriver,paste0(nc.vars[v],'_mean_glorys'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(out.mean,file.path(dir.out,paste0(nc.vars[v],"_mean")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
          rm(out.mean);gc()
        }
        
        if(nlayers(out.sum)>0){
          dir.out = file.path(dir.stdriver,paste0(nc.vars[v],'_sum_glorys'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(out.sum,file.path(dir.out,paste0(nc.vars[v],"_sum")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
          rm(out.sum); gc()
        }
      }
    }
  }
  #make static maps----
  print("Making static maps");flush.console()
  dirs.stvars = list.dirs(path=dir.stdriver,recursive=F)
  dir.static = file.path(dir.stdriver,"static")
  if(!dir.exists(dir.static)) dir.create(dir.static)
  
  for(i in 1:length(dirs.stvars)){
    #i=1
    stname = paste0(basename(dirname(dirs.stvars[i])),"_",basename(dirs.stvars[i]))
    print(paste0("now making ",i," of ",length(dirs.stvars),": ",stname)); flush.console()
    
    #read output--------------------------------------------------------------------
    #list ecospace output ascii files
    files.st = list.files(dirs.stvars[i],full.names = T,recursive=T)#[-c(1:6)]
    
    #read them into raster stack
    st.stack = stack(files.st)
    #names(st.stack) = paste0(rep(month.abb,24),rep(1997:2020,each=12))
    yrmo = substr(basename(files.st),nchar(basename(files.st))-9,nchar(basename(files.st))-4)
    names(st.stack) = paste0(month.abb[as.numeric(substr(yrmo,5,6))],substr(yrmo,1,4))
    
    #create average map
    st.mean_allyrs = calc(st.stack,fun=mean,na.rm=T)
    st.mean_yr1 = calc(st.stack[[1:12]],fun=mean,na.rm=T)
    st.mean_allyrs = mean(rast(st.stack),na.rm=T)
    st.mean_yr1 = mean(rast(st.stack[[1:12]]), na.rm=T)
    
    plot(st.mean_allyrs, colNA='black')
    plot(st.mean_yr1, colNA='black')
    
    #save static map
    writeRaster(raster(st.mean_allyrs),file.path(dir.static,paste0(gsub("_GLORYS","",stname),"_avg_allyrs")),format='ascii',overwrite=T)
    writeRaster(raster(st.mean_yr1),file.path(dir.static,paste0(gsub("_GLORYS","",stname),"_avg_yr1")),format='ascii',overwrite=T)
    rm(st.stack); gc()
  }
}
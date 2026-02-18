#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' @title Pull GLORYS
#' @description Pull monthly 1/4 degree physical and biogeochemical hindcast products from the 
#' Copernicus marine data store.  Requires the Copernicus marine toolbox CLI application.  
#' Check website https://data.marine.copernicus.eu/products for dataset id and variable names.
#' @param dir.out Output directory to save netcdf files.
#' @param bbox Numeric vector containing spatial domain (W,E,S,N).
#' @param startdate First day to pull data, YYYY-MM-DD. 
#' @param enddate Last day to pull data, YYYY-MM-DD. 
#' @param datid.phy Dataset id for physical data, visit website.
#' @param datid.bgc Dataset id for biogeochemical data, visit website.
#' @param vars.phy Character vector containing physical variable names to download.  Check website 
#' for correct names.
#' @param vars.bgc Character vector containing biogeochemical variable names to download.  Check 
#' website for correct names.
#' @param path.copernicusmarine File path location of copernicusmarine.exe file.
#' @param getbathygrid Logical.  Should bathymetric grids be downloaded?
#' @return One or more netcdf files containing the data.  
#' @examples
#' # example code:
#' \dontrun{fn.pull_glorys(dir.out=dir.glorys, bbox=bbox, startdate = startdate, enddate = enddate, vars.phy = c('thetao','so','bottomT'), vars.bgc = c('chl','nppv','o2','phyc'), path.copernicusmarine <- "C:/Users/dchagaris/OneDrive - University of Florida/copernicusmarine/copernicusmarine.exe")}
#' @export
fn.pull_glorys <- function(dir.out=dir.glorys, 
                           bbox=bbox,
                           startdate = startdate,
                           enddate = enddate,
                           datid.phy = "cmems_mod_glo_phy-all_my_0.25deg_P1M-m",
                           datid.bgc = "cmems_mod_glo_bgc_my_0.25deg_P1M-m",
                           vars.phy = c('thetao','so','bottomT'),
                           vars.bgc = c('chl','nppv','o2','phyc'),
                           path.copernicusmarine = "C:/~/copernicusmarine.exe",
                           getbathygrid=FALSE){
  
  varstring.phy <- paste(paste0("--variable ", vars.phy), collapse = " ")
  varstring.bgc <- paste(paste0("--variable ", vars.bgc), collapse = " ")
  
  #PHYSICAL-------------------------------------------------------
  ##get times----
  if(!is.null(datid.phy)){
  command <- paste(shQuote(path.copernicusmarine),
                   "describe --dataset-id",datid.phy,"--contains time --return-fields coordinates,maximum_value,minimum_value")
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
      "--dataset-id cmems_mod_glo_phy-all_my_0.25deg_P1M-m",
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
  }
  #BIOGEOCHEMICAL----------------------------------------------------------------------
  ##get times----
  if(!is.null(datid.bgc)){
  command <- paste(shQuote(path.copernicusmarine),
                     "describe --dataset-id",datid.bgc,"--contains time --return-fields coordinates,maximum_value,minimum_value")
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
  }
  
  #BATHYMETRY--------------------------------------------------------------
  if(getbathygrid){
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
  
  #15 min grid for bgc----
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
}
#Process Netcdf files-----------------------------------------------------------
#' @title Fill missing cells
#' @description Iteratively applies 3x3 nearest neighbor means to NA cells until all water cells 
#' have data.  Typically used to fill areas near the coastlilne with grids don't align perfectly.
#' @param r A raster layer.
#' @param w Neighborhood area, must be odd.
#' @param max_iter Maximum number to times to iterate until cells are filled. 
#' @param mask Raster layer with land cells NA. 
#' @return Raster layer of dim(r), with all land cells NA and all water cells containing data.  
#' @examples
#' \dontrun{
#' fill_coastal_cells(r, mask=depth)
#' }
#' @export
fill_coastal_cells <- function(r, w = 3, max_iter = 50, mask=depth) {
  #filled = ras.v[[t]]
  filled <- r
  kernel <- matrix(1, w, w)
  for (i in 1:max_iter) {
    focal_r <- focal(filled, w = kernel, fun = mean, na.rm = TRUE, fillvalue = NA)
    filled[is.na(filled)] <- focal_r[is.na(filled)]
    #if (global(is.na(filled), "sum") == 0) break
    if(length(which(!is.na(mask[]) & is.na(filled[])))==0) break
  }
  return(filled)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' @title Create ascii files from GLORYS netcdf
#' @description This function opens a netcdf file, loops through variables and timesteps, and 
#' creates ascii files for input to Ecospace. Includes vertical integration of chl-a and coastline filling.
#' @param nc.files A character variable containing the full file paths of the netcdf files to 
#' process.
#' @param dir.stdriver Parent directory for spatial temporal drivers.  A folder will be created for 
#' each variable in this directory to save the ascii files.
#' @param depth Raster grid of depth, matching the basemap of Ecospacem with land cells NA.  Saved ascii files will have
#' the same dimensions. 
#' @param make.monthly Create the monthly ST driver ascii file.  Default TRUE.
#' @param make.static Also create the static maps for initialization. Default TRUE.
#' @return A series of ascii files saved in dir.stdriver directory.  
#' @examples
#' # example code:
#' \dontrun{
#' files.nc = list.files(dir.glorys, pattern=".nc$", full.names=T)
#' fn.netcdf2ascii_glorys(nc.files = files.nc, depth=depth.15min,
#'                       dir.stdriver=file.path(dirname(getwd()),'ST drivers','15min'))}
#' @export
fn.netcdf2ascii_glorys <- function(nc.files=list.files(dir.glorys,pattern=".nc$",full.names=T), 
                                   dir.stdriver=dir.stdriver, depth=depth, make.monthly=TRUE, make.static=TRUE){
  # dir.in=dir.glorys
  # dir.stdriver=file.path(dirname(getwd()),'ST drivers','15min')
  # depth = depth.15min
  # nc.files = files.nc
  assign("depth", depth, envir = .GlobalEnv)
  #list netcdf files
  #nc.files <- list.files(dir.in,pattern=".nc$",full.names=T)
  #nc.files <- nc.files[which(!grepl("static",basename(nc.files)))]
  if(make.monthly){
  for(i in 1:length(nc.files)){  
    #i=1
    nc = nc_open(nc.files[i])
    nc.vars = names(nc$var)
    print(nc.vars)
    nc$dim$time
    if(i==1) nc.times = as.POSIXct('1950-01-01 00:00')+as.difftime(nc$dim$time$vals,units='hours')
    if(i==2) nc.times = as.POSIXct('1950-01-01 00:00')+as.difftime(nc$dim$time$vals,units='secs')
    zthickness <- nc$dim$depth$vals
    
    st.existing <- list.dirs(dir.stdriver,full.names=F,recursive=F)
    st.existing <- st.existing[which(!st.existing %in% c('hoard','plots','static'))]
    
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
          ras.t.filled <- fill_coastal_cells(ras.v[[t]])
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
        dim(nc.v)
        nc.tmp = aperm(nc.v,c(2,1,3,4))
        nc.tmp = nc.tmp[dim(nc.tmp)[1]:1,,,]
        
        #loop over timesteps----
        for(t in 1:ntimes){  #big monthly loop takes time, need to parallel
          #t=1
          print(paste0(nc.vars[v],"--month ",t," of ",ntimes));flush.console()
          brick1 = brick(nc.tmp[,,,t],crs=crs(depth))
          extent(brick1) = extent(depth)
          
          #surface layer----
          #if(!nc.vars[v] %in% c('chl','no3','nppv','phyc')){
          if(!nc.vars[v] %in% c('x')){
            ras.s = rast(brick1[[1]])
            ras.s = resample(ras.s,rast(depth)) #resample to output resolution
            ras.s.filled <- fill_coastal_cells(ras.s) #fill missing cells iteratively
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
            ras.b.filled <- fill_coastal_cells(ras.b)
            ras.b <- mask(ras.b.filled, rast(depth))
            out.bott = addLayer(out.bott,raster(ras.b))
            rm(brick2, ras.b, ras.b.filled); gc()
          }
          
          #mean over water column----
          
          # Should show a terra method for SpatRaster

          if(!nc.vars[v] %in% c('chl','no3','nppv','phyc')){
            ras.m = mean(rast(brick1),na.rm=T)
            #ras.m = terra::weighted.mean(x=rast(brick1),w=zthickness, na.rm=T)
            ras.m = resample(ras.m,rast(depth))
            ras.m.filled <- fill_coastal_cells(ras.m)
            ras.m <- mask(ras.m.filled, rast(depth))
            out.mean = addLayer(out.mean,raster(ras.m))
            rm(ras.m, ras.m.filled);gc()
          }
          
          #sum over water column----
          if(nc.vars[v] %in% c('chl','no3','nppv','phyc')){
            #ras.sum = sum(rast(brick1),na.rm=T)
            ras.sum <- app(rast(brick1), fun = 
                             function(v) {
                               if (all(is.na(v))) {
                                 NA               # keep land as NA
                               } else {
                                 sum(v * zthickness, na.rm = TRUE)
                               }
                             })
            ras.sum = resample(ras.sum, rast(depth))
            ras.sum.filled <- fill_coastal_cells(ras.sum, mask=depth)
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
          dir.out = file.path(dir.stdriver,gsub('_glor','_surf_glorys',nc.vars[v]))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(out.surf,file.path(dir.out,gsub('_glor','_surf',nc.vars[v])),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
          rm(out.surf);gc()
        }
        
        if(nlayers(out.bott)>0){
          dir.out = file.path(dir.stdriver,gsub('_glor','_bott_glorys',nc.vars[v]))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(out.bott,file.path(dir.out,gsub('_glor','_bott',nc.vars[v])),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
          rm(out.bott);gc()
        }
        
        if(nlayers(out.mean)>0){
          dir.out = file.path(dir.stdriver,gsub('_glor','_mean_glorys',nc.vars[v]))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(out.mean,file.path(dir.out,gsub('_glor','_mean',nc.vars[v])),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
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
  }
  #make static maps----
  if(make.static){
  print("Making static maps");flush.console()
  dirs.stvars = list.dirs(path=dir.stdriver,recursive=F)
  dir.static = file.path(dir.stdriver,"static")
  if(!dir.exists(dir.static)) dir.create(dir.static)
  dirs.stvars = dirs.stvars[-c(which(grepl("static",dirs.stvars)),which(grepl("plots",dirs.stvars)))]
  
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
    
    #plot(st.mean_allyrs, colNA='black')
    #plot(st.mean_yr1, colNA='black')
    
    #save static map
    writeRaster(raster(st.mean_allyrs),file.path(dir.static,paste0(gsub("_GLORYS","",stname),"_avg_allyrs")),format='ascii',overwrite=T)
    writeRaster(raster(st.mean_yr1),file.path(dir.static,paste0(gsub("_GLORYS","",stname),"_avg_yr1")),format='ascii',overwrite=T)
    rm(st.stack); gc()
  }
  }
}#eof


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' @title List variables and get URLs for CEFI data.
#' @description This function will get information on CEFI data and return a dataframe that can be 
#' subset for pulling data.
#' @param dir.cefi Parent directory for CEFI data.
#' @param file.cefivars Name of csv file to be created (or read) containing cefi variables. Default 
#' is CEFI_vars.csv.
#' @param region CEFI region nwa=northwest Atlantic, nep=northeast Pacific
#' @examples
#' # example code:
#' \dontrun{vars.cefi <- fn.get_cefi_vars(dir.cefi=dir.cefi, file.cefivars='CEFI_vars.csv', region='nwa')}
#' @export
fn.get_cefi_vars = function(dir.cefi=dir.cefi, file.cefivars='CEFI_vars.csv', region='nwa'){ 
  if(!file.exists(file.path(dir.cefi,file.cefivars))){
    #print("getting CEFI variable list");flush.consolse()
    if(region=='nwa'){
      vars.hcast =  read_html("https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.hindcast.html")
      vars.refore = read_html("https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.seasonal_reforecast.html")
      vars.fore =   read_html("https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.seasonal_forecast.html")
      vars.decade = read_html("https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.decadal_forecast.html")
    }
    vars.hcast = html_table(vars.hcast, fill=T)[[1]]
    vars.refore = html_table(vars.refore, fill=T)[[1]]
    vars.fore = html_table(vars.fore, fill=T)[[1]]
    vars.decade = html_table(vars.decade, fill=T)[[1]]
    
    vars.cefi <- rbind(vars.hcast, vars.refore, vars.fore, vars.decade)
    write.csv(vars.cefi,file.path(dir.cefi,file.cefivars), row.names=F)
  } else{
    vars.cefi <- read.csv(file.path(dir.cefi,file.cefivars))
  }
  return(vars.cefi)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' @title Pull CEFI data
#' @description Downloads the full netcdf for MOM6-COBALT variables from 
#' https://downloads.psl.noaa.gov/Projects/CEFI/regional_mom6 using the httr::GET function.  It 
#' applies the spatial subset based on the extent of depth and writes a new netcdf file.
#' @param dir.cefi Parent directory for CEFI data, where the netcdf files will be saved.
#' @param vars.cefi A dataframe returned by fn.get_cefi_vars(), likely subset to fewer variables. 
#' Must contain the variables "cefi_filename" and "cefi_rel_path".
#' @param reg CEFI region, either northwest_atlantic or northeast_pacific
#' @param depth Raster grid of depth, matching the basemap of Ecospace.  Saved ascii files will have
#' the same dimensions. 
#' @param delete.full.nc Logical. TRUE will delete the full nc files if they exist.  Default is False.
#' @examples
#' # example code:
#' \dontrun{fn.pull_cefi(vars.cefi=vars.cefi, dir.cefi=dir.cefi, reg='northwest_atlantic', depth=depth.5min, delete.full.nc=TRUE)}
#' @export
fn.pull_cefi = function(vars.cefi=vars.cefi, dir.cefi=dir.cefi, reg='northwest_atlantic', depth=depth.5min, delete.full.nc=FALSE){
  path.downloads <- "https://downloads.psl.noaa.gov/Projects/CEFI/regional_mom6"
  for(i in 1:nrow(vars.cefi)){
    #i=1
    path.thredds <- file.path(path.downloads,vars.cefi$cefi_rel_path[i])
    file.nc = vars.cefi$cefi_filename[i]
    file.out = file.path(dir.cefi,file.nc)
    file.nc.sub = gsub(".full.",".subset.",gsub(".nc",paste0("-",paste(abs(t(as.matrix(extent(depth)))),collapse="-"),".nc"), file.out))
    
    #download full netcdf if needed----
    if(file.exists(file.out)){
      message("file ",i," of ",nrow(vars.cefi)," already exists\n---- skipping download ",file.nc) 
    } else{
      message("downloading file ",i," of ",nrow(vars.cefi),"\n---- ",file.nc,"\n---- from ",path.thredds)
      url = file.path(path.thredds, file.nc )
      httr::GET(url, httr::write_disk(file.out, overwrite=T), httr::progress())
    }
    
    #spatial subset----
    if(file.exists(file.nc.sub)){
      message("file ",i," of ",nrow(vars.cefi)," subset already exists\n---- skipping spatial subset of ",file.nc) 
    } else{
      nc <- nc_open(file.out)
      nc.vars = names(nc$var)
      message("- Subsetting variable ",i," of ",nrow(vars.cefi),": ",nc.vars)
      
      vinfo <- nc$var[[1]]
      ndims = length(vinfo$dim)
      var_units <- vinfo$units
      
      #get timestamps----
      # refdate = as.Date(substr(nc$dim$time$units,nchar(nc$dim$time$units)-9,nchar(nc$dim$time$units)))
      # nc.times = seq(refdate, by='month', length.out=length(nc$dim$time$vals))
      time_units <- nc$dim$time$units       # original time units
      time_vals  <- nc$dim$time$vals        # original numeric time values
      
      #Get coordinate variables----
      lon <- ncvar_get(nc, "lon")
      lat <- ncvar_get(nc, "lat")
      zvals <- nc$dim$z_l$vals
      
      # Define lat/lon----
      lat_bounds <- extent(depth)[3:4]
      lon_bounds <- extent(depth)[1:2]
      
      # Find lat/lon indices for subsetting
      lon_idx <- which(lon >= lon_bounds[1] & lon <= lon_bounds[2])
      lat_idx <- which(lat >= lat_bounds[1] & lat <= lat_bounds[2])
      
      #sub arrays for centers
      lon_sub <- lon[lon_idx]
      lat_sub <- lat[lat_idx]
      
      # # Compute center spacings; if descending, use abs
      # dx <- abs(mean(diff(lon_sub)))
      # dy <- abs(mean(diff(lat_sub)))
      # 
      # # Extent from centers (pad by half a cell)
      # ext_nc <- extent(min(lon_sub) - dx/2, max(lon_sub) + dx/2,
      #                  min(lat_sub) - dy/2, max(lat_sub) + dy/2)
      # 
      # nrow_nc <- length(lat_idx)  # rows -> lat
      # ncol_nc <- length(lon_idx)  # cols -> lon
      # ntime   <- length(nc.times)
  
      # Define dimensions
      lon_dim  <- ncdim_def("lon", "degrees_east", lon_sub)
      lat_dim  <- ncdim_def("lat", "degrees_north", lat_sub)
      time_dim <- ncdim_def("time", time_units, time_vals)
      
      #read and subset variable----
      if(ndims==3){
        start <- c(min(lon_idx), min(lat_idx), 1)
        count <- c(length(lon_idx), length(lat_idx), -1)#length(zvals), 1)#length(time_idx))  # one time slice
        nc.sub <- ncvar_get(nc, var=nc.vars, start=start, count=count)
        
        #save subset as new netcdf----
        # Define variable
        var_def <- ncvar_def(nc.vars, var_units, dim = list(lon_dim, lat_dim, time_dim), compression = 4)
      
        # Create file and write
        nc_new <- nc_create(file.nc.sub, list(var_def))
        ncvar_put(nc_new, nc.vars, nc.sub)
      }
    
      nc_close(nc_new)
      nc_close(nc)
    
  }
      if(delete.full.nc) file.remove(file.out)
  }
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' @title Download full CEFI netcdf files
#' @description Downloads the full netcdf for MOM6-COBALT variables from 
#' https://downloads.psl.noaa.gov/Projects/CEFI/regional_mom6.  Uses the httr::GET function.
#' @param dir.cefi Parent directory for CEFI data, where the netcdf files will be saved.
#' @param vars.cefi A dataframe returned by fn.get_cefi_vars(), likely subset to fewer variables. 
#' Must contain the variables "cefi_filename" and "cefi_rel_path".
#' @param reg CEFI region, either northwest_atlantic or northeast_pacific
#' @examples
#' # example code:
#' \dontrun{vars.cefi <- fn.get_cefi_vars(dir.cefi=dir.cefi, file.cefivars='CEFI_vars.csv', region='nwa')}
#' @export
fn.download_cefi_full_nc = function(vars.cefi=vars.cefi, dir.cefi=dir.cefi, reg='northwest_atlantic'){
  path.downloads <- "https://downloads.psl.noaa.gov/Projects/CEFI/regional_mom6"
  for(i in 1:nrow(vars.cefi)){
    #i=1
    path.thredds <- file.path(path.downloads,vars.cefi$cefi_rel_path[i])
    file.nc = vars.cefi$cefi_filename[i]
    file.out = file.path(dir.cefi,file.nc)
    url = file.path(path.thredds, file.nc )
    if(!file.exists(file.out)){
      cat("downloading file",i,"of",nrow(vars.cefi),"\n----",file.nc,"\n---- from",path.thredds)
      httr::GET(url, httr::write_disk(file.out, overwrite=T), httr::progress())
    } else{
      cat("file",i,"of",nrow(vars.cefi),"already exists\n---- skipping",file.nc) 
    }
  }
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' @title Create ascii files from CEFI netcdf
#' @description This function processes netcdf files from local machine, 
#' downloaded with fn.download_cefi_full_nc().  It peforms spatial subset based on depth map, 
#' regrids and resamples to depth, fills missing water cells, and writes ascii files.
#' @param nc.files A character variable containing the file paths of the netcdf files to process.
#' @param dir.stdriver Parent directory for spatial temporal drivers.  A folder will be created for
#' each variable in this directory to save the ascii files.
#' @param depth Raster grid of depth, matching the basemap of Ecospace.  Saved ascii files will have
#' the same dimensions. 
#' @param make.monthly Logical.  TRUE will run all the code to create monthly raster stacks.  Use
#' FALSE if only want to make static maps.
#' @param make.static Logical.  TRUE will create a static map (average) for each variable to use as
#' Ecospace basemap layer.
#' @return A series of ascii files saved in dir.stdriver directory.  
#' @examples
#' \dontrun{
#' # example code:
#' files.nc = list.files(dir.cefi, pattern=".nc$", full.names=T)
#' fn.netcdf2ascii_cefi(nc.files = files.nc, depth=depth.15min,
#'                       dir.stdriver=file.path(dirname(getwd()),'ST drivers','5min'))}
#' @export

fn.netcdf2ascii_cefi <- function(nc.files=list.files(path=dir.cefi, pattern=".nc$", full.names=T), depth=depth.5min,
                                 dir.stdriver=dir.stdriver, make.monthly=TRUE, make.static=FALSE){
  # nc.files = files.nc
  # depth=depth.15min
  # dir.stdriver=file.path(dirname(getwd()),'ST drivers','15min')
  # make.monthly=TRUE
  # make.static=TRUE
  
  if(make.monthly){
  for(i in 1:length(nc.files)){
    #i=1
    nc <- nc_open(nc.files[i])
    on.exit(nc_close(nc), add = TRUE)
    nc.vars = names(nc$var)
    print(nc.vars)
    ndims = nc$ndims
    message("- Processing variable ",i," of ",length(nc.files),": ",nc.vars)
    nc.sub <- ncvar_get(nc, var=nc.vars)
    dim(nc.sub)
    message("-- Variable read into memory.")
    
    #get timestamps----
    refdate = as.Date(substr(nc$dim$time$units,nchar(nc$dim$time$units)-9,nchar(nc$dim$time$units)))
    nc.times = seq(refdate, by='month', length.out=length(nc$dim$time$vals))
    ntime   <- length(nc.times)
    
    #Get coordinate variables----
    lon <- ncvar_get(nc, "lon")
    lat <- ncvar_get(nc, "lat")
    zvals <- nc$dim$z_l$vals
    
    # Define lat/lon----
    lat_bounds <- extent(depth)[3:4]
    lon_bounds <- extent(depth)[1:2]
    
    # Find lat/lon indices for subsetting
    lon_idx <- which(lon >= lon_bounds[1] & lon <= lon_bounds[2])
    lat_idx <- which(lat >= lat_bounds[1] & lat <= lat_bounds[2])
    
    #sub arrays for centers
    lon_sub <- lon[lon_idx]
    lat_sub <- lat[lat_idx]
    
    # Compute center spacings; if descending, use abs
    dx <- abs(mean(diff(lon_sub)))
    dy <- abs(mean(diff(lat_sub)))
    
    # Extent from centers (pad by half a cell)
    ext_nc <- extent(min(lon_sub) - dx/2, max(lon_sub) + dx/2,
                     min(lat_sub) - dy/2, max(lat_sub) + dy/2)
    
    nrow_nc <- length(lat_idx)  # rows -> lat
    ncol_nc <- length(lon_idx)  # cols -> lon

    # ---- build stack on native grid ----
    message("-- Building stack on native grid...")
    var.stack = stack()
    for (t in 1:ntime) {
      #t=1
      
      mat <- nc.sub[lon_idx, lat_idx, t, drop = FALSE][,,1]
      
      # flip north-south
      mat <- mat[, ncol(mat):1]
      
      r_k <- raster(
        nrows = nrow_nc,
        ncols = ncol_nc,
        ext   = ext_nc,
        crs   = crs(depth)
      )
      
      values(r_k) <- as.vector(mat)
      #plot(r_k)
      
      # mat <- nc.sub[lon_idx, lat_idx,t] #subset for lat lon
      # mat <- t(mat)                 # swap lon,lat → lat,lon
      # if (lat[1] < lat[length(lat)]){mat <- mat[nrow(mat):1, ]}   # flip latitude only if ascending
      # 
      # #mat <- t(mat[, dim(mat)[2]:1])  # assumes [lat, lon] matrix, need to flip x-axis (lon) and transpose
      # r_k <- raster(nrows = nrow_nc, ncols = ncol_nc,
      #               ext = ext_nc, crs = crs(depth))
      # # raster fills values column-wise; t() keeps spatial layout
      # values(r_k) <- as.vector((mat))
      
      var.stack <- addLayer(var.stack, r_k)
    }
    
    names(var.stack) <- nc.times

    # ---- Regrid to depth (sq 1/12°) ----
    # Project CRS if needed
    message("-- Regridding to depth...")
    depth_proj <- depth
    if (!compareCRS(depth_proj, var.stack)) {
      depth_proj <- projectRaster(depth_proj, crs = crs(var.stack), method = "bilinear")
    }
    
    # Resample the stack to the depth grid (now square cells)
    var_on_depth <- resample(var.stack, depth_proj, method = "bilinear")
    
    # Crop & mask by depth coverage
    var_crop   <- crop(var_on_depth, extent(depth_proj))
    var_masked <- mask(var_crop, !is.na(depth_proj))
    var.regrid <- readAll(var_masked)
    
    # ---- quick checks ----
    message("-- Native res: ", paste(round(res(var.stack),4), collapse = ", "),
            " | Native grid: ", paste(dim(var.stack)[1:2],collapse="x"),
            " | Depth res: ", paste(round(res(depth_proj),4), collapse = ", "),
            " | Depth grid: ", paste(dim(depth)[1:2],collapse="x"),
            " | Resampled res: ", paste(round(res(var.regrid),4), collapse = ", "),
            " | Resampled grid: ", paste(dim(var.regrid)[1:2],collapse="x"))
    
    
    #---- fill NA water cells----
    message('-- Filling water cells...')
    var.out = stack()
    for(s in 1:nlayers(var.regrid)){
      #s=1
      message("----",s," of ",nlayers(var.regrid))
      ras.s = rast(var.regrid[[s]])
      ras.s = resample(ras.s,rast(depth)) #resample to output resolution
      ras.s.filled <- fill_coastal_cells(ras.s, mask='depth') #fill missing cells iteratively
      ras.s <- mask(ras.s.filled, rast(depth))  #mask the land
      #plot(ras.s, colNA='black',main=paste(nc.vars,var_units))
      var.out = addLayer(var.out,raster(ras.s)) #add to stack
      rm(ras.s, ras.s.filled);gc()
    }  

    #write ascii------------------------------------------------------------
    message("-- Writing ascii files...")
    if(nlayers(var.out)>0){
      dir.out = file.path(dir.stdriver,paste0(nc.vars,'_cefi'))
      if(!dir.exists(dir.out)) dir.create(dir.out)
      writeRaster(var.out,file.path(dir.out,nc.vars),bylayer=T,suffix=format(nc.times,"%Y%m"),format='ascii',overwrite=T)
      rm(var.out);gc()
    }
  }
  }
  
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' @title Make static maps from monthly ST files.
#' @description This stacks all the ascii in stdriver directory and compute long term mean and 
#' year-1 mean layers, and saves them as ascii files in the 'static' folder.
#' @param dir.stdriver Directory for ST driver variable.  A 'static' folder will be created one 
#' level up if it doesn't exist.
#' @return Saves an ascii file to the static folder.  
#' @examples
#' # example code:
#' \dontrun{fn.make_static_maps(dir.stdriver=file.path(dir.stdriver,"tos_cefi"))}
#' @export
fn.make_static_maps <- function(dir.stdriver=dir.stdriver){
  #dir.stdriver = list.dirs(dir.stdriver,full.names=T, recursive=F)[1]
  message("- Making static maps for ",basename(dir.stdriver))
  dir.static = file.path(dirname(dir.stdriver),"static")
  if(!dir.exists(dir.static)) dir.create(dir.static)
  
  stname = paste0(basename(dirname(dir.stdriver)),"_",basename(dir.stdriver))
  #read output--------------------------------------------------------------------
  #list ecospace output ascii files
  files.st = list.files(dir.stdriver,full.names = T,recursive=T)#[-c(1:6)]

  #read them into raster stack
  st.stack = stack(files.st)
  #names(st.stack) = paste0(rep(month.abb,24),rep(1997:2020,each=12))
  yrmo = substr(basename(files.st),nchar(basename(files.st))-9,nchar(basename(files.st))-4)
  names(st.stack) = paste0(month.abb[as.numeric(substr(yrmo,5,6))],substr(yrmo,1,4))
  
  #create average map
  #st.mean_allyrs = calc(st.stack,fun=mean,na.rm=T)
  #st.mean_yr1 = calc(st.stack[[1:12]],fun=mean,na.rm=T)
  st.mean_allyrs = mean(rast(st.stack),na.rm=T)
  st.mean_yr1 = mean(rast(st.stack[[1:12]]), na.rm=T)
  
  #plot(st.mean_allyrs, colNA='black')
  #plot(st.mean_yr1, colNA='black')
  
  #save static map
  writeRaster(raster(st.mean_allyrs),file.path(dir.static,paste0(gsub("_GLORYS","",stname),"_avg_allyrs")),format='ascii',overwrite=T)
  writeRaster(raster(st.mean_yr1),file.path(dir.static,paste0(gsub("_GLORYS","",stname),"_avg_yr1")),format='ascii',overwrite=T)
  rm(st.stack); gc()
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' @title Make monthly climatology maps.
#' @description This stacks all the ascii in stdriver directory and compute long term mean for each
#' month.
#' @param dir.stdriver Directory for ST driver variable.
#' @return Saves 12 ascii files to the variable folder with suffix '9999mm', where mm is numeric
#' month.  
#' @examples
#' \dontrun{
#' # example code:
#' for(i in 1:length(dirs.stvars)){
#' fn.make_monthly_climatology_maps(dir.stdriver=dirs.stvars[i])
#' }}
#' @export
fn.make_monthly_climatology_maps <- function(dir.stdriver=dir.stdriver){
  #dir.stdriver = list.dirs(dir.stdriver,full.names=T, recursive=F)[1]
  message("- Making monthly climatology maps for ",basename(dir.stdriver))

  stname = paste0(basename(dirname(dir.stdriver)),"_",basename(dir.stdriver))
  #read output--------------------------------------------------------------------
  #list ecospace output ascii files
  files.st = list.files(dir.stdriver,full.names = T,recursive=T)#[-c(1:6)]
  vars<-unique(basename(dirname(files.st)))
  
  for (ivar in vars) {
      
    #ivar<-vars[1]
    ifiles.st<-files.st[grepl(ivar,files.st)]
    month.id = substr(basename(ifiles.st),nchar(basename(ifiles.st))-5,nchar(basename(ifiles.st))-4)
    #read them into raster stack
    st.stack = raster::stack(ifiles.st)
    #names(st.stack) = paste0(rep(month.abb,24),rep(1997:2020,each=12))
    yrmo = substr(basename(ifiles.st),nchar(basename(ifiles.st))-9,nchar(basename(ifiles.st))-4)
    names(st.stack) = paste0(month.abb[as.numeric(substr(yrmo,5,6))],substr(yrmo,1,4))
    
    monthly_means <- stackApply(st.stack, indices = month.id, fun = mean, na.rm = TRUE)
    
    writeRaster(monthly_means,file.path(dir.stdriver,ivar),bylayer=T,suffix=paste0("9999",sort(unique(month.id))),format='ascii',overwrite=T)
  
    rm(st.stack); gc()
  }
}





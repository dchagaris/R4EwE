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
  startdate.d <- startdate[which(names(startdate)==datid.phy)]
  enddate.d <- enddate[which(names(enddate)==datid.phy)]
  if(as.Date(startdate.d)<time.range[2]){
    startdate1 = max(as.Date(startdate.d),time.range[1])
    enddate1 = min(as.Date(enddate.d),time.range[2])
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
  startdate.d <- startdate[which(names(startdate)==datid.bgc)]
  enddate.d <- enddate[which(names(enddate)==datid.bgc)]
  if(as.Date(startdate.d)<time.range[2]){
    startdate1 = max(as.Date(startdate.d),time.range[1])
    enddate1 = min(as.Date(enddate.d),time.range[2])
    
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



#' Append two CMEMS GLORYS monthly NetCDF files along time
#'
#' @description
#' Appends a newer CMEMS GLORYS monthly NetCDF file to an older one along the
#' \code{time} dimension and writes out a single merged NetCDF file.
#'
#' The two input files are assumed to be identical in structure (variables,
#' dimensions, metadata) and to differ only in their time coverage.
#'
#' All data variables are appended automatically; coordinate variables
#' (longitude, latitude, depth, time) are preserved from the original file.
#'
#' @param old_file Character string giving the file path to the older NetCDF file.
#' @param new_file Character string giving the file path to the newer NetCDF file.
#'
#' @details
#' The function performs basic compatibility checks on variable and dimension names.
#' Data are appended along the time dimension, which is assumed to be the last
#' dimension for all variables (as is standard for CMEMS GLORYS monthly products).
#'
#' Variables are processed one at a time to limit memory usage.
#'
#' @return
#' Invisibly returns the path to the output file (\code{out_file}).
#'
#' @seealso
#' \code{\link[ncdf4]{nc_open}}, \code{\link[ncdf4]{nc_create}},
#' \code{\link[abind]{abind}}
#'
#' @examples
#' \dontrun{
#' append_glorys_monthly_nc(
#'   old_file = "glorys_monthly_1993_2015.nc",
#'   new_file = "glorys_monthly_2016_2024.nc",
#'   out_file = "glorys_monthly_1993_2024.nc"
#' )
#' }
#'
#' @export
fn.append_netcdf_glorys <- function(old_file, new_file) {
  
  #old_file <- file.path(dir.glorys,inventory$file[1])
  #new_file <- file.path(dir.glorys,files.nc[2])
  out_file <- file.path(dirname(old_file),
                        paste0(substr(basename(old_file),1,nchar(basename(old_file))-13),
                               substr(basename(new_file),nchar(basename(new_file))-12,nchar(basename(new_file)))))
  
  stopifnot(file.exists(old_file), file.exists(new_file))
  
  library(ncdf4)
  library(abind)
  
  # Open files----
  nc_old <- nc_open(old_file)
  nc_new <- nc_open(new_file)
  
  on.exit({
    try(nc_close(nc_old), silent = TRUE)
    try(nc_close(nc_new), silent = TRUE)
  })
  
  # Basic compatibility checks----
  if (!identical(names(nc_old$var), names(nc_new$var))) {
    stop("Variable names do not match between NetCDF files.")
  }
  
  if (!identical(names(nc_old$dim), names(nc_new$dim))) {
    stop("Dimension names do not match between NetCDF files.")
  }
  
  # Time----
  time_old <- ncvar_get(nc_old, "time")
  time_new <- ncvar_get(nc_new, "time")
  
  if (min(time_new) <= max(time_old)) {
    warning("New file time overlaps old file time.")
  }
  
  time_all <- c(time_old, time_new)
  
  # Define dimensions (copied from OLD file)----
  dims_old <- nc_old$dim
  
  dim_defs <- vector("list", length(dims_old))
  names(dim_defs) <- names(dims_old)
  
  for (d in names(dims_old)) {
    if (d == "time") {
      dim_defs[[d]] <- ncdim_def(
        "time",
        dims_old$time$units,
        time_all,
        unlim = TRUE
      )
    } else {
      dim_defs[[d]] <- ncdim_def(
        d,
        dims_old[[d]]$units,
        dims_old[[d]]$vals
      )
    }
  }
  
  # Define variables----
  coord_vars <- names(dim_defs)
  var_names  <- setdiff(names(nc_old$var), coord_vars)
  
  var_defs <- vector("list", length(var_names))
  names(var_defs) <- var_names
  
  for (v in var_names) {
    vinfo <- nc_old$var[[v]]
    
    dims <- lapply(vinfo$dimids + 1, function(i) {
      dim_defs[[ names(dims_old)[i] ]]
    })
    
    var_defs[[v]] <- ncvar_def(
      name    = v,
      units   = vinfo$units,
      dim     = dims,
      missval = vinfo$missval,
      prec    = ifelse(vinfo$prec == "double", "double", "float")
    )
  }
  
  # Create output file----
  nc_out <- nc_create(out_file, var_defs)
  
  # Append variables----
  for (v in var_names) {
    message("Appending variable: ", v)
    
    v_old <- ncvar_get(nc_old, v)
    v_new <- ncvar_get(nc_new, v)
    
    # append along TIME (assumed last dimension)
    time_dim <- length(dim(v_old))
    v_all <- abind(v_old, v_new, along = time_dim)
    
    ncvar_put(nc_out, v, v_all)
  }
  
  nc_close(nc_out)
  
  invisible(out_file)
  
  
  rm(nc_old, nc_new, nc_out)
  gc()
  
  
  dir.archive <- file.path(dirname(old_file),'archive')
  dir.create(dir.archive)
  file.copy(from=old_file, to=dir.archive, overwrite=T)
  #file.copy(from=new_file, to=dir.archive, overwrite=T)
  file.remove(old_file)
  file.remove(new_file)
} #eof


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
                                   dir.stdriver=dir.stdriver, depth=depth, scalePP=TRUE,do.vars=NULL,
                                   make.monthly=TRUE, make.static=TRUE, make.climatology=TRUE){
  library('terra')
  library('raster')
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
    #print(nc.vars)
    nc$dim$time
    time_units = nc$dim$time$units
    time_refdate = as.POSIXct(paste0(gsub("[a-zA-Z ]", "", time_units)," 00:00"))
    time_stepsize = ifelse(grepl("hours",time_units),"hours",ifelse(grepl("sec",time_units),'secs',ifelse(grepl("day",time_units),'days',NA)))
    nc.times = time_refdate+as.difftime(nc$dim$time$vals,units=time_stepsize)
    zvals <- nc$dim$depth$vals
    zthickness <- c(diff(zvals-zvals[1]),tail(diff(zvals-zvals[1]),1))
    
    st.existing <- list.dirs(dir.stdriver,full.names=F,recursive=F)
    st.existing <- st.existing[which(!st.existing %in% c('hoard','plots','static'))]
    
    for(v in 1:length(nc.vars)){
      #v=1
      if(!is.null(do.vars) & !nc.vars[v] %in% do.vars) next
      
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
        plot(ras.v)
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
        dir.out = file.path(dir.stdriver,nc.vars[v])
        if(!dir.exists(dir.out)) dir.create(dir.out)
        writeRaster(out.v,file.path(dir.out,gsub("_glorys","",nc.vars[v])),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
        rm(out.v);gc()
      }
      
      if(ndims==4){
        ntimes = dim(nc.v)[4]
        dimnames(nc.v)[[4]] <- gsub("-","_",substr(nc.times,1,10))
        out.surf = out.mean = out.bott = out.sum = chlos.zint = out.surf.norm = out.sum.norm = chlos.zint.norm = stack() 
        
        #the next 2 lines transpose the array
        dim(nc.v)
        nc.tmp = aperm(nc.v,c(2,1,3,4))
        nc.tmp = nc.tmp[dim(nc.tmp)[1]:1,,,]
        
        #loop over timesteps----
        for(t in 1:ntimes){  #big monthly loop takes time, need to parallel
          #t=1
          message(paste0(nc.vars[v],"--month ",t," of ",ntimes));flush.console()
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
          #plot(out.surf)
          if(nc.vars[v] %in% c('chl','no3','nppv','phyc')){
            #ras.sum = sum(rast(brick1),na.rm=T)
 
            # zthickness50 <- zthickness
            # zthickness50[floor(zthickness/10)>5] <- 0
            # ras.sum <- app(rast(brick1), fun = 
            #                  function(v) {
            #                    if (all(is.na(v))) {
            #                      NA               # keep land as NA
            #                    } else {
            #                      sum(v * zthickness50, na.rm = TRUE)
            #                    }
            #                  })
            ras.sum <- app(rast(brick1), fun = function(v) {
               if (all(is.na(v))) {
                 NA               # keep land as NA
               } else {
                 #sum(v * zthickness, na.rm = TRUE) #old way, entire watercolumn
                 # Calculate cumulative depth for each layer
                 # Find layers within first 50m
                 within_50m <- zvals <= 50
                 # For layers that partially extend beyond 50m, 
                 # adjust thickness to only include the portion within 50m
                 adjusted_thickness <- zthickness
                 if(any(zvals > 50 & c(0, zvals[-length(zvals)]) < 50)) {
                   partial_layer <- which(zvals > 50 & c(0, zvals[-length(zvals)]) < 50)
                   adjusted_thickness[partial_layer] <- 50 - c(0, zvals[-length(zvals)])[partial_layer]
                 }
                 cbind(zvals,zthickness,adjusted_thickness)
                 
                 # Sum only the layers within 50m
                 sum(v[within_50m] * adjusted_thickness[within_50m], na.rm = TRUE)
               }
             })
            #plot(ras.sum)

            ras.sum = resample(ras.sum, rast(depth))
            ras.sum.filled <- fill_coastal_cells(ras.sum, mask=depth)
            ras.sum <- mask(ras.sum.filled, rast(depth))
            out.sum <- addLayer(out.sum, raster(ras.sum))
            rm(ras.sum, ras.sum.filled); gc()
          }
        }
        rm(brick1);gc()
        
        #change units
        if(nc.vars[v]=='nppv'){
          message('changing nppv units from mg C m-3 day-1 to g C m-3 day-1')
          out.surf <- out.surf/1000
          out.sum <- out.sum/1000
        }
        
        if(nc.vars[v]=='phyc'){
          message('changing phyc units from mmol C m-3 to g C m-3')
          out.surf <- out.surf*.012
          out.sum <- out.sum*.012
        }
        
        if(nlayers(out.surf)>0) out.surf <- round(out.surf,4)
        if(nlayers(out.bott)>0) out.bott <- round(out.bott,4)
        if(nlayers(out.mean)>0) out.mean <- round(out.mean,4)
        if(nlayers(out.sum)>0) out.sum <- round(out.sum,4)
        
        if(nlayers(out.surf)>0) names(out.surf) <- paste0("X",format(nc.times,"%Y%m"))
        if(nlayers(out.bott)>0) names(out.bott) <- paste0("X",format(nc.times,"%Y%m"))
        if(nlayers(out.mean)>0) names(out.mean) <- paste0("X",format(nc.times,"%Y%m"))
        if(nlayers(out.sum)>0) names(out.sum) <- paste0("X",format(nc.times,"%Y%m"))
        

        #normalize PP to mean of first year----
        if(nc.vars[v] %in% c('chl','nppv','phyc') & scalePP){
          message("Normalizing PP drivers to mean=1 in first year of data")
          surf.yr1mean <- calc(out.surf[[1:12]],mean,na.rm=T)
          sum.yr1mean <-  calc(out.sum[[1:12]],mean,na.rm=T) 
          
          out.surf.norm <- out.surf/cellStats(surf.yr1mean,mean,na.rm=T)
          out.sum.norm <- out.sum/cellStats(sum.yr1mean,mean,na.rm=T)
        }
        
        #write ascii------------------------------------------------------------
        if(nlayers(out.surf)>0){
          dir.out = file.path(dir.stdriver,paste0(gsub("_glorys","",nc.vars[v]),"_surf"))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(round(out.surf,4),file.path(dir.out,paste0(gsub("_glorys","",nc.vars[v]),'_surf')),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
          rm(out.surf);gc()
        }
        
        if(nlayers(out.bott)>0){
          dir.out = file.path(dir.stdriver,paste0(gsub("_glorys","",nc.vars[v]),"_bott"))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(round(out.bott,4),file.path(dir.out,paste0(gsub("_glorys","",nc.vars[v]),'_bott')),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
          rm(out.bott);gc()
        }
        
        if(nlayers(out.mean)>0){
          dir.out = file.path(dir.stdriver,paste0(gsub("_glorys","",nc.vars[v]),"_mean"))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(round(out.mean,4),file.path(dir.out,paste0(gsub("_glorys","",nc.vars[v]),'_mean')),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
          rm(out.mean);gc()
        }
        
        if(nlayers(out.sum)>0){
          dir.out = file.path(dir.stdriver,paste0(gsub("_glorys","",nc.vars[v]),'_sum'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(round(out.sum,4),file.path(dir.out,paste0(gsub("_glorys","",nc.vars[v]),"_sum")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T, options = c("DECIMAL_PRECISION=4"))
          rm(out.sum); gc()
        }
        
        # if(nlayers(chlos.zint)>0){
        #   dir.out = file.path(dir.stdriver,paste0(gsub("_glorys","",nc.vars[v]),"_surf_Zint"))
        #   if(!dir.exists(dir.out)) dir.create(dir.out)
        #   writeRaster(round(chlos.zint,4),file.path(dir.out,paste0(gsub("_glorys","",nc.vars[v]),'_surf_Zint')),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
        #   rm(chlos.zint);gc()
        # }
        
        if(nlayers(out.sum.norm)>0 & scalePP){
          dir.out = file.path(dir.stdriver,paste0(gsub("_glorys","",nc.vars[v]),'_sum_norm'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(round(out.sum.norm,4),file.path(dir.out,paste0(gsub("_glorys","",nc.vars[v]),"_sum_norm")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T, options = c("DECIMAL_PRECISION=4"))
          rm(out.sum.norm); gc()
        }
        
        if(nlayers(out.surf.norm)>0 & scalePP){
          dir.out = file.path(dir.stdriver,paste0(gsub("_glorys","",nc.vars[v]),'_surf_norm'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(round(out.surf.norm,4),file.path(dir.out,paste0(gsub("_glorys","",nc.vars[v]),"_surf_norm")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T, options = c("DECIMAL_PRECISION=4"))
          rm(out.surf.norm); gc()
        }
        
        # if(nlayers(chlos.zint.norm)>0 & scalePP){
        #   dir.out = file.path(dir.stdriver,paste0(gsub("_glorys","",nc.vars[v]),'_surf_Zint_norm'))
        #   if(!dir.exists(dir.out)) dir.create(dir.out)
        #   writeRaster(round(chlos.zint.norm,4),file.path(dir.out,paste0(gsub("_glorys","",nc.vars[v]),"_surf_Zint_norm")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T, options = c("DECIMAL_PRECISION=4"))
        #   rm(chlos.zint.norm); gc()
        # }
      }
    }
  }
  }
  #make static maps----
  if(make.static | make.climatology){
  message("Making static and monthly climatology maps")
  dirs.stvars = list.dirs(path=dir.stdriver,recursive=F)
  dir.static = file.path(dir.stdriver,"static")
  if(!dir.exists(dir.static)) dir.create(dir.static)
  dirs.stvars = dirs.stvars[-c(which(grepl("static",dirs.stvars)),which(grepl("plots",dirs.stvars)))]
  basename(dirs.stvars)
  for(i in 1:length(dirs.stvars)){
    #i=1
    stvar.i <- unlist(strsplit(basename(dirs.stvars[i]),"_"))[1]
    if(!is.null(do.vars) & !stvar.i %in% c(do.vars,gsub("_glor","",do.vars))) next
    stname = basename(dirs.stvars[i]) #paste0(basename(dirname(dirs.stvars[i])),"_",basename(dirs.stvars[i]))
    message(paste0("now making ",i," of ",length(dirs.stvars),": ",stname))
    
    #read output--------------------------------------------------------------------
    #list ecospace output ascii files
    files.st = list.files(dirs.stvars[i],full.names = T,recursive=T)#[-c(1:6)]
    files.st = files.st[!grepl("_0000_",basename(files.st))]
    
    #read them into raster stack
    st.stack = stack(files.st)
    #names(st.stack) = paste0(rep(month.abb,24),rep(1997:2020,each=12))
    yrmo = substr(basename(files.st),nchar(basename(files.st))-9,nchar(basename(files.st))-4)
    names(st.stack) = paste0(month.abb[as.numeric(substr(yrmo,5,6))],substr(yrmo,1,4))
    
    #create static map s
    st.mean_yr1 = calc(st.stack[[1:12]],fun=mean,na.rm=T)
    st.mean_allyrs = mean(st.stack,na.rm=T)
    #st.mean_allyrs = calc(st.stack,fun=mean,na.rm=T)
    #st.mean_yr1 = mean(rast(st.stack[[1:12]]), na.rm=T)
    
    #save static map
    writeRaster(st.mean_allyrs,file.path(dir.static,paste0(stname,"_avg_allyrs")),format='ascii',overwrite=T)
    writeRaster(st.mean_yr1,file.path(dir.static,paste0(stname,"_avg_yr1")),format='ascii',overwrite=T)
    
    #no need to normalize, since the monthly maps were already normalized and in their own folder.
    # if(substr(basename(stname),1,3) %in% c('chl','npp','phy')){
    #   st.mean_yr1.norm = st.mean_yr1/cellStats(st.mean_yr1,mean,na.rm=T)
    #   st.mean_allyrs.norm = st.mean_allyrs/cellStats(st.mean_allyrs,mean,na.rm=T)
    #   
    #   writeRaster(st.mean_allyrs.norm,file.path(dir.static,paste0(stname,"_avg_allyrs_norm")),format='ascii',overwrite=T)
    #   writeRaster(st.mean_yr1.norm,file.path(dir.static,paste0(stname,"_avg_yr1_norm")),format='ascii',overwrite=T)
    # }
    
    #monthly climatology
    month_indices <- match(substr(names(st.stack),1,3),month.abb)
    monthly_averages <- stackApply(st.stack, indices = month_indices, fun = mean, na.rm = TRUE)

    writeRaster(monthly_averages,file.path(dirs.stvars[i],paste0(basename(dirs.stvars[i]),'_0000')),
                suffix=formatC(1:12,width=2,flag=0), format='ascii',overwrite=T, bylayer=T)
    
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

fn.netcdf2ascii_cefi <- function(nc.files=list.files(path=dir.cefi, pattern=".nc$", full.names=T), depth=depth,
                                 dir.stdriver=dir.stdriver, scalePP=TRUE, do.vars=NULL, 
                                 make.monthly=TRUE, make.static=TRUE, make.climatology=TRUE){
  
  if(make.monthly){
  for(i in 1:length(nc.files)){
    #i=3
    nc <- nc_open(nc.files[i])
    nc.vars = names(nc$var)
    message(nc.vars)
    ndims = nc$ndims
    message("- Processing variable ",i," of ",length(nc.files),": ",nc.vars)
    
    #Get coordinate variables----
    lon <- ncvar_get(nc, "lon")
    lat <- ncvar_get(nc, "lat")
    zvals <- nc$dim$z_l$vals
    zthickness <- c(diff(zvals-zvals[1]),tail(diff(zvals-zvals[1]),1))
    #zthickness <- diff(zvals-zvals[1])
    #zthickness <- zvals
    time <- ncvar_get(nc,'time')
    ntime <- length(time)
    refdate <- as.Date(gsub("[A-Za-z ]", "",nc$dim$time$units))
    nc.times <- seq.Date(from=refdate,by='month',length.out=ntime)
    # Find spatial indices for bbox
    lon_idx <- which(lon >= bbox[1] & lon <= bbox[2])
    lat_idx <- which(lat >= bbox[3] & lat <= bbox[4])
    cat("  Spatial subset: lon", range(lon[lon_idx]), "lat", range(lat[lat_idx]), "\n")

    # Extract data for spatial subset
    if(ndims==4){
    nc.v <- ncvar_get(nc, nc.vars, 
                         start = c(min(lon_idx), min(lat_idx), 1,1),
                         count = c(length(lon_idx), length(lat_idx), -1,-1))
    nc.v <- aperm(nc.v,c(2,1,3,4))
    
    } else{
    nc.v <- ncvar_get(nc, nc.vars, 
                         start = c(min(lon_idx), min(lat_idx), 1),
                         count = c(length(lon_idx), length(lat_idx), -1))
    nc.v <- aperm(nc.v,c(2,1,3))
    }
    nc_close(nc)
    
    if(make.monthly){
      ndims = length(dim(nc.v))
  
      if(ndims==3){ 
        ntimes = dim(nc.v)[3]
        dimnames(nc.v)[[3]] <- gsub("-","_",substr(nc.times,1,10))
        #nc.tmp = aperm(nc.v,c(2,1,3))
        nc.tmp = nc.v[dim(nc.v)[1]:1,,]
        brick1 = brick(nc.tmp, crs=crs(depth))
        extent(brick1) = extent(depth)
        ras.v <- resample(rast(brick1),rast(depth))
        #plot(ras.v)
        out.v = stack()
        for(t in 1:nlyr(ras.v)){
          #t=1
          message(paste0(nc.vars,"--month ",t," of ",ntimes))
          ras.t.filled <- fill_coastal_cells(ras.v[[t]])
          ras.t <- mask(ras.t.filled,rast(depth))
          out.v <- addLayer(out.v, raster(ras.t))
          rm(ras.t, ras.t.filled); gc()
        }
        names(out.v) <- paste0("X",format(nc.times,"%Y%m"))[1:dim(out.v)[3]]
        
        #----unit conversions----
        #convert units of NPP from mol N m-2 s-1 to gC m-2 month-1
        if(nc.vars=='wc_vert_int_npp'){
          message("Converting units of NPP from mol N m-2 s-1 to gC m-2 month-1")
          #mol N m-2 s-1 * C:N ratio of 106:16 * molar weight of C 12 g * seconds in month
          out.v <- out.v*6.625*12*30.5*24*60*60
        }
        
        #convert units of chlos from kg m-3 to mg m-3
        if(nc.vars=='chlos'){
          message("Converting units of chlos from kg m-3 to mg m-3")
          #mol N m-2 s-1 * C:N ratio of 106:16 * molar weight of C 12 g * seconds in month
          out.v <- out.v*1e6
        }
        
        #convert units of chlo from ug/kg to mg m-3
        if(nc.vars=='chl'){
          message("Converting units of chl from ug kg-1 to mg m-3")
          #multiply by seawater density 1.025 kg/L
          out.v <- out.v*1.025
        }
        
        #convert units of oxygen from mol kg-1 to  mmol m-3
        if(nc.vars=='btm_o2'){
          message("Converting units of btm_o2 from mol kg-1 to  mmol m-3")
          #mol N m-2 s-1 * C:N ratio of 106:16 * molar weight of C 12 g * seconds in month
          out.v <- out.v*1000*1025
        }
        
        if(nlayers(out.v)>0) out.v <- round(out.v,4)
        if(nlayers(out.v)>0) names(out.v) <- paste0("X",format(nc.times,"%Y%m"))[1:dim(out.v)[3]]
        
        #normalize PP to mean of first year
        if(nc.vars %in% c('chlos','wc_vert_int_npp') & scalePP){
          message("Normalizing PP drivers to mean=1 in first year of data")
          surf.yr1mean <-  calc(out.v[[1:12]],mean,na.rm=T) 
          out.v.norm <- out.v/cellStats(surf.yr1mean,mean,na.rm=T)

          dir.out = file.path(dir.stdriver,paste0(nc.vars,'_norm'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(round(out.v.norm,4),file.path(dir.out,paste0(nc.vars,"_norm")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T, options = c("DECIMAL_PRECISION=4"))
          rm(out.v.norm); gc()
          
        }
        
        
        #write ascii
        dir.out = file.path(dir.stdriver,nc.vars)
        if(!dir.exists(dir.out)) dir.create(dir.out)
        writeRaster(out.v,file.path(dir.out,nc.vars),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T)
        rm(out.v);gc()
      }
      
      if(ndims==4){
        ntimes = dim(nc.v)[4]
        dimnames(nc.v)[[4]] <- gsub("-","_",substr(nc.times,1,10))
        out.sum = out.sum.norm = stack() 
        
        #loop over timesteps----
        for(t in 1:ntimes){  #big monthly loop takes time, need to parallel
          #t=1
          print(paste0(nc.vars,"--month ",t," of ",ntimes));flush.console()
          brick1 = brick(nc.v[,,,t],crs=crs(depth))
          dim(brick1)
          extent(brick1) = extent(depth)
          brick1 <- flip(brick1, direction='y')
          nlayers(brick1)
  
          #sum over water column----
          if(nc.vars %in% c('chl')){
            #ras.sum = sum(rast(brick1),na.rm=T)
            ras.sum <- app(rast(brick1), fun = function(v) {
                           if (all(is.na(v))) {
                             NA               # keep land as NA
                           } else {
                             #sum(v * zthickness, na.rm = TRUE) #old way, entire watercolumn
                             # Calculate cumulative depth for each layer
                             #cumul_depth <- cumsum(zthickness)
                             # Find layers within first 50m
                             within_50m <- zvals <= 50
                             # For layers that partially extend beyond 50m, 
                             # adjust thickness to only include the portion within 50m
                             adjusted_thickness <- zthickness
                             if(any(zvals > 50 & c(0, zvals[-length(zvals)]) < 50)) {
                               partial_layer <- which(zvals > 50 & c(0, zvals[-length(zvals)]) < 50)
                               adjusted_thickness[partial_layer] <- 50 - c(0, zvals[-length(zvals)])[partial_layer]
                             }
                             
                             # Sum only the layers within 50m
                             sum(v[within_50m] * adjusted_thickness[within_50m], na.rm = TRUE)
                           }
                         })
            #plot(ras.sum)
            ras.sum = resample(ras.sum, rast(depth))
            ras.sum.filled <- fill_coastal_cells(ras.sum, mask=depth)
            ras.sum <- mask(ras.sum.filled, rast(depth))
            out.sum <- addLayer(out.sum, raster(ras.sum))
            rm(ras.sum, ras.sum.filled); gc()
          }
        }
        rm(brick1);gc()
        
        #----unit conversions----
        #convert units of chl from ug/kg to mg m-3
        if(nc.vars=='chl'){
          message("Converting units of chl from ug kg-1 to mg m-3")
          #multiply by seawater density 1.025 kg/L
          out.sum <- out.sum*1.025
        }
        
        if(nlayers(out.sum)>0) out.sum <- round(out.sum,4)
        if(nlayers(out.sum)>0) names(out.sum) <- paste0("X",format(nc.times,"%Y%m"))
        
        #normalize PP to mean of first year
        if(nc.vars %in% c('chl') & scalePP){
          message("Normalizing PP drivers to mean=1 in first year of data")
          sum.yr1mean <-  calc(out.sum[[1:12]],mean,na.rm=T) 
          out.sum.norm <- out.sum/cellStats(sum.yr1mean,mean,na.rm=T)
        }
        
        #write ascii------------------------------------------------------------
        if(nlayers(out.sum)>0){
          dir.out = file.path(dir.stdriver,paste0(nc.vars,'_sum'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(round(out.sum,4),file.path(dir.out,paste0(nc.vars,"_sum")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T, options = c("DECIMAL_PRECISION=4"))
          rm(out.sum); gc()
        }
 
        if(nlayers(out.sum.norm)>0 & scalePP){
          dir.out = file.path(dir.stdriver,paste0(nc.vars,'_sum_norm'))
          if(!dir.exists(dir.out)) dir.create(dir.out)
          writeRaster(round(out.sum.norm,4),file.path(dir.out,paste0(nc.vars,"_sum_norm")),bylayer=T,suffix=format(nc.times[1:ntimes],"%Y%m"),format='ascii',overwrite=T, options = c("DECIMAL_PRECISION=4"))
          rm(out.sum.norm); gc()
        }
      }
          }
        }
      }
      #make static maps----
      if(make.static | make.climatology){
        message("Making static and monthly climatology maps")
        dirs.stvars = list.dirs(path=dir.stdriver,recursive=F)
        dir.static = file.path(dir.stdriver,"static")
        if(!dir.exists(dir.static)) dir.create(dir.static)
        if(length(c(which(grepl("static",dirs.stvars)),which(grepl("plots",dirs.stvars)))>0)){
        dirs.stvars = dirs.stvars[-c(which(grepl("static",dirs.stvars)),which(grepl("plots",dirs.stvars)))]
        }
        basename(dirs.stvars)
        for(i in 1:length(dirs.stvars)){
          #i=1
          stvar.i <- unlist(strsplit(basename(dirs.stvars[i]),"_"))[1]
          if(!is.null(do.vars) & !stvar.i %in% c(do.vars,gsub("_glor","",do.vars))) next
          stname = basename(dirs.stvars[i]) #paste0(basename(dirname(dirs.stvars[i])),"_",basename(dirs.stvars[i]))
          message(paste0("now making ",i," of ",length(dirs.stvars),": ",stname))
          
          #read output--------------------------------------------------------------------
          #list ecospace output ascii files
          files.st = list.files(dirs.stvars[i],full.names = T,recursive=T)#[-c(1:6)]
          files.st = files.st[!grepl("_0000_",basename(files.st))]
          
          #read them into raster stack
          st.stack = stack(files.st)
          #names(st.stack) = paste0(rep(month.abb,24),rep(1997:2020,each=12))
          yrmo = substr(basename(files.st),nchar(basename(files.st))-9,nchar(basename(files.st))-4)
          names(st.stack) = paste0(month.abb[as.numeric(substr(yrmo,5,6))],substr(yrmo,1,4))
          
          #create static map s
          st.mean_yr1 = calc(st.stack[[1:12]],fun=mean,na.rm=T)
          st.mean_allyrs = mean(st.stack,na.rm=T)
          #st.mean_yr1 = mean(rast(st.stack[[1:12]]), na.rm=T)
          #st.mean_allyrs = calc(st.stack,fun=mean,na.rm=T)
          
          
          #save static map
          writeRaster(st.mean_allyrs,file.path(dir.static,paste0(stname,"_avg_allyrs")),format='ascii',overwrite=T)
          writeRaster(st.mean_yr1,file.path(dir.static,paste0(stname,"_avg_yr1")),format='ascii',overwrite=T)
          
          #no need to normalize, since the monthly maps were already normalized and in their own folder.
          # if(substr(basename(stname),1,3) %in% c('chl','npp','phy')){
          #   st.mean_yr1.norm = st.mean_yr1/cellStats(st.mean_yr1,mean,na.rm=T)
          #   st.mean_allyrs.norm = st.mean_allyrs/cellStats(st.mean_allyrs,mean,na.rm=T)
          #   
          #   writeRaster(st.mean_allyrs.norm,file.path(dir.static,paste0(stname,"_avg_allyrs_norm")),format='ascii',overwrite=T)
          #   writeRaster(st.mean_yr1.norm,file.path(dir.static,paste0(stname,"_avg_yr1_norm")),format='ascii',overwrite=T)
          # }
          
          #monthly climatology
          month_indices <- match(substr(names(st.stack),1,3),month.abb)
          monthly_averages <- stackApply(st.stack, indices = month_indices, fun = mean, na.rm = TRUE)
          
          writeRaster(monthly_averages,file.path(dirs.stvars[i],paste0(basename(dirs.stvars[i]),'_0000')),
                      suffix=formatC(1:12,width=2,flag=0), format='ascii',overwrite=T, bylayer=T)
          
          rm(st.stack); gc()
        }
      }
}#eof
    
    
    
    



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' Download and Process ESA CCI Ocean Data
#'
#' Downloads and processes Ocean Colour, Sea Surface Temperature, and Sea Surface Salinity
#' data from ESA Climate Change Initiative (CCI) datasets. Creates spatially and temporally
#' subsetted NetCDF files for use in ocean modeling applications.
#'
#' @param dir.out Character. Output directory path for downloaded netcdf files. Default: dir.esacci
#' @param bbox Numeric vector of length 4. Bounding box as c(lon_min, lat_min, lon_max, lat_max)
#' @param startdate Date. Start date for data extraction. Default: "1980-01-01"
#' @param enddate Date. End date for data extraction. Default: current system date
#' @param oc.url Character. Base URL for Ocean Colour CCI data via OPeNDAP
#' @param cds_user Character. Copernicus Climate Data Store username
#' @param cds_key Character. Copernicus Climate Data Store API key
#' @param sst.output_file Character. Base filename for SST output (without extension)
#' @param sss.url Character. Complete URL for Sea Surface Salinity CCI data via OPeNDAP
#'
#' @return Invisible NULL. Function creates NetCDF files as side effects
#'
#' @details
#' This function processes three ESA CCI datasets:
#' \itemize{
#'   \item Ocean Colour: Monthly chlorophyll-a (1997-present)
#'   \item Sea Surface Temperature: Monthly SST via Copernicus CDS (any years)
#'   \item Sea Surface Salinity: Bi-monthly data aggregated to monthly (2010-2023)
#' }
#'
#' @examples
#' \dontrun{
#' # Download data for Gulf of Mexico region
#' bbox <- c(-98, 18, -80, 31)
#' fn.pull_esacci(
#'   bbox = bbox,
#'   startdate = as.Date("2010-01-01"),
#'   enddate = as.Date("2020-12-31")
#' )
#' }
#'
#' @export


fn.pull_esacci <- function(dir.out = NULL,
    bbox = NULL,
    startdate = as.Date("1980-01-01"),
    enddate = Sys.Date(),
    oc.url = "https://www.oceancolour.org/thredds/dodsC/cci/v6.0-release/geographic/monthly/chlor_a/",
    sss.url = "https://data.cci.ceda.ac.uk/thredds/dodsC/esacci.SEASURFACESALINITY.15-days.L4.SSS.multi-sensor.multi-platform.GLOBAL-MERGED_OI_Monthly_CENTRED_15Day_0-25deg.5-5.r1",
    cds_user = NULL,
    cds_key = NULL,
    sst.output_file = 'esacci.SST.L4.combined.3_0.monthly',
    chl.output_file = 'OC_CCI_chlora',
    sss.output_file = 'esacci.SSS.L4.combined.5_5.monthly'){
  # dir.out = dir.esacci
  # bbox = bbox
  # startdate = as.Date("1980-01-01")
  # enddate = Sys.Date()
  # oc.url = "https://www.oceancolour.org/thredds/dodsC/cci/v6.0-release/geographic/monthly/chlor_a/"
  # cds_user = 'David Chagaris'
  # cds_key = '1519c994-451b-4247-b0ba-44896a3d38a2'
  # sst.output_file <- file.path(dir.out,'esacci.SST.L4.combined.3_0.monthly')
  # sss.url = "https://data.cci.ceda.ac.uk/thredds/dodsC/esacci.SEASURFACESALINITY.15-days.L4.SSS.multi-sensor.multi-platform.GLOBAL-MERGED_OI_Monthly_CENTRED_15Day_0-25deg.5-5.r1"
  
  # Setup function for CDS credentials (call this once)
  setup_cds_credentials(user = cds_user, key = cds_key)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(dir.out)) dir.create(dir.out)
  
  ##OC-CCI-----------------------------------------------------
  if(!is.null(chl.output_file)){
  # Initialize storage
  all_data <- list()
  
  year_start <- max(year(startdate),1997)
  year_end <- year(enddate)
  
  ###Download data year, month loop----
  for(year in year_start:year_end){
    for(month in 1:12){
      #year=1997;month=09
      month_str <- sprintf("%02d",month)
      
      filename <- sprintf("ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-%d%s-fv6.0.nc",year,month_str)  
      file_url <- paste0(oc.url,year,"/",filename)
      tryCatch({
        cat("  Downloading:", filename, "\n")
        
        # Open NetCDF via OPeNDAP
        nc <- nc_open(file_url)
        
        # Get coordinates (only once)
        if (length(all_data) == 0) {
          lon <- ncvar_get(nc, "lon")
          lat <- ncvar_get(nc, "lat")
          
          # Find spatial indices for bbox
          lon_idx <- which(lon >= bbox[1] & lon <= bbox[2])
          lat_idx <- which(lat >= bbox[3] & lat <= bbox[4])
          
          cat("  Spatial subset: lon", range(lon[lon_idx]), "lat", range(lat[lat_idx]), "\n")
        }
        
        # Get time
        time <- ncvar_get(nc, "time")
        
        # Extract chlorophyll-a data for spatial subset
        chlor_a <- ncvar_get(nc, "chlor_a", 
                             start = c(min(lon_idx), min(lat_idx), 1),
                             count = c(length(lon_idx), length(lat_idx), 1))
        
        # Store data
        time_key <- paste(year, month_str, sep = "-")
        all_data[[time_key]] <- list(
          chlor_a = chlor_a,
          time = time,
          year = year,
          month = month
        )
        nc_close(nc)
        
      }, error = function(e) {
        cat("  Error with", filename, ":", e$message, "\n")
      })
      
    } #end month loop
  } #end year loop
  
  ###save netcdf----
  if (length(all_data) > 0) {
    # Create combined NetCDF file
    output_file <- file.path(dir.out,paste0(chl.output_file, sprintf("_%s_%s_subset.nc", names(all_data)[1], names(all_data)[length(all_data)])))
    
    # Prepare data arrays
    n_times <- length(all_data)
    subset_lon <- lon[lon_idx]
    subset_lat <- lat[lat_idx]
    
    # Combined data array — shape MUST match the (lon, lat, time) order of the
    # variable declared below; ncvar_put does not enforce shape and will silently
    # blast the buffer in column-major order, scrambling the grid otherwise.
    combined_chlora <- array(NA, dim = c(length(subset_lon), length(subset_lat), n_times))
    time_values <- numeric(n_times)

    for (i in seq_along(all_data)) {
      combined_chlora[,,i] <- all_data[[i]]$chlor_a   # ncvar_get already returned [lon, lat]
      time_values[i] <- all_data[[i]]$time
    }
    
    # Create NetCDF file
    lon_dim <- ncdim_def("longitude", "degrees_east", subset_lon)
    lat_dim <- ncdim_def("latitude", "degrees_north", subset_lat)
    time_dim <- ncdim_def("time", "days since 1970-01-01", time_values)
    
    chlor_var <- ncvar_def("chlor_a", "mg m^-3", 
                           list(lon_dim, lat_dim, time_dim), 
                           -999,
                           longname = "Chlorophyll-a concentration")
    
    nc_out <- nc_create(output_file, list(chlor_var))
    ncvar_put(nc_out, chlor_var, combined_chlora)
    
    # Add attributes
    ncatt_put(nc_out, "chlor_a", "long_name", "Chlorophyll-a concentration")
    ncatt_put(nc_out, "chlor_a", "standard_name", "mass_concentration_of_chlorophyll_a_in_sea_water")
    ncatt_put(nc_out, 0, "title", "Ocean Colour CCI Chlorophyll-a Monthly Data")
    ncatt_put(nc_out, 0, "source", "ESA Ocean Colour Climate Change Initiative")
    ncatt_put(nc_out, 0, "bbox", paste(bbox, collapse = ","))
    ncatt_put(nc_out, 0, "creation_date", as.character(Sys.time()))
    
    nc_close(nc_out)
    
    cat("Combined file saved:", output_file, "\n")
  } else {
    cat("No OC data downloaded\n")
    return(NULL)
  }
  }
#SST-CCI--------------------------------------------------------------------------------------------
  if(!is.null(sst.output_file)){
  #monthly SST data from CCI are provided by compernicus climate data store
  start_year <- year(startdate)
  end_year <- year(enddate)
  ###Authenticate user----
  #Set up CDS API credentials if provided
  if (!is.null(cds_user) && !is.null(cds_key)) {
    wf_set_key(user = cds_user, key = cds_key)
  }
  
  # Check if credentials are set
  if (is.null(wf_get_key(user=cds_user))) {
    stop("CDS API credentials not found. Please set them using wf_set_key() or provide cds_user and cds_key arguments.")
  }
  ###Download data
  # Format years and months as character vectors
  years <- as.character(start_year:end_year)
  months <- 1:12
  months_formatted <- sprintf("%02d", months)
  
  #setup temporary folder to hold single month netcdf files
  temp_folder <- file.path(dir.esacci,'sst_temp')
  #temp_folder <- "C:\\sst_temp"
  if(!dir.exists(temp_folder)) dir.create(temp_folder)
  
  ### Download monthly netcdf in 1-year chunks----
  #loop through years
  for(y in 1:length(years)){
    #y=1
    # Define the request using the correct dataset and parameters
    request <- list(
      dataset_short_name = "satellite-sea-surface-temperature",
      variable = "sst",
      processinglevel = "level_4",
      sensor_on_satellite = "combined_product",
      version = "3_0",
      temporal_resolution = "monthly",
      year = years[y],
      month = months_formatted,
      area = c(bbox[4],bbox[1],bbox[3],bbox[2]),
      target = sst.output_file
    )
    # Submit and download the request
    cat("Submitting request to CDS API...\n")
    cat("Dataset: satellite-sea-surface-temperature\n")
    cat("Years:", paste(years[y], collapse = " to "), "\n")
    cat("Months:", paste(months_formatted, collapse = ", "), "\n")
    
    result <- wf_request(
      request = request,
      user = cds_user,
      transfer = TRUE,
      path = temp_folder
    )

    unzip(result, exdir=temp_folder)
    unlink(result)
    if(y %in% seq(10,length(years),10)) Sys.sleep(20)
  } #end year loop
  
  ### combine and save single netcdf----
  # Get list of NetCDF files
  nc_files <- list.files(temp_folder, pattern = "\\.nc$", full.names = TRUE)
  nc_files <- sort(nc_files)  # Sort to ensure chronological order
  cat("Found", length(nc_files), "NetCDF files\n")
  
  # Initialize storage
  all_sst <- list()
  all_times <- c()
  
  # Process each file
  for (i in seq_along(nc_files)) {
    cat("Processing file", i, "of", length(nc_files), "\n")
    nc <- nc_open(nc_files[i])
    
    # Get coordinates from first file
    if (i == 1) {
      lon <- ncvar_get(nc, "lon")
      lat <- ncvar_get(nc, "lat")
      
      # Get variable info
      var_names <- names(nc$var)
      sst_var <- var_names[grepl("sst|temperature", var_names, ignore.case = TRUE)][1]
      if (is.na(sst_var)) sst_var <- "analysed_sst"  # fallback
      
      cat("Using SST variable:", sst_var, "\n")
    }
    
    # Extract SST data
    sst_data <- ncvar_get(nc, sst_var)
    
    # Get time
    time_var <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")$value
    
    # Store data — keep source [lon, lat] orientation; lat coord saved below is
    # in the same source order, so labels and pixels stay aligned.
    all_sst[[i]] <- sst_data - 273.15  # convert K -> degC
    all_times[i] <- time_var

    nc_close(nc)
  }

  # Combine into 3D array. Buffer shape MUST match the (lon, lat, time) order of
  # the variable declared below; ncvar_put does not enforce shape and will
  # silently blast the buffer in column-major order, scrambling the grid.
  dims <- dim(all_sst[[1]])  # [lon, lat]
  combined_sst <- array(NA, dim = c(dims[1], dims[2], length(all_sst)))
  for (i in seq_along(all_sst)) {
    combined_sst[,,i] <- all_sst[[i]]
  }
  
  # Create output NetCDF file
  # Define dimensions
  lon_dim <- ncdim_def("longitude", "degrees_east", lon)
  lat_dim <- ncdim_def("latitude", "degrees_north", lat)
  time_dim <- ncdim_def("time", time_units, all_times)
  
  # Define variable
  sst_var_def <- ncvar_def("sst", "celsius", 
                           list(lon_dim, lat_dim, time_dim), 
                           -999,
                           longname = "Sea Surface Temperature")
  
  # Create file
  outfile <- paste0(sst.output_file,"_",substr(basename(nc_files)[1],1,6),"_",substr(basename(nc_files)[length(nc_files)],1,6),".nc")
  nc_out <- nc_create(outfile, list(sst_var_def))
  
  # Write data
  ncvar_put(nc_out, sst_var_def, combined_sst)
  
  # Add attributes
  ncatt_put(nc_out, "sst", "long_name", "Sea Surface Temperature")
  ncatt_put(nc_out, "sst", "standard_name", "sea_surface_temperature")
  ncatt_put(nc_out, 0, "title", "ESA CCI SST Monthly Data")
  ncatt_put(nc_out, 0, "source", "ESA Climate Change Initiative")
  ncatt_put(nc_out, 0, "creation_date", as.character(Sys.time()))
  
  nc_close(nc_out)
  
  # Clean up
  unlink(temp_folder, recursive = TRUE, force=T)
  
  cat("Combined NetCDF file saved to:", outfile, "\n")
  
  
  # Optional: return some basic info about the downloaded file
  if (file.exists(outfile)) {
    nc <- nc_open(outfile)
    cat("File dimensions:\n")
    for (dim_name in names(nc$dim)) {
      cat("  ", dim_name, ":", nc$dim[[dim_name]]$len, "\n")
    }
    cat("Variables:\n")
    for (var_name in names(nc$var)) {
      cat("  ", var_name, "\n")
    }
    nc_close(nc)
  }
  }
  #SSS-CCI--------------------------------------------------------------------------------------------
  if(!is.null(sss.output_file)){
  startyear <- max(year(startdate),2010)
  endyear <- min(year(enddate), 2023)
  
  # OPeNDAP URL
  #base_url <- paste0(sss.url,"?lat[0:1:719],lon[0:1:1439],time[0:1:334],sss[0:1:334][0:1:719][0:1:1439]")
  ### get coordinate data----
  base_url <- paste0(sss.url,"?lat,lon,time,sss")
  
  # Open the NetCDF file
  nc <- nc_open(base_url)
  
  # Get coordinate variables
  lons <- ncvar_get(nc, "lon")
  lats <- ncvar_get(nc, "lat")
  times <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  
  # Convert time to dates
  time_origin <- as.Date(sub("days since ", "", time_units))
  dates <- time_origin + times
  
  # Find spatial indices
  lon_idx <- which(lons >= bbox[1] & lons <= bbox[2])
  lat_idx <- which(lats >= bbox[3] & lats <= bbox[4])
  
  lon_vals <- lons[lon_idx]
  lat_vals <- lats[lat_idx]
  
  # Find temporal indices for the specified years
  start_date <- as.Date(paste0(startyear, "-01-01"))
  end_date <- as.Date(paste0(endyear, "-12-31"))
  time_idx <- which(dates >= start_date & dates <= end_date)
  date_vals <- dates[time_idx]
  
  
  ### Extract salinity data----
  cat("Extracting salinity data...\n")
  sss_data <- ncvar_get(nc, "sss", 
                        start = c(min(lon_idx), min(lat_idx), min(time_idx)),
                        count = c(length(lon_idx), length(lat_idx), length(time_idx)))
  
  sss_data <- aperm(sss_data,c(2,1,3))
  sss_data <- sss_data[dim(sss_data)[1]:1,,]
  dimnames(sss_data) <- list(lats[lat_idx],lons[lon_idx],as.character(dates[time_idx]))
  
  # average by month
  month.grps <- sort((substr(dimnames(sss_data)[[3]],1,7)))
  sss_mon <- apply(sss_data, c(1,2), function(x) tapply(x, month.grps, mean, na.rm=TRUE))
  sss_mon <- aperm(sss_mon,c(2,3,1))
  
  # Close NetCDF connection
  nc_close(nc)
  
  ###save as netcdf----
  # Create dimensions
  month_vals <- dimnames(sss_mon)[[3]]
  
  # Convert time strings to numeric (days since origin)
  month_dates <- as.Date(paste0(month_vals, "-15"))  # Use 15th of each month
  time_numeric <- as.numeric(month_dates - as.Date("1970-01-01"))
  
  # Define dimensions
  lat_dim <- ncdim_def("lat", "degrees_north", lat_vals)
  lon_dim <- ncdim_def("lon", "degrees_east", lon_vals)
  time_dim <- ncdim_def("time", "days since 1970-01-01 00:00:00 UTC", time_numeric, unlim=TRUE)
  
  # Define the SSS variable
  sss_var <- ncvar_def("sss", "pss", list(lat_dim, lon_dim, time_dim), 
                       missval = -999,
                       longname = "Sea Surface Salinity",
                       prec = "float")
  
  # Create the NetCDF file
  output_file <- file.path(dir.out,paste0(sss.output_file,".",month_vals[1],"-",month_vals[length(month_vals)],'.nc'))
  nc_out <- nc_create(output_file, list(sss_var), force_v4=TRUE)
  
  # Write the data (need to transpose back to lon,lat,time order)
  #sss_out <- aperm(sss_mon, c(2,1,3))
  ncvar_put(nc_out, sss_var, sss_mon)
  
  # Add global attributes
  ncatt_put(nc_out, 0, "title", "Monthly Sea Surface Salinity from ESA CCI")
  ncatt_put(nc_out, 0, "source", "ESA CCI Sea Surface Salinity v5.5")
  ncatt_put(nc_out, 0, "created", as.character(Sys.time()))
  
  # Close the file
  nc_close(nc_out)
  
  cat("Saved monthly SSS data to", output_file, "\n")
  }
}#eof




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
  month.id = substr(basename(files.st),nchar(basename(files.st))-5,nchar(basename(files.st))-4)
  #read them into raster stack
  st.stack = stack(files.st)
  #names(st.stack) = paste0(rep(month.abb,24),rep(1997:2020,each=12))
  yrmo = substr(basename(files.st),nchar(basename(files.st))-9,nchar(basename(files.st))-4)
  names(st.stack) = paste0(month.abb[as.numeric(substr(yrmo,5,6))],substr(yrmo,1,4))
  
  monthly_means <- stackApply(st.stack, indices = month.id, fun = mean, na.rm = TRUE)
  
  writeRaster(monthly_means,file.path(dir.stdriver,gsub("_cefi","",basename(dir.stdriver))),bylayer=T,suffix=paste0("9999",sort(unique(month.id))),format='ascii',overwrite=T)

  rm(st.stack); gc()
}





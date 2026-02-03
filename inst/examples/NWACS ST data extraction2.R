rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
library('raster')
source("C:\\Users\\dchagaris\\Github\\R4EwE\\R\\ST_data_extraction.R")
getwd()

#SETUP-----------------------------------------------------------------------
maxdepth = 5000
depth.30min = raster(file.path(dirname(getwd()),"maps","depth","depth 26x35 30min.asc"))
depth.15min = raster(file.path(dirname(getwd()),"maps","depth","depth 52x70 15min.asc"))
depth.5min = raster(file.path(dirname(getwd()),"maps","depth","depth 156x210 5min.asc"))
#bbox = c(-80,-62.5,33,46) #(W,E,S,N)
bbox = as.numeric(as.matrix(extent(depth.15min)))[c(1,3,2,4)]

#GLORYS-------------------------------------------------------------------------
#there is a CopernicusMarine package, but I couldn't get it working for salinity, so still doing it the other CLI way.
startdate = "1993-01-01"
enddate = Sys.Date()
dir.glorys = paste0(getwd(),"/GLORYS")  #set directory to save netcdf file.
path.copernicusmarine <- "C:/Users/dchagaris/OneDrive - University of Florida/copernicusmarine/copernicusmarine.exe"

#pull the GLORYS data and create netcdf files.  This can take awhile for many dates and large grids
fn.pull_glorys(dir.out=dir.glorys, bbox=bbox, startdate=startdate, enddate=enddate, 
               path.copernicusmarine=path.copernicusmarine,
               datid.phy = "cmems_mod_glo_phy-all_my_0.25deg_P1M-m",
               datid.bgc = "cmems_mod_glo_bgc_my_0.25deg_P1M-m",
               vars.phy = c('thetao_glor','so_glor'),
               vars.bgc = c('chl','nppv','o2','phyc'))

#save asciis
files.nc = list.files(dir.glorys, pattern=".nc$", full.names=T)
basename(files.nc)
files.nc = files.nc[c(1,3)]
fn.netcdf2ascii_glorys(nc.files = files.nc, depth=depth.15min,
                      dir.stdriver=file.path(dirname(getwd()),'ST drivers','15min'),
                      make.monthly=TRUE,
                      make.static=TRUE)

#CEFI-------------------------------------------------------------------------
startdate = "1993-01-01"
enddate = Sys.Date()
dir.cefi = paste0(getwd(),"/CEFI")  #set directory to save raw cefi data pull
if(!dir.exists(dir.cefi)) dir.create(dir.cefi)

vars.cefi <- fn.get_cefi_vars(dir.cefi=dir.cefi, file.cefivars='CEFI_vars.csv', region='nwa')

##subset for variables of interest----
do.vars = c('btm_o2','btm_temp','sos','tob','tos','wc_vert_int_npp')
vars.cefi <- subset(vars.cefi, cefi_output_frequency=='monthly' & 
                      cefi_grid_type=='regrid' & 
                      cefi_experiment_type=='hindcast' &
                      cefi_experiment_name=='nwa12_cobalt_v2' &
                      cefi_variable %in% do.vars)

#download the full netcdf files
fn.download_cefi_full_nc(vars.cefi=vars.cefi, dir.cefi=dir.cefi, reg='northwest_atlantic')

#save asciis
files.nc = list.files(dir.cefi, pattern=".nc$", full.names=T)
basename(files.nc)
#files.nc = files.nc[c(1,3)]  #subset
fn.netcdf2ascii_cefi(nc.files = files.nc, depth=depth.15min,
                       dir.stdriver=file.path(dirname(getwd()),'ST drivers','15min'),
                       make.monthly=TRUE,
                       make.static=TRUE)



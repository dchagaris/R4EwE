rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
library('raster')
library('terra')
library('ncdf4')
source("C:\\Users\\dchagaris\\Github\\R4EwE\\R\\ST_data_extraction.R")
getwd()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#SETUP-----------------------------------------------------------------------
maxdepth = 500
res <- 15

## Directories ----
wd = getwd()
dir.stdriver <- file.path(dirname(wd), "WFS EwE", "Ecospace", "ST drivers",paste0(res,'min'))
dir.depth <- file.path(dirname(dirname(dir.stdriver)), "maps/bathymetry/")
dir.cefi <- file.path(wd,"CEFI")
dir.modis <- file.path(wd, "MODIS")
dir.glorys <- file.path(wd, "GLORYS")
path.copernicusmarine <- "C:/Users/dchagaris/OneDrive - University of Florida/copernicusmarine/copernicusmarine.exe"
for(d in c(dir.modis,dir.glorys,dir.cefi,dir.stdriver)){if(!dir.exists(d)) dir.create(d)}

## Depth map ----
file.depth <- list.files(dir.depth,pattern=".asc$",full.names=T)
file.depth <- file.depth[grep(paste0(" ",res,'min'),basename(file.depth))]
file.depth <- file.depth[-grep('excl',basename(file.depth))]
file.excl = gsub("/depth ","/excl layer ",file.depth)

#get depth raster
depth <- raster(file.depth)
plot(depth, colNA='black')

#geographic extent
bbox <- matrix(extent(depth))[,1]

startdate = "1993-01-01"
enddate = Sys.Date()

#Pull CEFI----------------------------------------------------------------------
#full domain netcdf exist below, need to move that to stand alone cefi folder
#dir.cefi <- "C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\ST data extraction\\CEFI"
if(!dir.exists(dir.cefi)) dir.create(dir.cefi)

vars.cefi <- fn.get_cefi_vars(dir.cefi=dir.cefi, file.cefivars='CEFI_vars.csv', region='nwa')

##subset for variables of interest----
do.vars = c('btm_o2','btm_temp','sos','tob','tos','wc_vert_int_npp','chlos')
do.vars = 'chl'
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
fn.netcdf2ascii_cefi(nc.files = files.nc[1], depth=depth,
                     dir.stdriver=dir.stdriver,
                     make.monthly=TRUE,
                     make.static=TRUE)

#make static maps
dir.stdriver
dirs.stvars = list.dirs(dir.stdriver, recursive=F, full.names=T)
dirs.stvars = dirs.stvars[!grepl('hoard|plots|static',basename(dirs.stvars))]
dirs.stvars = dirs.stvars[2]
for(i in 1:length(dirs.stvars)){
  #i=1
  message(paste0('making static map ',i,' of ',length(dirs.stvars), ': ',basename(dirs.stvars[1])))
  fn.make_static_maps(dir.stdriver=dirs.stvars[i])
} #eof

##monthly climatology maps----
for(i in 1:length(dirs.stvars)){
  fn.make_monthly_climatology_maps(dir.stdriver=dirs.stvars[i])
}


#Pull GLORYS--------------------------------------------------------------------
##pull from Copernicus----
# variable    units
# chl         mg m-3
# phyc        mmol C m-3
# nppv        mg C m-3 day-1

#need to list available variables
fn.pull_glorys(dir.out=dir.glorys, 
               bbox=bbox,
               startdate = startdate,
               enddate = enddate,
               datid.phy = "cmems_mod_glo_phy-all_my_0.25deg_P1M-m",
               datid.bgc = "cmems_mod_glo_bgc_my_0.25deg_P1M-m",
               vars.phy = c('thetao','so','bottomT'),
               vars.bgc = c('chl','nppv','o2','phyc'),
               path.copernicusmarine = path.copernicusmarine,
               getbathygrid=FALSE)

##write to ascii----
files.nc = list.files(dir.glorys, pattern=".nc$", full.names=T)
basename(files.nc)
fn.netcdf2ascii_glorys(nc.files = files.nc, depth=depth,
                       dir.stdriver=dir.stdriver,
                       make.monthly=TRUE,
                       make.static=TRUE)

#Pull pre-GLORYS CORA-----------------------------------------------------------
#Pull MODIS--------------------------------------------------------------------
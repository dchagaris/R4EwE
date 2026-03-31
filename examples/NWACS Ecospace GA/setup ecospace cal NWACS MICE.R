# rm(list=ls())
library('dplyr'); library('reshape2'); library('tidyr'); library('stringr'); library('foreach')
library('PBSmodelling'); library('snowfall'); library('parallel'); library('snow');
library('doSNOW'); library('openxlsx'); library('tcltk')
library('plyr')
library('digest')
library('raster')
library('terra')

#source("C:/Users/dchagaris/GitHub/R-EwE-console/R_EwE_console_functions.R")
#source("C:/Users/dchagaris/GitHub/EcospaceCal/R/R_Ecospace_GA_calibration.R")
#source("C:/Users/dchagaris/GitHub/EcospaceCal/R/R_Ecospace_calibration_functions_cl.R")

#source("C:/dchagaris/GitHub/R-EwE-console/R_EwE_console_functions.R")
#source("C:/dchagaris/GitHub/EcospaceCal/R/R_Ecospace_calibration_functions_cl.R")
#source("C:/dchagaris/GitHub/EcospaceCal/R/R_Ecospace_GA_calibration.R")

source("C:\\dchagaris\\GitHub\\R4EwE\\R\\Ecospace_objective_functions.R")
source("C:\\dchagaris\\GitHub\\R4EwE\\R\\R_EwE_console_functions.R")
source("C:\\dchagaris\\GitHub\\R4EwE\\R\\input_output.R")
source("C:\\dchagaris\\GitHub\\R4EwE\\R\\utils.R")
source("C:\\dchagaris\\GitHub\\R4EwE\\R\\plotting.R")
source("C:\\dchagaris\\GitHub\\R4EwE\\R\\R_Ecospace_GA_calibration.R")
source("C:\\dchagaris\\GitHub\\R4EwE\\R\\fn.GA.R")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#USER INPUT---------------------------------------------------------------------
##directories----
dir.main <- getwd()
dir.ewe <- dirname(dir.main)
dir.out <- paste0("C:/NWACS MICE/GA output ",Sys.Date())
dir.runs <- list.dirs(dir.out, recursive=F)
dir.ga_results <- file.path(dir.main,'ga_results')
dir.obsmaps <- file.path(dir.main,'obs spatial data')

file.console <- "C:/Program Files/EcoSpace Console 1.4.4.0 64-bit/EwEClientConsole.exe"
file.obsts <- paste0(dir.ewe,"/ts1_1985_2023_annual_effort forcing.csv")

#30 min: start with base run pars
file.parsens <- paste0(dir.ewe,"/sensitivity/runlistLL 2026-02-19_30min.rds")
file.cmdbase <- paste0(dir.main,"/NWACS MICE CommandFile_template_30min.txt")
file.basevul <- paste0(dir.ewe,"/NWACS MICE v3.3-Vulnerabilities_base.csv")
file.basedisp <- paste0(dir.ewe,"/NWACS MICE v3.3-Dispersal_base.csv")
file.envpars <- paste0(dir.ewe,"/NWACS MICE v3.3-Functional responses pars.csv")
file.fleetdyn <- paste0(dir.ewe,"/NWACS MICE v3.3-Fleet dynamics_base.csv")


#15 min: start with 30 min estimates
# file.parsens <- paste0(dir.ewe,"/sensitivity/runlistLL 2026-02-04_15min.rds")
# file.cmdbase <- paste0(dir.main,"/NWACS MICE CommandFile_template_15min.txt")
# file.basevul <- paste0(dir.main,"/vuls_ga_final_20260203_182641_30min.csv")
# file.basedisp <- paste0(dir.main,"/disp_ga_final_20260203_182641_30min.csv")

# #5 min: start with 15 min estimates
# file.parsens <- paste0(dir.ewe,"/sensitivity/runlistLL 2026-02-04_5min.rds")
# file.cmdbase <- paste0(dir.main,"/NWACS MICE CommandFile_template_5min.txt")
# file.basevul <- paste0(dir.main,"/vuls_ga_final_20260204_132813_15min.csv")
# file.basedisp <- paste0(dir.main,"/disp_ga_final_20260204_132813_15min.csv")


#file.basevul <- paste0(dir.ewe,"/NWACS MICE v3.3-Vulnerabilities_gafit_15min.csv")
#file.basevul <- paste0(dir.main,"/vuls_ga_final_20260115_182355_30min.csv")
#file.basedisp <- paste0(dir.ewe,"/NWACS MICE v3.3-Dispersal_gafit_15min.csv")
#file.basedisp <- paste0(dir.main,"/disp_ga_final_20260115_182355_30min.csv")


#clear contents of output folders
del.oldout <- F  #deletes output from previous days' runs.
del.out1 <- T
del.out2 <- T
del.runs <- T

##years----
#years to run at each phase (shorter years for sensitivity test to speed up)
styear <- 1985
enyear <- 2023
nyrs <- enyear - styear + 1

##n parameters to estimate----
n_parsens = 34

##types of parameters to estimate----
est_par_types <- c('vul','disp','med','env','flt')
est_par_types <- est_par_types[c(1,5)]

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#PROCESSING---------------------------------------------------------------------
##setup folder structure-----
if(!dir.exists(dir.main)) dir.create(dir.main)
if(!dir.exists(dirname(dir.out))) dir.create(dirname(dir.out))
if(!dir.exists(dir.out)) dir.create(dir.out)
dir.oldout <- list.dirs(dirname(dir.out),full.names=T, recursive=F)
if(del.oldout) unlink(dir.oldout[-which(dir.oldout==dir.out)], recursive=T)
if(del.runs) unlink(dir.runs, recursive=T)

##get base command file and vulnerability starting values----
cmd_base <- as.character(read.delim(file.cmdbase, header=F, blank.lines.skip = T)[,1])
#vuls.base <- read.csv(file.basevul,row.names = 1)
vuls.base <- read.csv(file.basevul,row.names=1)
#dispersal rates
disp.base <- read.csv(file.basedisp, row.names=1)#[,1:5]
if('disp_final' %in% names(disp.base)){
  disp.base$disp = disp.base$disp_final
  disp.base = disp.base[,1:5]
} else{
  disp.base <- disp.base[,1:5]
}
names(disp.base) <- c('group','disp','ibm_hab_mult','rel_vul_mult','fitness_mult')

##Standardize FG names----
group.names <- gsub(" ","_",vuls.base[,1]); group.names
df.names <- data.frame(group.names = group.names, num = 1:length(group.names)) ## 2 column df with standardized name and pool code
write.csv(df.names, paste0(dir.ewe,"/FG_names_standardized.csv"), row.names = FALSE)

##Rename columns of vul matrix----
n_consumers <- ncol(vuls.base)-1; n_consumers
names(vuls.base)[1:n_consumers+1] <- group.names[1:n_consumers]
names(vuls.base)[1] <- "prey"
vuls.base[1] <- gsub(" ","_",vuls.base[,1])
vuls.base$prey <- as.factor(vuls.base$prey)

##melt vul matrix----
preprey <- fn.longvuls(vuls.base)
# predprey = as.matrix(vuls.base[,-1])
# predprey = melt(predprey, na.rm=T)
# names(predprey) = c('prey','pred','basevul')
# predprey$pred = match(predprey$pred, df.names$group.names) #as.numeric(gsub("X","",predprey$pred))
# predprey$basevul = NULL
# predprey = predprey[,2:1]

#fleet dynamics
fleetdyn.base <- read.csv(file.fleetdyn)
names(fleetdyn.base) <- c('fleet.idx','fleet','effective.power','effort.mult')

#environmental response pars
envpars.base <- read.csv(file.envpars, row.names=1)
envpars.base$shape.idx <- ifelse(envpars.base$shape=='linear',1,
                          ifelse(envpars.base$shape=='sigmoid_legacy',2,
                          ifelse(envpars.base$shape=='hyperbolic',3,
                          ifelse(envpars.base$shape=='exponential',4,
                          ifelse(envpars.base$shape=='betapdf',5,
                          ifelse(envpars.base$shape=='normal',6,
                          ifelse(envpars.base$shape=='rightshoulder',7,
                          ifelse(envpars.base$shape=='leftshoulder',8,
                          ifelse(envpars.base$shape=='trapezoid',9,
                          ifelse(envpars.base$shape=='sigmoid',10,
                          ifelse(envpars.base$shape=='logistic4params',11,0)))))))))))


# ##get most sensitive pars
# names(runlist.LL)
# if(n_parsens>0){
# runlist.LL = readRDS(file.parsens)
# names(runlist.LL)[which(names(runlist.LL)=='type')] <- 'shape'
# runlist.LL$shape.idx = 9
# #runlist.LL[,c(1:3,5,which(names(runlist.LL)=='LL.total'):ncol(runlist.LL))] = lapply(runlist.LL[,c(3:6,9,which(names(runlist.LL)=='LL.total'):ncol(runlist.LL))],as.numeric)
# runlist.LL$delta.LL = abs(runlist.LL$LL.total-runlist.LL$LL.total[runlist.LL$run.num==0])
# runlist.LL = runlist.LL[order(-runlist.LL$delta.LL),]
# names(runlist.LL)
# runlist.LL$dup = duplicated(runlist.LL[,c('prey','pred','tag.type','fleet.idx','resp.idx')])
# if(length(which(runlist.LL$dup))>0) runlist.LL = runlist.LL[-which(runlist.LL$dup),]
# parsens.topN = unique(runlist.LL[1:n_parsens,c('run.num','tag.type','prey','pred','fleet.idx','resp.idx','resp.name','shape','shape.idx','base.val','par1','par2','par3','par4','LL.total','delta.LL')])#c(4,5,1:3)])
# }


##load reference data-----
#read observed timeseries csv file input to Ecosim
obs.ts <- fn.read_ecosim_timeseries(file.obsts)
obs.ts$obsB.head$weight[obs.ts$obsB.head$weight==10] = 1
obs.ts$obsC.head$weight[obs.ts$obsC.head$weight==10] = 1

#read observed raster stacks
# files.rasters = gsub(".gri","",list.files(dir.obsmaps, pattern='.gri$', full.names=T, recursive=F))
# obs.maps <- stack(files.rasters)
# names(obs.maps)
# names(obs.maps) <- paste0(names(obs.maps),"_",rep(c(12,5,8,6,3,10),2))
# 
# strsplit(names(obs.maps),split="_")
# season = matrix(unlist(strsplit(names(obs.maps),split="_")),ncol=3,byrow=T)[,2]
# grpnum = matrix(unlist(strsplit(names(obs.maps),split="_")),ncol=3,byrow=T)[,3]
# 
# obs.maps.meta <- data.frame(layername = names(obs.maps),
#                             source = basename(files.rasters),
#                             yr.start = ifelse(basename(files.rasters)=='NEFSC sdm rasters', 2010 ,NA),
#                             yr.end = ifelse(basename(files.rasters)=='NEFSC sdm rasters', 2019 ,NA),
#                             mo.start = ifelse(season=='FALL',9,ifelse(season=='SPRING',3,NA)),
#                             mo.end = ifelse(season=='FALL',11,ifelse(season=='SPRING',5,NA)),
#                             grp.num = grpnum)


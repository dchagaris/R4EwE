#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Script to run the Ecospace GA calibration
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls()); graphics.off(); rm(.SavedPlots); windows(record=T)

# --- Setup----
#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
#remotes::install_github("dchagaris/R4EwE", dependencies=T, force=T)
#library("R4EwE")
#pkgload::load_all("C:/dchagaris/GitHub/R4EwE")
#ls("package:R4EwE")

# This file must define cmd_base, fn.runEwE, and other required variables
source("./setup ecospace cal NWACS MICE.R")

#---- Check base run ----
dir.pred = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\output\\sp01_30min_res"
#dir.pred = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\output\\sp02_15min_res"
#dir.pred = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\output\\sp03_5min_res"

predB = fn.ecospace_predB_ts2array(dir.out=dir.pred, timestep='annual',n.reg=0)
predC = fn.ecospace_predC_ts2array(dir.out=dir.pred, timestep='annual',n.reg=0)
dimnames(predB)
fn.ecospace_plot_ts(predB=predB, predC=predC, timestep='annual', obs.ts=obs.ts, scale2run=1, 
                    pltB.dims=c(3,3), pltC.dims=c(4,3), scaleCatch=T, plot2pdf=F,plt.cols=c('black','salmon','seagreen'))

#fn.objfxn2(dir.pred=dir.pred, obs.ts=obs.ts, obs.maps=obs.maps, obs.maps.meta=obs.maps.meta, autoweight.LL=TRUE)

#---- Prepare parameters ----
fit.vuls = TRUE
fit.env = TRUE
fit.disp = TRUE
fit.med = FALSE
fit.fltdyn = TRUE
# #vulnerability parameters to estimate
# predcol_vuls = data.frame(pred=1:15, prey=NA, baseval=c(36,74,714,25,4,7,5,2,5,3,11,9,2,2,2))
# 
# #environmental response parameters
#envpars <- read.csv("C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\Aquamaps\\NWACS MICE import response functions parameters.csv")
#envpars$Function.number = 1:nrow(envpars)
##get most sensitive pars
n_parsens <- 34
runlist.LL = readRDS(file.parsens)
runlist.LL$delta.LL = abs(runlist.LL$LL.total-runlist.LL$LL.total[runlist.LL$run.num==0])
runlist.sub = data.frame()
if(fit.vuls) runlist.sub = rbind(runlist.sub, runlist.LL[which(runlist.LL$tag.type=='predprey'),])
if(fit.env) runlist.sub = rbind(runlist.sub, runlist.LL[which(runlist.LL$tag.type=='environmental response'),])
if(fit.disp) runlist.sub = rbind(runlist.sub, runlist.LL[which(runlist.LL$tag.type=='dispersal'),])
if(fit.med) runlist.sub = rbind(runlist.sub, runlist.LL[which(runlist.LL$tag.type=='mediation'),])
if(fit.fltdyn) runlist.sub = rbind(runlist.sub, runlist.LL[which(runlist.LL$tag.type%in%c('effective power','effort multiplier')),])
runlist.sub = runlist.sub[order(-runlist.sub$delta.LL),]
runlist.sub$dup = duplicated(runlist.sub[,c('prey','pred','tag.type','fleet.idx','resp.idx')])
if(length(which(runlist.sub$dup))>0) runlist.sub = runlist.sub[-which(runlist.sub$dup),]
parsens.topN =  runlist.sub[1:n_parsens,c('run.num','tag.type','prey','pred','fleet.idx','resp.idx','resp.name','shape','shape.idx','base.val','par1','par2','par3','par4','LL.total','delta.LL')]



#make the parameter vector
if(fit.vuls) predprey_pars <- parsens.topN[parsens.topN$tag.type=='predprey',c('pred','prey','base.val')]
if(fit.disp) disp_pars <- parsens.topN[parsens.topN$tag.type=='dispersal',c('pred','base.val')]
if(fit.fltdyn) fltdyn_pars <- parsens.topN[parsens.topN$tag.type%in%c('effective power','effort multiplier'),c('fleet.idx','tag.type','base.val')]
if(fit.env) env_pars <- parsens.topN[parsens.topN$tag.type=='environmental response',c('resp.idx','shape.idx','par1','par2','par3','par4')] 

fn.makeparvec(do.vuls=fit.vuls, vul_pars=predprey_pars, vul.min=1.01, vul.max=1e6, vul.cv=0.4, 
              do.env=fit.env, env_pars=env_pars, env.min=0.1, env.max=2, env.cv=0.3,
              do.disp=fit.disp, disp_pars=disp_pars, disp.cv=0.2,
              do.med = fit.med, med_pars=NULL, med.xbase.cv = .1,
              do.fltdyn = fit.fltdyn, fltdyn_pars=fltdyn_pars, fltdyn.min=0.1, fltdyn.max=5, fltdyn.cv=0.3)

data.frame(L.bounds, est_par_vec, U.bounds, par_cv_vec, par.groups)

# create output directory------------------------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <<- file.path(dir.out, paste0("GA_Run_", timestamp))
dir.create(run_dir, recursive = TRUE)
message(sprintf("Output will be saved in: %s\n", run_dir))
unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)

# --- cluster setup with logging (outside GA) ---
#need to throttle back the number of workers, at least on the Titan
#8 workers took 10 minutes for 1 gen at popsize 340.
#24 workers took 5 minutes for 1 gen at popsize 340, which is about the same time as 96 workers
try(silent = TRUE, stopCluster(cl))
rm(cl); gc()
closeAllConnections()
workers <- floor(detectCores() / 4)
cl <- parallel::makePSOCKcluster(workers, outfile = "cluster_workers.log")
doParallel::registerDoParallel(cl)



myconfig <- list(popSize = 10*length(unique(par.groups)), #10x number of parameters
                 pmutation = 0.1, #1/length(log_par_vec), 
                 elitism=ceiling(.01*10*length(est_par_vec)), #1% of population size
                 n_gen=30,
                 do.penalty=FALSE,
                 pen.wt.mult=0.001,
                 gapop.pardist='uniform',
                 gapop.vuldist='uniform',
                 mutate.margin=0.2,
                 dir.ga_results=dir.ga_results)

# Export only what’s needed, once:
clusterExport(cl, c("file.console",'cmd_base','myconfig','run_dir',
  "safe_runEwE", "fn.runEwE", "fn.objfxn1", "fn.objfxn2",
  "fn.ecospace_predB_ts2array", "fn.ecospace_predC_ts2array",
  "styear", "enyear", "group.names", "df.names","obs.ts"))


#run the GA function
fn.GA(myconfig)

#save bestfit pars
ga.results = read.csv(file.path(dir.ga_results,paste0("ga_results_",timestamp,".csv")), skip=6, header=T)
bestfit_pars = data.frame(par_label = names(est_par_vec),
                          pred=ifelse(substr(names(est_par_vec),1,3)%in%c('vul','dis'),as.numeric(sapply(strsplit(names(est_par_vec),"_"),function(x)x[2])),NA),
                          prey=ifelse(substr(names(est_par_vec),1,3)=='vul',as.numeric(sapply(strsplit(names(est_par_vec),"_"),function(x)x[3])),NA),
                          envresp= ifelse(substr(names(est_par_vec),1,3)=='env',as.numeric(gsub("_","",gsub('env','',substr(names(est_par_vec),1,5)))),NA),
                          fleet= ifelse(substr(names(est_par_vec),1,3)=='flt',as.numeric(gsub("_","",gsub('flt','',substr(names(est_par_vec),1,5)))),NA),
                          
                                                    base_pars= est_par_vec,
                          final_pars = as.numeric((ga.results[nrow(ga.results),-c(1:3)])))
vuls.final = vuls.base
for(i in vul.par.idx){
  pd.i = bestfit_pars$pred[i]
  py.i = bestfit_pars$prey[i]
  vul.i = bestfit_pars$final_pars[i]
  vuls.final[py.i,pd.i+1] <- round(vul.i,2)
}
write.csv(vuls.final,file.path(dir.ga_results,paste0('vuls_ga_final_',timestamp,'.csv')))

disp.final = disp.base
disp.final$disp_final = disp.final$disp
disp.final$disp_final[bestfit_pars$pred[substr(bestfit_pars$par_label,1,4)=='disp']] <- bestfit_pars$final_pars[substr(bestfit_pars$par_label,1,4)=='disp']
write.csv(disp.final, file.path(dir.ga_results,paste0('disp_ga_final_',timestamp,'.csv')))


#explore runs-------------------------------------------------------------------
#ga.results = ga.results[1:13,]
gapop.fit <- ga.results[c(1,nrow(ga.results)),-c(1:3)]

files.cmd <- lapply(1:nrow(gapop.fit),function(i) fn.parvec2cmd(par_vec=gapop.fit[i,], g=0, idx=i)) 
files.cmd <- unlist(files.cmd)
fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1, delete.output=F)

plot.runs = c(list.dirs(run_dir,recursive=F))#,"C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\output\\sp02_15min_res")
#predB <- fn.ecospace_predB_ts2array(dir.out=list.dirs(run_dir,recursive=F), timestep='annual', n.reg=0)
#predC <- fn.ecospace_predC_ts2array(dir.out=list.dirs(run_dir,recursive=F), timestep='annual', n.reg=0)
predB <- fn.ecospace_predB_ts2array(dir.out=plot.runs, timestep='annual', n.reg=0)
predC <- fn.ecospace_predC_ts2array(dir.out=plot.runs, timestep='annual', n.reg=0)


dimnames(predB)[[3]] <- dimnames(predC)[[3]] <- substr(dimnames(predB)[[3]],1,14) #c('base','ga_fit2')#,'ga_fit_saved')#,'ga_fit_gui','15min')
# fn.ecospace_plot_ts(predB=predB, predC=predC, timestep='annual', obs.ts=obs.ts, scale2run=2,dir.plts=getwd(),
#                     pltB.dims=c(3,3), pltC.dims=c(4,3), plt.cols=c('black','salmon','seagreen','steelblue'),scaleCatch=T, plot2pdf=F)
graphics.off();rm(.SavedPlots);windows(record=T)
fn.ecospace_plot_ts(predB=predB, predC=predC, timestep='annual', obs.ts=obs.ts, scale2run=2,dir.plts=dir.ga_results,
                    pltB.dims=c(3,3), pltC.dims=c(4,3), plt.cols=1:7,scaleCatch=T, plot2pdf=T, 
                    run.label=gsub("GA_Run_","",basename(dirname(plot.runs))))


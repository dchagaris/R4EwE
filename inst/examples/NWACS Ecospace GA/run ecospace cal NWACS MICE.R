#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Script to run the Ecospace GA calibration
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls()); graphics.off(); rm(.SavedPlots); windows(record=T)
# --- Setup----
# This file must define cmd_base, fn.runEwE, and other required variables
source("./setup ecospace cal NWACS MICE.R")

#---- Check base run ----
#dir.pred = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\output\\sp01_30min_res"
#dir.pred = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\output\\sp02_15min_res"
dir.pred = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\output\\sp03_5min_res"

predB = fn.ecospace_predB_ts2array(dir.out=dir.pred, timestep='annual',n.reg=0)
predC = fn.ecospace_predC_ts2array(dir.out=dir.pred, timestep='annual',n.reg=0)
dimnames(predB)
fn.ecospace_plot_ts(predB=predB, predC=predC, timestep='annual', obs.ts=obs.ts, scale2run=1, 
                    pltB.dims=c(3,3), pltC.dims=c(4,3), scaleCatch=T, plot2pdf=F,plt.cols=c('black','salmon','seagreen'))
matrix(unlist(lapply(1:dim(predB)[3],function(i) fn.objfxn1(dir.pred=dir.pred, obs.ts=obs.ts, run.idx=i))),nrow=dim(predB)[3],byrow=T)
#fn.objfxn2(dir.pred=dir.pred, obs.ts=obs.ts, obs.maps=obs.maps, obs.maps.meta=obs.maps.meta, autoweight.LL=TRUE)

#---- Prepare parameters ----
# #vulnerability parameters to estimate
# predcol_vuls = data.frame(pred=1:15, prey=NA, baseval=c(36,74,714,25,4,7,5,2,5,3,11,9,2,2,2))
# 
# #environmental response parameters
envpars <- read.csv("C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\Aquamaps\\NWACS MICE import response functions parameters.csv")
envpars$Function.number = 1:nrow(envpars)

#make the parameter vector
predprey_pars <- parsens.topN[parsens.topN$tag.type=='predprey',c('pred','prey','base.val')]
disp_pars <- parsens.topN[parsens.topN$tag.type=='dispersal',c('pred','base.val')]

fn.makeparvec(do.vuls=TRUE, vul_pars=predprey_pars, vul.min=1.01, vul.max=1e6, vul.cv=0.4, 
              do.env=FALSE, env_pars=NULL, env.min=-1, env.max=1, env.cv=0.2,
              do.disp=TRUE, disp_pars=disp_pars, disp.cv=0.2)

cbind(L.bounds, log_par_vec, U.bounds, exp(L.bounds), exp(log_par_vec), exp(U.bounds), L.bounds*.8, U.bounds*1.2, exp(L.bounds*.5), exp(U.bounds*1.5))
cbind(exp(L.bounds), exp(log_par_vec), exp(U.bounds))


# create output directory------------------------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <<- file.path(dir.out, paste0("GA_Run_", timestamp))
dir.create(run_dir, recursive = TRUE)
message(sprintf("Output will be saved in: %s\n", run_dir))
#unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)


# --- cluster setup with logging (outside GA) ---
#need to throttle back the number of workers, at least on the Titan
#8 workers took 10 minutes for 1 gen at popsize 340.
#24 workers took 5 minutes for 1 gen at popsize 340, which is about the same time as 96
try(silent = TRUE, stopCluster(cl))
rm(cl); gc()
closeAllConnections()
workers <- floor(detectCores() / 4)
cl <- parallel::makePSOCKcluster(workers, outfile = "cluster_workers.log")
doParallel::registerDoParallel(cl)



myconfig <- list(popSize = 10*length(log_par_vec), #10x number of parameters
                 pmutation = 0.1, #1/length(log_par_vec), 
                 elitism=ceiling(.01*10*length(log_par_vec)), #1% of population size
                 n_gen=30,
                 do.penalty=FALSE,
                 pen.wt.mult=0.001)

# Export only what’s needed, once:
clusterExport(cl, c("file.console",'cmd_base','myconfig','run_dir',
  "safe_runEwE", "fn.runEwE", "fn.objfxn1", "fn.objfxn2",
  "fn.ecospace_predB_ts2array", "fn.ecospace_predC_ts2array",
  "styear", "enyear", "group.names", "df.names","obs.ts"))


#run the GA function
fn.GA(myconfig)

#explore runs-------------------------------------------------------------------
ga.results = read.csv(file.path(getwd(),"ga_results_20260116_165410.csv"), skip=6, header=T)
#ga.results = ga.results[1:13,]
gapop.fit <- log(ga.results[c(1,nrow(ga.results)),-c(1:3)])

files.cmd <- lapply(1:nrow(gapop.fit),function(i) fn.parvec2cmd(log_par_vec=gapop.fit[i,], g=0, idx=i)) 
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

fn.ecospace_plot_ts(predB=predB, predC=predC, timestep='annual', obs.ts=obs.ts, scale2run=which.min(fitness),dir.plts=getwd(),
                    pltB.dims=c(3,3), pltC.dims=c(4,3), plt.cols=1:7,scaleCatch=T, plot2pdf=T)
#c('black','salmon','seagreen')

#save bestfit pars
bestfit_pars = data.frame(par_label = names(log_par_vec),
                          pred=as.numeric(sapply(strsplit(names(log_par_vec),"_"),function(x)x[2])),
                          prey=as.numeric(sapply(strsplit(names(log_par_vec),"_"),function(x)x[3])),
                          base_log_pars= log_par_vec,
                          base_pars = exp(log_par_vec),
                          final_log_pars = as.numeric(log(ga.results[which.min(fitness),-c(1:3)])),
                          final_pars = (as.numeric(ga.results[which.min(fitness),-c(1:3)])))
vuls.final = vuls.base
for(i in vul.par.idx){
  pd.i = bestfit_pars$pred[i]
  py.i = bestfit_pars$prey[i]
  vul.i = bestfit_pars$final_pars[i]
  vuls.final[py.i,pd.i+1] <- round(vul.i,2)
}
write.csv(vuls.final,file.path(getwd(),paste0('vuls_ga_final_',timestamp,'_5min.csv')))

disp.final = disp.base
disp.final$disp_final = disp.final$disp
disp.final$disp_final[bestfit_pars$pred[substr(bestfit_pars$par_label,1,4)=='disp']] <- bestfit_pars$final_pars[substr(bestfit_pars$par_label,1,4)=='disp']
write.csv(disp.final, file.path(getwd(),paste0('disp_ga_final_',timestamp,'_5min.csv')))

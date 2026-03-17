source("C:\\dchagaris\\GitHub\\R4EwE\\R\\Ecospace_objective_functions.R")
#source("C:\\dchagaris\\GitHub\\R4EwE\\R\\R_EwE_console_functions.R")
source("C:\\dchagaris\\GitHub\\R4EwE\\R\\input_output.R")
#source("C:\\dchagaris\\GitHub\\R4EwE\\R\\utils.R")
source("C:\\dchagaris\\GitHub\\R4EwE\\R\\plotting.R")
#source("C:\\dchagaris\\GitHub\\R4EwE\\R\\R_Ecospace_GA_calibration.R")
#source("C:\\dchagaris\\GitHub\\R4EwE\\R\\fn.GA.R")


#---- Check base run ----
dir.pred = "C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\output\\sp01_30min_res"
predB = fn.ecospace_predB_ts2array(dir.out=dir.pred, timestep='annual',n.reg=0)
predC = fn.ecospace_predC_ts2array(dir.out=dir.pred, timestep='annual',n.reg=0)
fn.ecospace_plot_ts(predB=predB, predC=predC, timestep='annual', obs.ts=obs.ts, scale2run=1, 
                    pltB.dims=c(3,3), pltC.dims=c(4,3), scaleCatch=T, plot2pdf=F,plt.cols=c('black','salmon','seagreen'))

#objective fxn for multiple runs
matrix(unlist(lapply(1:dim(predB)[3],function(i) fn.objfxn1(dir.pred=dir.pred, obs.ts=obs.ts, run.idx=i))),nrow=dim(predB)[3],byrow=T)
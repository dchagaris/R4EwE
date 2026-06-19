#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.runEwE
# fn.runEwE <- function(cmdfile, do.obj = 1) {
#   #runlist=runlist_sens; i=3; do.obj=T
#   #create command line and run Ecospace
#   #dir.cmdfile <- runlist$cmd_file[1]
#   #cmdfile <- files.cmd[1]
#   cmd = paste(paste0('"', file.console, '"'),
#               paste0('"', cmdfile, '"'))
#   system(cmd,intern = F)
#   
# 
#   #calculate objective function
#   if(do.obj==1){
#     #objout <- fn.objfxn1(i,runlist)
#     objout <- fn.objfxn1(dir.pred=dirname(cmdfile),obs.ts=obs.ts)
#     return(objout)
#   }
#   if(do.obj==2){
#     objout <- fn.objfxn2(dir.pred=dirname(cmdfile), obs.ts=obs.ts, obs.map=obs.map)
#     return(objout)
#   }
# }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.runEwE-----------------------------------------------------------------------------------------
#' Run one Ecospace job (system2 version)
#' @param cmdfile Character: full path to the command file Ecospace consumes
#' @param do.obj Integer: 0=none, 1=\code{fn.objfxn1}, 2=\code{fn.objfxn2},
#' 3=\code{fn.ecospace_objfxn}
#' @param console Character: full path to the Ecospace console executable (defaults to global
#' file.console)
#' @param stdout_target NULL/FALSE to suppress, or a file path to capture stdout
#' @param stderr_target NULL/FALSE to suppress, or a file path to capture stderr
#' @param fit.abs.catch Logical, passed to \code{fn.ecospace_objfxn} when \code{do.obj == 3}.
#' If \code{TRUE} (default), absolute-type series (BiomassAbs, Catches, Landings, Discards,
#' DiscardsTotalAbs) are fit on their literal scale (\code{q = 1}); if \code{FALSE}, every
#' series is q-rescaled. Ignored for \code{do.obj} 1 or 2.
#' @param timeout Numeric seconds; passed to \code{system2}. If the Ecospace child
#' process hasn't exited within this many seconds it is killed and the run is treated as
#' a failure (returns \code{NA}). Guards against hung EwE.exe processes that would
#' otherwise block a worker indefinitely (e.g., a Windows Error Reporting dialog from
#' a child crash). Default 300 (5 minutes).
#' @return Integer exit status (0 = success); invisibly
fn.runEwE <- function(cmdfile,
                      do.obj   = 1,
                      console  = file.console,
                      stdout_target = FALSE,
                      stderr_target = FALSE,
                      fit.abs.catch = TRUE,
                      timeout = 300) {
  
  # cmdfile=files.cmd[1]
  # do.obj=1
  # console=file.console
  # stdout_target=file.path(dirname(cmdfile),"stdout.txt")
  # stderr_target=FALSE

  # Normalize paths (do NOT shQuote() for system2 on Windows)
  console_path <- normalizePath(console, winslash = "\\", mustWork = TRUE)
  cmdfile_path <- normalizePath(cmdfile, winslash = "\\", mustWork = TRUE)
  
  # Build argument vector: the command file is passed as an argument
  args <- shQuote(cmdfile_path)

  # Run the Ecospace console; block until it finishes (wait = TRUE) or until `timeout`
  # seconds elapse, at which point the child is killed and a non-zero status is returned.
  # The timeout guards against an EwE.exe child that hangs (e.g., on a Windows Error
  # Reporting dialog after a crash); without it a single bad parameter vector can stall
  # a worker indefinitely.
  status <- tryCatch(
    system2(command = console_path,
            args    = args,
            stdout  = stdout_target,    # FALSE to suppress, or a file path
            stderr  = stderr_target,    # FALSE to suppress, or a file path
            wait    = TRUE,
            timeout = timeout),
    error = function(e) {
      # Could not launch the process (bad path, permissions, etc.)
      warning(sprintf("Ecospace call failed to start: %s", conditionMessage(e)))
      999L
    }
  )
  
  
  # # If only running the external job, return its status invisibly
  # if (isTRUE(do.obj == 0L)) {
  #   if (!identical(status, 0L)) {
  #     warning(sprintf("Ecospace returned non-zero exit status: %s", status))
  #   }
  #   return(invisible(status))
  # }
  # 
  # # If we are computing an objective, do it only for successful runs
  # if (!identical(status, 0L)) {
  #   warning(sprintf("Ecospace returned non-zero exit status: %s", status))
  #   return(NA_real_)
  # }
  # 
  # 
  # # Optionally warn if non-zero exit status
  #  if (!is.numeric(status) || length(status) != 1L || status != 0L) {
  #    warning(sprintf("Ecospace returned non-zero exit status: %s", status))
  #  }

  # Compute objective (must return a single numeric or NA_real_)
  pred_dir <- dirname(cmdfile_path)
  
  out <- tryCatch({
    if (do.obj == 1L) {
      fn.objfxn1(dir.pred = pred_dir, obs.ts = obs.ts)
    } else if (do.obj == 2L) {
      fn.objfxn2(dir.pred = pred_dir, obs.ts = obs.ts, obs.map = obs.map)
    } else if (do.obj == 3L) {
      fn.ecospace_objfxn(dir.pred = pred_dir, obs.ts = obs.ts,
                         group.names = group.names, fleet.names = fleet.names,
                         fit.abs.catch = fit.abs.catch,
                         return = "vector")
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  
  if (is.numeric(out) && length(out) >= 1L && is.finite(out[1])){
    return(out)
  } else {
    return(NA_real_)
  }
  
} #eof


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.runEwE.parallel-----------------------------------------------------------------------------------------
#' @title Run EwE CLI app in parallel
#' @description Calls fn.runEwE in parallel framework.
#' @param runlist A dataframe that must contain a variable named 'cmd_file', which is the full file 
#' path of the command files to run.
#' @param obj.fxn A single number indicating the likelihood function to use. 0=none, 1=timeseries
#' only (\code{fn.objfxn1}), 2=timeseries and spatial data (\code{fn.objfxn2}), 3=Ecospace
#' timeseries with fleet-specific landings and discards (\code{fn.ecospace_objfxn}).
#' \code{obj.fxn = 3} additionally requires \code{group.names} and \code{fleet.names} to be in
#' the worker's environment (add them to \code{cl.export}).
#' @param cl.export A list of objects exported to each worker.
#' @param fit.abs.catch Logical, passed to \code{fn.ecospace_objfxn} when \code{obj.fxn == 3}.
#' If \code{TRUE} (default), absolute-type series are fit on their literal scale; if
#' \code{FALSE}, every series is q-rescaled. Ignored for other \code{obj.fxn} values.
#' @return Model output is saved according to command file. If do.obj!=0 then a vector likelihoods
#' is returned and added to the runlist dataframe.
#' @examples
#' # example code:
#' \dontrun{result <- fn.runEwE.parallel(runlist=myrunlist, obj.fxn=1, cl.export=list('myrunlist','obs.ts'))}
#' @export
fn.runEwE.parallel <-  function(
    runlist=runlist,
    obj.fxn=1,
    cl.export = list("runlist", "obs.ts"),
    fit.abs.catch = TRUE,
    timeout = 1800,
    capture_log = TRUE,
    max_retries = 5,
    n_cores = NULL
){
  t1 <- Sys.time()
  if(is.null(n_cores)) n_cores <- max(1L, detectCores()-1L)
  cl <- makeSOCKcluster(n_cores)
  registerDoSNOW(cl)
  clusterExport(cl, append(cl.export, list(
    "file.console", "fn.runEwE",
    "fn.objfxn1", "fn.objfxn2", "fn.ecospace_objfxn",
    "fn.read_pred_ecospace_wide", "fn.read_pred_ecospace_discards_split", "fn.fg_meta",
    ".find_year_skip", ".detect_ecospace_regions", ".ecospace_dmort_by_year",
    "styear", "enyear", "group.names", "fleet.names", "df.names")))
  print(paste('Setup',n_cores,'Clusters: Overhead time',round(as.numeric(Sys.time()-t1),2)))

  pbar <- winProgressBar("Running Ecospace Console",label=paste0("Simulation 0 of ",nrow(runlist)),max=100)
  prog <- function(n) setWinProgressBar(pbar,(n/nrow(runlist)*100),label=paste("Simulation Run", n,"of", nrow(runlist),"Completed"))
  opts <- list(progress=prog)

  print(paste('Running',nrow(runlist),'Ecospace simulations'))
  t1 <- Sys.time()
  runs <- foreach(i = 1:nrow(runlist), .errorhandling = 'pass', .options.snow = opts) %dopar% {
    stdo <- if(capture_log) file.path(dirname(runlist$cmd_file[i]), "stdout.txt") else FALSE
    stde <- if(capture_log) file.path(dirname(runlist$cmd_file[i]), "stderr.txt") else FALSE
    fn.runEwE(cmdfile = runlist$cmd_file[i], do.obj = obj.fxn,
              stdout_target = stdo, stderr_target = stde,
              fit.abs.catch = fit.abs.catch, timeout = timeout)
  }
  close(pbar)
  print(paste('Run time',round(as.numeric(Sys.time()-t1),2)))

  ##missing runs----
  # Treat a run as failed if its dir contains nothing but the cmd file +/- log files
  # (i.e., no Ecospace CSV output). Just counting files <=1 misses cases where the
  # stdout/stderr logs exist but EwE didn't write any CSVs.
  has_ecospace_output <- function(d){
    f <- list.files(d, pattern = "^Ecospace_.*\\.csv$", ignore.case = TRUE)
    length(f) > 0L
  }
  filecheck <- sapply(runlist$dir.out, has_ecospace_output)
  erruns <- which(!filecheck)

  tries <- 0L
  while(length(erruns)>=1 && tries < max_retries){
    tries <- tries + 1L
    message(paste0('Redo missing runs (try ', tries, '/', max_retries, '): n=', length(erruns)))

    pbar <- winProgressBar("Running Ecospace Console: Missing Runs",label=paste0("Simulation 0 of ",length(erruns)),max=100)
    prog <- function(n) setWinProgressBar(pbar,(n/length(erruns)*100),label=paste("Simulation Run", n,"of", length(erruns),"Completed"))
    opts <- list(progress=prog)

    runs.erruns <- foreach(i=1:length(erruns),.errorhandling='pass',.options.snow=opts) %dopar% {
      stdo <- if(capture_log) file.path(dirname(runlist$cmd_file[erruns[i]]), "stdout.txt") else FALSE
      stde <- if(capture_log) file.path(dirname(runlist$cmd_file[erruns[i]]), "stderr.txt") else FALSE
      fn.runEwE(cmdfile=runlist$cmd_file[erruns[i]], do.obj=obj.fxn,
                stdout_target = stdo, stderr_target = stde,
                fit.abs.catch=fit.abs.catch, timeout = timeout)
    }
    close(pbar)

    for(k in 1:length(erruns)) runs[[erruns[k]]] <- unlist(runs.erruns[k])

    filecheck <- sapply(runlist$dir.out, has_ecospace_output)
    erruns <- which(!filecheck)
  }
  if(length(erruns) > 0L){
    warning(sprintf("fn.runEwE.parallel: %d run(s) still missing Ecospace output after %d retries: %s",
                    length(erruns), tries, paste(erruns, collapse=", ")))
  }
  
  runlist.out = cbind(runlist, do.call(rbind, runs))
  print('All runs completed')
  stopCluster(cl);
  closeAllConnections()
  return(runlist.out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
#fn.vul_testval_tags-----
#' @title Make vulnerability parameter taglines for sensitivity runs.
#' @description Create parameter taglines for the command file to test sensitivity of vulnerability 
#' parameters.
#' @param predprey A dataframe containing 3 numeric columns, "pred", "prey", and "basevul", 
#' containing the group numbers of predator prey pairs to evaluate and starting values for the 
#' vulnerabilities.  Missing (NA) in the prey column will set the vulnerability by predator column.
#' @param test.type Either 'absolute' in which case maxvul, minvul, and maxvul.mult are applied, or
#' 'percentage', where pct.change and (1+pct.change) is multiplied by the base value.
#' @param pct.change The percent change to use when test.type='percentage'.  Applied as pct.change 
#' x base and (1+pct.change) x base. Default is 0.5.
#' @param maxvul A high vulnerability number to test
#' @param minvul A low vulnerability number to test.
#' @param is.maxvul.mult Logical indicating whether the maxvul is treated as a multiplier on the 
#' basevul (TRUE) or as the absolute value (FALSE).
#' @return A character vector containing the parameter taglines for the command file.
#' @examples
#' # example code:
#' \dontrun{tags.vul <- fn.vul_testval_tags(predprey=data.frame(pred=1:5,prey=6:10, basevul=2), maxvul=1000, minvul=1.01, is.maxvul.mult=T)}
#' @export
fn.vul_testval_tags = function(predprey, pct.change=0.5, maxvul=100, minvul=1.01, is.maxvul.mult=T, test.type='absolute'){
  #predprey = data.frame(pred=1:5,prey=6:10, basevul=2)
  #predprey = predprey.vuls
  if(test.type=='absolute'){
  predprey.low = predprey.hi = predprey
  predprey.low$vul = minvul
  if(is.maxvul.mult) predprey.hi$vul = maxvul*predprey$basevul
  if(!is.maxvul.mult) predprey.hi$vul = maxvul
  predprey = rbind(predprey.low,predprey.hi)
  predprey = predprey[order(predprey$pred,predprey$prey),]
  tags <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", predprey$pred, " ", ifelse(is.na(predprey$prey),"",predprey$prey),
                     "), ", sprintf("%.2f", predprey$vul), ", Indexed.Single")
  return(tags)
  }
  
  if(test.type=='percentage'){
    predprey.low = predprey.hi = predprey
    predprey.low$vul = ifelse(pct.change*predprey$basevul<1.01,1.01,pct.change*predprey$basevul)
    predprey.hi$vul = (1+pct.change)*predprey$basevul
    predprey.out = rbind(predprey.low,predprey.hi)
    predprey.out = predprey.out[order(predprey.out$pred,predprey.out$prey),]
    tags <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", predprey.out$pred, " ", ifelse(is.na(predprey.out$prey),"",predprey.out$prey),
                   "), ", sprintf("%.2f", predprey.out$vul), ", Indexed.Single")
    return(tags)
  }
} #end of function

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
#fn.disp_testval_tags-----
#' @title Make dispersal parameter taglines for sensitivity runs.
#' @description Create parameter lines for the command file to test sensitivity of dispersal rates.
#' @param disp.base A dataframe with baseline dispersal parameters, exported from EwE and read and 
#' renamed in the setup file.
#' @param pct.change The percent change to use, 
#' applied as pct.change x base and (1+pct.change) x base. Default is 0.5 (+/-50%).
#' @return A character vector containing the parameter taglines for the command file.
#' @examples
#' # example code:
#' \dontrun{tags.vul <- fn.vul_testval_tags(predprey=data.frame(pred=1:5,prey=6:10, basevul=2), maxvul=1000, minvul=1.01, is.maxvul.mult=T)}
#' @export
fn.disp_testval_tags = function(disp.base, pct.change=0.5){
  disp.low = disp.hi = disp.base
  disp.low$pred = disp.hi$pred = as.numeric(row.names(disp.base))
  disp.low$disp = disp.low$disp*pct.change
  disp.hi$disp = disp.hi$disp*(1+pct.change)
  disp = rbind(disp.low[,c('pred','disp')],disp.hi[,c('pred','disp')])
  disp = disp[order(disp$pred),]
  disp$disp = round(disp$disp,2)
  tags <- paste0("<ECOSPACE_DISPERSAL_RATE_INDEXED>(",disp$pred,"), ",disp$disp,", Indexed.Single")
  return(tags)
} #end of function 

#fn.fleetdyn_testval_tags-----
#' @title Make fleet dynamics parameter taglines for sensitivity runs.
#' @description Create parameter lines for the command file to test sensitivity of effective power and effort multiplier parameters.
#' @param fleetdyn.base A dataframe with baseline parameters, exported from EwE and read 
#' in the setup file.
#' @param pct.change The percent change to use, 
#' applied as pct.change x base and (1+pct.change) x base. Default is 0.5 (+/-50%).
#' @return A character vector containing the parameter taglines for the command file.
#' @examples
#' # example code:
#' \dontrun{tags.vul <- fn.vul_testval_tags(predprey=data.frame(pred=1:5,prey=6:10, basevul=2), maxvul=1000, minvul=1.01, is.maxvul.mult=T)}
#' @export
fn_fleetdyn_testval_tags <- function(fleetdyn.base, pct.change=0.5){
  pow.low = pow.hi = fleetdyn.base
  pow.low$fleet = pow.hi$fleet = as.numeric(row.names(fleetdyn.base))
  pow.low$effective.power = pow.low$effective.power*pct.change
  pow.hi$effective.power = pow.hi$effective.power*(1+pct.change)
  pow = rbind(pow.low[,c('fleet','effective.power')],pow.hi[,c('fleet','effective.power')])
  pow = pow[order(pow$fleet),]
  pow$effective.power = round(pow$effective.power,2)
  tags.pow <- paste0("<ECOSPACE_FLEET_EFFECTIVEPOWER_INDEXED>(",pow$fleet,"), ",pow$effective.power,", Indexed.Single")
  
  mult.low = mult.hi = fleetdyn.base
  mult.low$fleet = mult.hi$fleet = as.numeric(row.names(fleetdyn.base))
  mult.low$effort.mult = mult.low$effort.mult*pct.change
  mult.hi$effort.mult = mult.hi$effort.mult*(1+pct.change)
  mult = rbind(mult.low[,c('fleet','effort.mult')],mult.hi[,c('fleet','effort.mult')])
  mult = mult[order(mult$fleet),]
  mult$effort.mult = round(mult$effort.mult,2)
  tags.mult <- paste0("<ECOSPACE_FLEET_TOTAL_EFFORTMULT_INDEXED>(",mult$fleet,"), ",mult$effort.mult,", Indexed.Single")
  
  tags <- c(tags.pow,tags.mult)
  return(tags)
} #eof

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
#fn.mediation_testval_tags-----
#' @title Make mediation shape parameter taglines for sensitivity runs.
#' @description Create parameter lines for the command file to test sensitivity to mediation shapes.
#' @param meds.base A dataframe with baseline mediation parameters, this has to be created manually.
#' @param do.xbase Logical.  Should the x_base parameter also be tested? (default=TRUE).
#' @return A character vector containing the parameter taglines for the command file.
#' @examples
#' # example code:
#' \dontrun{tags.vul <- fn.vul_testval_tags(predprey=data.frame(pred=1:5,prey=6:10, basevul=2), maxvul=1000, minvul=1.01, is.maxvul.mult=T)}
#' @export
fn.mediation_testval_tags = function(meds.base, do.xbase=TRUE){
  meds.low = meds.hi = meds.base
  meds.out = data.frame()
  #meds.low$pred = meds.hi$pred = as.numeric(meds.base$group_index)
  tags <- character()
  for(i in 1:nrow(meds.base)){
    #i=1
    #LINEAR SHAPE....
    if(meds.base$shape_type[i]==1){
      meds.low$slope[i] = meds.base$slope[i]-1.96*meds.base$stderr[i]
      meds.low$par1[i] = round(1-50*meds.low$slope[i],2)
      meds.low$par2[i] = round(100*meds.low$slope[i]+meds.low$par1[i],2)
      
      meds.hi$slope[i] = meds.base$slope[i]+1.96*meds.base$stderr[i]
      meds.hi$par1[i] = round(1-50*meds.hi$slope[i],2)
      meds.hi$par2[i] = round(100*meds.hi$slope[i]+meds.hi$par1[i],2)
      
      tag.low.i = paste("<MEDIATION_FUNCTION_INDEXED>(",meds.low$med_index[i],"),",meds.low$x_base[i],meds.low$shape_type[i], meds.low$par1[i],meds.low$par2[i],",Indexed.Single[]")
      tag.hi.i = paste("<MEDIATION_FUNCTION_INDEXED>(",meds.hi$med_index[i],"),",meds.hi$x_base[i],meds.hi$shape_type[i],meds.hi$par1[i],meds.hi$par2[i],",Indexed.Single[]")
      tags <- c(tags,tag.low.i,tag.hi.i) 
      
      meds.out <- rbind(meds.out,meds.low[i,],meds.hi[i,])
      
      if(do.xbase){
        x_base.ci = ceiling(1.96*(.017*meds.base$x_base[i]))
        meds.xbase = rbind(meds.base[i,],meds.base[i,])
        meds.xbase$x_base = c(meds.base$x_base[i]-x_base.ci, meds.base$x_base[i]+x_base.ci)
        tag.xbaselow.i = paste("<MEDIATION_FUNCTION_INDEXED>(",meds.base$med_index[i],"),",meds.base$x_base[i]-x_base.ci, meds.base$shape_type[i],meds.base$par1[i],meds.base$par2[i],",Indexed.Single[]")
        tag.xbasehi.i = paste("<MEDIATION_FUNCTION_INDEXED>(",meds.base$med_index[i],"),",meds.base$x_base[i]+x_base.ci,meds.base$shape_type[i],meds.base$par1[i],meds.base$par2[i],",Indexed.Single[]")
        tags <- c(tags,tag.xbaselow.i,tag.xbasehi.i) 
        meds.out <- rbind(meds.out,meds.xbase)
      }
    } #end linear
    
    #HYPERBOLIC SHAPE....
    if(meds.base$shape_type[i]==3){
      #low = flat curve
      meds.low$par1[i] = round(meds.base$par1[i]-1.96*meds.base$par1[i]*meds.base$stderr[i],3)
      meds.low$par2[i] = round(meds.base$par2[i]+1.96*meds.base$par2[i]*meds.base$stderr[i],3)
      meds.low$par3[i] = round(meds.base$par3[i]+1.96*meds.base$par3[i]*meds.base$stderr[i],3)
      #hi = steep curve
      meds.hi$par1[i] = round(meds.base$par1[i]+1.96*meds.base$par1[i]*meds.base$stderr[i],3) 
      meds.hi$par2[i] = round(meds.base$par2[i]-1.96*meds.base$par2[i]*meds.base$stderr[i],3)
      meds.hi$par3[i] = round(meds.base$par3[i]-1.96*meds.base$par3[i]*meds.base$stderr[i],3)

      tag.low.i = paste("<MEDIATION_FUNCTION_INDEXED>(",meds.low$med_index[i],"),",meds.low$x_base[i],meds.low$shape_type[i], meds.low$par1[i],meds.low$par2[i],meds.low$par3[i],meds.low$par4[i],",Indexed.Single[]")
      tag.hi.i = paste("<MEDIATION_FUNCTION_INDEXED>(",meds.hi$med_index[i],"),",meds.hi$x_base[i],meds.hi$shape_type[i],meds.hi$par1[i],meds.hi$par2[i],meds.hi$par3[i],meds.hi$par4[i],",Indexed.Single[]")
      tags <- c(tags,tag.low.i,tag.hi.i) 
      
      meds.out <- rbind(meds.out,meds.low[i,],meds.hi[i,])
      
      if(do.xbase){
        x_base.ci = ceiling(1.96*(.017*meds.base$x_base[i]))
        meds.xbase = rbind(meds.base[i,],meds.base[i,])
        meds.xbase$x_base = c(meds.base$x_base[i]-x_base.ci, meds.base$x_base[i]+x_base.ci)
        tag.xbaselow.i = paste("<MEDIATION_FUNCTION_INDEXED>(",meds.base$med_index[i],"),",meds.base$x_base[i]-x_base.ci, meds.base$par1[i],meds.base$par2[i],meds.base$par3[i],meds.base$par4[i],",Indexed.Single[]")
        tag.xbasehi.i = paste("<MEDIATION_FUNCTION_INDEXED>(",meds.base$med_index[i],"),",meds.base$x_base[i]+x_base.ci,meds.base$shape_type[i],meds.base$par1[i],meds.base$par2[i],meds.base$par3[i],meds.base$par4[i],",Indexed.Single[]")
        tags <- c(tags,tag.xbaselow.i,tag.xbasehi.i) 
        meds.out <- rbind(meds.out,meds.xbase)
      }
    } #end hyperbolic
  } #end loop
  meds.out$tag <- tags
  #mediation <<- meds.out
  #return(tags)
  return(meds.out)
} #eof


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
#fn.envresp_testval_tags-----
#' @title Make environmental response shape parameter taglines for sensitivity runs.
#' @description Create parameter lines for the command file to test sensitivity to environmental response shapes.
#' @param env.base A dataframe with baseline mediation parameters, this has to be created manually.
#' @return A character vector containing the parameter taglines for the command file.
#' @examples
#' # example code:
#' \dontrun{tags.vul <- fn.vul_testval_tags(predprey=data.frame(pred=1:5,prey=6:10, basevul=2), maxvul=1000, minvul=1.01, is.maxvul.mult=T)}
#' @export
fn.envresp_testval_tags = function(env.base, fact=2){
  idx_vec <- .env_resp_idx(env.base)
  env.narrow = env.wide = env.base
  env.out = data.frame()
  tags <- character()
  for(i in 1:nrow(env.base)){
    #TRAPEZOID....
    if(env.base$shape.idx[i]==9){
      pref.range = env.base$par3[i]-env.base$par2[i]
      #wide
      env.wide$par2[i] <- (env.base$par2[i]-env.base$par1[i])/fact
      env.wide$par3[i] <- env.base$par3[i]+(env.base$par4[i]-env.base$par3[i])/fact
      #narrow
      env.narrow$par1[i] <- (env.base$par2[i]-env.base$par1[i])/fact
      env.narrow$par2[i] <- env.base$par2[i]+(pref.range*(1/fact)^2)
      env.narrow$par3[i] <- env.base$par3[i]-(pref.range*(1/fact)^2)
      env.narrow$par4[i] <- env.base$par3[i]+(env.base$par4[i]-env.base$par3[i])/fact
      tag.wide.i = paste("<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(",idx_vec[i],"),",env.wide$shape.idx[i], env.wide$par1[i],env.wide$par2[i],env.wide$par3[i],env.wide$par4[i],",Indexed.Single[]")
      tag.narrow.i = paste("<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(",idx_vec[i],"),",env.narrow$shape.idx[i], env.narrow$par1[i],env.narrow$par2[i],env.narrow$par3[i],env.narrow$par4[i],",Indexed.Single[]")
      tags <- c(tags,tag.wide.i,tag.narrow.i)
      env.out <- rbind(env.out,env.base[i,],env.base[i,])
    } #end trapezoid

   } #end loop
  env.out$tag <- tags
  return(env.out)
} #eof

# Internal: resolve the env-response index column.
# Accepts either 'fxn_num' (current convention) or 'resp.idx' (legacy).
#' @keywords internal
#' @noRd
.env_resp_idx <- function(env.base){
  if("fxn_num"  %in% names(env.base)) return(env.base$fxn_num)
  if("resp.idx" %in% names(env.base)) return(env.base$resp.idx)
  stop("env.base needs a 'fxn_num' (preferred) or 'resp.idx' column for the EcoSpace response index.")
}

# Internal: normalize an env-response parameters data frame to the canonical
# schema (shape.idx + par1..par6). Accepts either the setup-script schema
# (shape.idx, par1..par6) or the sensitivity-results schema (shape, Param.1..Param.6).
# Unspecified par columns are left absent.
#' @keywords internal
#' @noRd
.normalize_env_pars <- function(env_pars){
  if(is.null(env_pars)) return(env_pars)
  out <- env_pars
  if(!"shape.idx" %in% names(out) && "shape" %in% names(out))
    out$shape.idx <- suppressWarnings(as.integer(out$shape))
  for(n in 1:6){
    canon <- paste0("par", n)
    alias <- paste0("Param.", n)
    if(!canon %in% names(out) && alias %in% names(out))
      out[[canon]] <- out[[alias]]
  }
  out
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.envresp_offval_tags-----
#' @title Make environmental-response "off" parameter taglines for sensitivity runs.
#' @description Creates one sensitivity tagline per environmental response function that
#'   neutralizes that response by replacing its shape with a flat \strong{linear} response.
#'   For foraging responses both linear endpoints are set to 1 (flat response of 1, no
#'   spatial preference). For mortality responses both linear endpoints are set to 0
#'   (flat response of 0, no extra mortality). The original shape and parameters of the
#'   response are discarded in the emitted tag; \code{par3} and \code{par4} of the linear
#'   shape are written as 0 placeholders. Use this instead of
#'   \code{fn.envresp_testval_tags} when you want to evaluate the effect of each env response
#'   by turning it off (one run per response) rather than perturbing its parameters.
#' @param env.base A data frame with one row per env response function. Must include
#'   the EcoSpace response index in either \code{fxn_num} (preferred) or \code{resp.idx}
#'   (legacy), and \code{resp_type} (character; either \code{"foraging"} or
#'   \code{"mortality"}). Other columns are passed through unchanged so the result can be
#'   \code{rbind.fill}'d into a runlist alongside other parameter types.
#' @return A data frame with one row per env response, with \code{shape.idx} set to 1
#'   (linear), \code{par1..par4} overwritten with the off values, and a \code{tag} column
#'   containing the formatted \code{<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>} tagline.
#'   Rows with unrecognized \code{resp_type} are dropped with a warning.
#' @examples
#' \dontrun{tags.off <- fn.envresp_offval_tags(envpars.base)}
#' @export
fn.envresp_offval_tags <- function(env.base){
  if(!"resp_type" %in% names(env.base))
    stop("env.base is missing required column 'resp_type'.")
  idx_vec <- .env_resp_idx(env.base)             # accepts fxn_num or resp.idx

  out <- env.base
  out$shape.idx <- 1L                            # linear "off" shape

  rt <- tolower(trimws(as.character(out$resp_type)))
  off_val <- ifelse(rt == "foraging", 1,
                    ifelse(rt == "mortality", 0, NA_real_))

  if(any(is.na(off_val))){
    bad <- unique(out$resp_type[is.na(off_val)])
    warning("Dropping ", sum(is.na(off_val)),
            " env-response row(s) with unrecognized resp_type: ",
            paste(bad, collapse = ", "),
            " (expected 'foraging' or 'mortality').")
    keep    <- !is.na(off_val)
    out     <- out[keep, , drop = FALSE]
    idx_vec <- idx_vec[keep]
    off_val <- off_val[keep]
  }

  out$par1 <- off_val
  out$par2 <- off_val
  out$par3 <- 0
  out$par4 <- 0

  # Drop any par5/par6 (or Param.5/Param.6) — linear uses only 2 params and the CLI tag
  # carries 4 slots, so anything beyond par4 is irrelevant to the off tag.
  for(col in c("par5", "par6", "Param.5", "Param.6"))
    if(col %in% names(out)) out[[col]] <- NULL

  out$tag <- paste("<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(", idx_vec, "),",
                   out$shape.idx, out$par1, out$par2, out$par3, out$par4,
                   ",Indexed.Single[]")
  rownames(out) <- NULL
  out
} #eof

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.make_cmd_files----------------------------------------------------------------------------------
#' @title Make command files.
#' @description This function adds taglines to the base command file and saves them in run folders.
#' @param runlist A dataframe where each row is a separate model run.  It must include, at a 
#' minimum, the 3 variables 'dir.out', which provides the output directory for each run; 
#' 'cmd_file' which is the full file path of the command file to be saved in the folder; and 'tag', 
#' which is the tagline to be appended to the command file. 
#' @param nyrs A single number representing the number of years to run ecospace.
#' @return Create command .txt files and saves them in their respecitive run folders.
#' @examples
#' # example code:
#' \dontrun{fn.make_cmd_files(runlist)}
#' @export
fn.make_cmd_files = function(runlist=runlist,nyrs=nyrs){
  #runlist = runlist_sens; iter=1
  for(i in 1:nrow(runlist)){
    #i=1
    cmd_i = cmd_base
    #set output directory
    param = "<ECOSPACE_OUTPUT_DIR>"
    update = paste0(runlist$dir.out[i])
    n = which(substr(cmd_i,1,nchar(param))==param)
    cmd_i[n]= paste(param, update, "System.String, Updated", sep = ", ")

    #set run length
    param = "<N_ECOSPACE_YEARS>"
    update = nyrs
    n = which(substr(cmd_i,1,nchar(param))==param)
    cmd_i[n]= paste(param, update, "System.Int32, Updated", sep = ", ")
    
    #add taglines
    #cmd_i = rbind(cmd_i,runlist$tag[i])
    cmd_i = c(cmd_i,runlist$tag[i])
    
    #save command file to run folder
    write.table(cmd_i, runlist$cmd_file[i],  row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}#eof
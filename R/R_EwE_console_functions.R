#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.runEwE-----------------------------------------------------------------------------------------
#' @title Run EwE CLI app
#' @description Function to run EwE using the CLI
#' @param cmdfile A file path pointing to the command file
#' @param do.obj Toggle swith indicating the likelihood function to use. 0=none, 1=timeseries only, 2=timeseries and spatial data
#' @return Model output is saved according to command file. If do.obj!=0 then a vector likelihoods is returned. 
#' @examples
#' # example code:
#' result <- fn.runEwE(cmdfile=C://input//EwEcmd.txt, do.obj=1)
#' @export
fn.runEwE <- function(cmdfile, do.obj = 1) {
  #runlist=runlist_sens; i=3; do.obj=T
  #create command line and run Ecospace
  #dir.cmdfile <- runlist$cmd_file[1]
  #cmdfile <- files.cmd[1]
  cmd = paste(paste0('"', file.console, '"'),
              paste0('"', cmdfile, '"'))
  system(cmd,intern = F)
  

  #calculate objective function
  if(do.obj==1){
    #objout <- fn.objfxn1(i,runlist)
    objout <- fn.objfxn1(dir.pred=dirname(cmdfile),obs.ts=obs.ts)
    return(objout)
  }
  if(do.obj==2){
    objout <- fn.objfxn2(dir.pred=dirname(cmdfile), obs.ts=obs.ts, obs.map=obs.map)
    return(objout)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.runEwE.parallel-----------------------------------------------------------------------------------------
#' @title Run EwE CLI app in parallel
#' @description Calls fn.runEwE in parallel framework.
#' @param runlist A dataframe that must contain a variable named 'cmd_file', which is the full file path of the command files to run.
#' @param obj.fxn A single number indicating the likelihood function to use. 0=none, 1=timeseries only, 2=timeseries and spatial data.
#' @param cl.export A list of objects exported to each worker.
#' @return Model output is saved according to command file. If do.obj!=0 then a vector likelihoods is returned and added to the runlist dataframe. 
#' @examples
#' # example code:
#' result <- fn.runEwE.parallel(runlist=myrunlist, obj.fxn=1, cl.export=list('myrunlist','obs.ts'))
#' @export
fn.runEwE.parallel <-  function(
    runlist=runlist, 
    obj.fxn=1, 
    cl.export = list("runlist", "obs.ts")
){
  #source(file.setup)
  t1 <- Sys.time()
  cl <- makeSOCKcluster(detectCores()-1)
  registerDoSNOW(cl)
  clusterExport(cl,append(cl.export,list("file.console", "fn.runEwE", "fn.objfxn1", "fn.objfxn2","startyear","endyear_sens","group.names","df.names")))
  print(paste('Setup',detectCores()-1,'Clusters: Overhead time',round(as.numeric(Sys.time()-t1),2)))
  
  #runlist=runlist[1:20,]
  pbar <- winProgressBar("Running Ecospace Console",label=paste0("Simulation 0 of ",nrow(runlist)),max=100)
  prog <- function(n) setWinProgressBar(pbar,(n/nrow(runlist)*100),label=paste("Simulation Run", n,"of", nrow(runlist),"Completed"))
  opts <- list(progress=prog)
  
  print(paste('Running',nrow(runlist),'Ecospace simulations'))
  t1 <- Sys.time()
  runs <- foreach(i = 1:nrow(runlist), .errorhandling = 'pass', .options.snow = opts) %dopar% {
    fn.runEwE(cmdfile = runlist$cmd_file[i], do.obj = obj.fxn)
  }
  close(pbar)
  print(paste('Run time',round(as.numeric(Sys.time()-t1),2)))
  
  ##missing runs----
  filecheck <- sapply(runlist$dir.out,FUN=function(x)length(list.files(x)))
  erruns <- which(filecheck<=1)  
  #erruns <- c(3,5,7)

  while(length(erruns)>=1){
    print(paste0('Redo missing runs: n=',length(erruns)))
    
    pbar <- winProgressBar("Running Ecospace Console: Missing Runs",label=paste0("Simulation 0 of ",length(erruns)),max=100)
    prog <- function(n) setWinProgressBar(pbar,(n/length(erruns)*100),label=paste("Simulation Run", n,"of", length(erruns),"Completed"))
    opts <- list(progress=prog)
    
    runs.erruns <- foreach(i=1:length(erruns),.errorhandling='pass',.options.snow=opts) %dopar% {
      fn.runEwE(dir.cmdfile=runlist$cmd_file[erruns[i]], do.obj=obj.fxn)
    }
    close(pbar)
    
    for(k in 1:length(erruns)) runs[[erruns[k]]] <- unlist(runs.erruns[k])
    
    filecheck <- sapply(runlist$dir.out,FUN=function(x)length(list.files(x)))
    erruns <- which(filecheck<=1)  #if there are many missing runs, then need to do this in parallel
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
#' @description Create parameter taglines for the command file to test sensitivity of vulnerability parameters.
#' @param predprey A dataframe containing 3 numeric columns, "pred", "prey", and "basevul", containing the group numbers of predator prey pairs to evaluate and starting values for the vulnerabilities.  Missing (NA) in the prey column will set the vulnerability by predator column.
#' @param test.type Either 'absolute' in which case maxvul, minvul, and maxvul.mult are applied, or 'percentage', where pct.change and (1+pct.change) is multiplied by the base value.
#' @param pct.change The percent change to use when test.type='percentage'.  Applied as pct.change x base and (1+pct.change) x base. Default is 0.5.
#' @param maxvul A high vulnerability number to test
#' @param minvul A low vulnerability number to test.
#' @param is.maxvul.mult Logical indicating whether the maxvul is treated as a multiplier on the basevul (TRUE) or as the absolute value (FALSE).
#' @return A character vector containing the parameter taglines for the command file.
#' @examples
#' # example code:
#' tags.vul <- fn.vul_testval_tags(predprey=data.frame(pred=1:5,prey=6:10, basevul=2), maxvul=1000, minvul=1.01, is.maxvul.mult=T)
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
  
  
}
  
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#fn.make_cmd_files----------------------------------------------------------------------------------
#' @title Make command files.
#' @description This function adds taglines to the base command file and saves them in run folders.
#' @param runlist A dataframe where each row is a separate model run.  It must include, at a minimum, the 3 variables 'dir.out', which provides the output directory for each run; 'cmd_file' which is the full file path of the command file to be saved in the folder; and 'tag', which is the tagline to be appended to the command file. 
#' @param nyrs A single number representing the number of years to run ecospace.
#' @return Create command .txt files and saves them in their respecitive run folders.
#' @examples
#' # example code:
#' fn.make_cmd_files(runlist)
#' @export
fn.make_cmd_files = function(runlist=runlist,nyrs=nyrs){
  #runlist = runlist_sens; iter=1
  for(i in 1:nrow(runlist)){
    #i=1
    cmd_i = cmd_base
    #set output directory
    param = "<ECOSPACE_OUTPUT_DIR>"
    update = paste0(runlist$dir.out[i])
    n = which(substr(cmd_i[,1],1,nchar(param))==param)
    cmd_i[n,1]= paste(param, update, "System.String, Updated", sep = ", ")

    #set run length
    param = "<N_ECOSPACE_YEARS>"
    update = nyrs
    n = which(substr(cmd_i[,1],1,nchar(param))==param)
    cmd_i[n,1]= paste(param, update, "System.Int32, Updated", sep = ", ")
    
    #add taglines
    cmd_i = rbind(cmd_i,runlist$tag[i])
    
    #save command file to run folder
    write.table(cmd_i, runlist$cmd_file[i],  row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
#fn.envpref_testvals-----
fn.envpref_testvals <- function(runlist) {
  # Created an empty tag column
  runlist$tag <- NA
  
  for (i in 1:nrow(runlist)) {
    if (runlist$shape.type[i] == "trapezoid") {
      runlist$tag[i] <- paste0(
        "<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(",
        runlist$responsefxn.number[i], "), ",
        runlist$shape.type.number[i], " ",
        runlist$abs.min[i], " ",
        runlist$pref.min[i], " ",
        runlist$pref.max[i], " ",
        runlist$abs.max[i],
        ", Indexed.Single[]"
      )
    } else if (runlist$shape.type[i] == "linear") {
      runlist$tag[i] <- paste0(
        "<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(",
        runlist$responsefxn.number[i], "), ",
        runlist$shape.type.number[i], " ",
        runlist$abs.min[i], " ",
        runlist$abs.max[i],
        ", Indexed.Single[]"
      )
    }
  }
  
  return(runlist)
}




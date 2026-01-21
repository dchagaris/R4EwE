
#prepare parameter vector----
fn.makeparvec <- function(
  do.vuls = TRUE, 
  vul_pars = predprey_pars,
  vul.min = 1.01,
  vul.max = 1e6,
  vul.cv = 0.3,
  do.env = FALSE, 
  env_pars = NULL,
  env.min=0.25,
  env.max=2,
  env.cv=0.2,
  do.disp = TRUE,
  disp_pars = disp_pars,
  disp.min=3,
  disp.max=3000,
  disp.cv=0.1){
  
  # do.vuls = TRUE
  # vul_pars = predprey_pars
  # vul.min = 1.01
  # vul.max = 1e6
  # vul.cv = 0.4
  # do.env = FALSE
  # env_pars = NULL
  # env.min=0.25
  # env.max=2
  # env.cv=0.2
  # do.disp = TRUE
  # disp_pars = disp_pars
  # disp.min=3
  # disp.max=3000
  # disp.cv=0.4
  
  #validate setup----
  if(do.vuls){
    pdcol = unique(vul_pars$pred[is.na(vul_pars$prey)])
    pdpy = c(0,unique(vul_pars$pred[!is.na(vul_pars$prey)]))
    if(length(which(pdcol %in% pdpy))>0){
      stop("Trying to estimate ki and kij for at least on predator. Check vul_pars input.")
    }
  }
  
  #build parameter vector----
  n_vuls <- n_env <- n_disp <- 0
  if(do.vuls) n_vuls <- nrow(vul_pars)
  if(do.env) n_env <- nrow(env_pars)  #need to melt the env resp fxn dataframe
  if(do.disp) n_disp <- nrow(disp_pars)
  #medations functions
  
  n_pars <<- n_vuls + n_env + n_disp
  if (n_pars == 0) stop("n_pars==0: No predator-prey pairs or environmental responses parameters to estimate.  Check inputs.")
  
  log_vuln_vec <- log_env_vec <- log_disp_vec <- numeric()
  par.idx <- par.labels <- character()
  if(do.vuls) {
    log_vuln_vec <- log(vul_pars$base.val)
    par.idx <- c(par.idx,rep('vul',n_vuls))
    par.labels <- c(par.labels,paste('vul',vul_pars$pred,vul_pars$prey,sep="_"))
  }
  if(do.env) {
    log_env_vec = rep(0,n_env)
    par.idx = c(par.idx,rep('env',n_env))
    par.labels = c(par.labels,paste('env',env_pars$Function.number,sep="_"))
  }
  if(do.disp){
    log_disp_vec = log(disp_pars$base.val)
    par.idx = c(par.idx,rep('disp',n_disp))
    par.labels = c(par.labels,paste('disp',disp_pars$pred,sep="_"))
  }
  log_par_vec <- c(log_vuln_vec,log_env_vec, log_disp_vec)
  names(log_par_vec) <- par.labels
  par.idx
  # index parameter types
  #need to think about how to index the parameter vector as it continues to grow
  vul.par.idx <<- which(par.idx=='vul')
  env.par.idx <<- which(par.idx=='env')
  disp.par.idx <<- which(par.idx=='disp')

  # index pred-prey vulnerabilities
  vul_pars <<- vul_pars
  disp_pars <<- disp_pars
  
  # set parameter bounds--------------------------------------------------------
  lower.vuls <- upper.vuls <- lower.env <- upper.env <- lower.disp <- upper.disp <- numeric()
  #lower.vuls <- rep(log(vul.min+1),n_vuls)
  #upper.vuls <- rep(log(vul.max+1),n_vuls)
  #lower.vuls <- vul_pars$base.val-1.96*vul_pars$base.val*vul.cv
  #upper.vuls <- log(vul_pars$base.val+1.96*vul_pars$base.val*vul.cv)
  lower.vuls <- exp(log(vul_pars$base.val)-1.96*sqrt(log(vul.cv^2+1)))
  lower.vuls <- ifelse(lower.vuls<=1.0,1.01,lower.vuls)
  lower.vuls <- log(lower.vuls)
  upper.vuls <- log(vul_pars$base.val)+1.96*sqrt(log(vul.cv^2+1))
  #cbind(lower.vuls,vul_pars$base.val,upper.vuls)
  
  lower.env <- rep(env.min,n_env)
  upper.env <- rep(env.max,n_env)
  #lower.disp <- rep(log(disp.min),n_disp)
  #upper.disp <- rep(log(disp.max),n_disp)
  #lower.disp <- log(disp_pars$base.val-1.96*disp_pars$base.val*disp.cv)
  #upper.disp <- log(disp_pars$base.val+1.96*disp_pars$base.val*disp.cv)
  lower.disp <- (log(disp_pars$base.val)-1.96*sqrt(log(disp.cv^2+1)))
  upper.disp <- (log(disp_pars$base.val)+1.96*sqrt(log(disp.cv^2+1)))
  #cbind(lower.disp,disp_pars$base.val,upper.disp)
  
  
  #super assignment
  log_par_vec <<- log_par_vec
  L.bounds <<- c(lower.vuls,lower.env,lower.disp)
  U.bounds <<- c(upper.vuls, upper.env, upper.disp)
  n_vuls <<- n_vuls
  n_env <<- n_env
  n_disp <<- n_disp
  
  #create vector of SD for rnorm() draws
  par_cv_vec <- log_par_vec
  par_cv_vec[vul.par.idx] <- vul.cv
  par_cv_vec[env.par.idx] <- env.cv
  par_cv_vec[disp.par.idx] <- disp.cv
  par_cv_vec <<- par_cv_vec
  
} #end function

#make GA populations----
fn.GApop = function(run_config=myconfig, usedist='unif'){
  #run_config <- myconfig
  if(usedist=='unif'){
    mat <- matrix(runif(run_config$popSize*n_pars, L.bounds, U.bounds), nrow=run_config$popSize, byrow=T)
  }
  if(usedist=='ln'){
    log_par_sd <- sqrt(log(par_cv_vec^2+1))
    mat <- matrix(rlnorm(run_config$popSize*n_pars, meanlog=log_par_vec, sdlog=log_par_sd), nrow=run_config$popSize, byrow=T)
    idx <- mat[, vul.par.idx] < 1.0
    mat[,vul.par.idx][idx] <- 1.01
  }
  return(log(mat))
  #matrix(rep(log_par_vec, run_config$popSize), nrow=run_config$popSize, byrow=T)
}


# Alternative minimalistic approach:
# safe_make_run_dir_tempfile <- function(run_dir, log_par_vec) {
#   h <- digest::digest(log_par_vec, algo = "md5", serialize = TRUE)
#   run_path <- tempfile(pattern = paste0("run_", h, "_"), tmpdir = run_dir)
#   dir.create(run_path, showWarnings = TRUE, recursive = TRUE)
#   run_path
# }


safe_make_run_dir_tempfile <- function(base_dir,
                                       par_vec,
                                       gen = gen,
                                       idx = idx,
                                       prefix = "run",
                                       random_len = 8) {

  
  # Create a compact, reproducible tag for the parameter vector
  # Use a hash to avoid very long names; include gen/idx if provided
  par_tag <- tryCatch({
    # stringify with limited precision to keep it small and stable
    s <- paste(round(par_vec, 6), collapse = "_")
    # hash the string; if digest not available, fallback to simple checksum
    if (requireNamespace("digest", quietly = TRUE)) {
      digest::digest(s, algo = "xxhash64")
    } else {
      # simple fallback: integer checksum
      sprintf("%08X", as.integer(sum(charToRaw(s))))
    }
  }, error = function(e) {
    sprintf("%08X", as.integer(runif(1, 0, .Machine$integer.max)))
  })
  
  # Add generation/index if provided
  gid <- if (!is.null(gen)) sprintf("g%03d", as.integer(gen)) else NULL
  iid <- if (!is.null(idx)) sprintf("i%04d", as.integer(idx)) else NULL
  
  # Add process id for extra uniqueness
  pid <- tryCatch(Sys.getpid(), error = function(e) sample.int(1e6, 1))
  
  # Random slug to avoid race collisions (tempfile-like behavior)
  rand <- paste(sample(c(letters, LETTERS, 0:9), random_len, replace = TRUE), collapse = "")
  
  # Assemble a safe folder name: prefix + tags
  parts <- c(prefix, gid, iid, paste0("pid", pid), paste0("h", par_tag), rand)
  folder_name <- paste(parts[!sapply(parts, is.null)], collapse = "_")
  
  # Sanitize (just in case)
  folder_name <- gsub("[^A-Za-z0-9_.-]", "_", folder_name)
  
  run_path <- file.path(base_dir, folder_name)
  
  # Create (retry if rare collision)
  if (!dir.exists(run_path)) {
    ok <- dir.create(run_path, recursive = TRUE, showWarnings = FALSE)
    if (!ok) {
      # one more try with a new random slug
      rand2 <- paste(sample(c(letters, LETTERS, 0:9), random_len, replace = TRUE), collapse = "")
      run_path <- file.path(base_dir, paste0(folder_name, "_", rand2))
      dir.create(run_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(run_path)
}



#write command files----
fn.parvec2cmd <- function(log_par_vec,g=gen,idx=i, out_dir=run_dir){
  #par_vec <- log_par_vec
  #log_par_vec=gapop.init[1,]
  
  #output directory----
  #run_id <- paste0("run_", digest(log_par_vec, algo=c("md5")))
  #run_path <- file.path(run_dir, run_id)
  #dir.create(run_path, showWarnings = TRUE)
  run_path <- safe_make_run_dir_tempfile(base_dir=out_dir,
                                         par_vec=log_par_vec,
                                         gen = g,
                                         idx = idx,
                                         prefix = "run",
                                         random_len = 8)
  
  #parameter tags----
  tags.vul = tags.env = tags.disp = character()
  if(length(vul.par.idx)>0){
    vuln_vec = exp(log_par_vec[vul.par.idx])
    tags.vul <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", vul_pars$pred, " ", ifelse(is.na(vul_pars$prey),"",vul_pars$prey),
                       "), ", sprintf("%.5f", vuln_vec), ", Indexed.Single")
  }

  if(length(disp.par.idx)>0){
    disp_vec = exp(log_par_vec[disp.par.idx])
    tags.disp <- paste0("<ECOSIM_DISPERSAL_RATE_INDEXED>(", disp_pars$pred,"), ",disp_vec, ", Indexed.Single")
  }
  
  if(length(env.par.idx)>0){
    env_vec = exp(log_par_vec[env.par.idx])
    respfxn_num = envpars$Function.number
    respfxn_type = envpars$Function.type
    pars1 <- pars2 <- pars3 <- pars4 <- pars5 <- numeric(length=n_env)
    for(p in 1:n_env){
      #p=1
      if(respfxn_type[p]==9){
        pars1[p] = ifelse(envpars$Param.1[p]==0 & envpars$Param.2[p]==0,0,0.5*(1-env_vec[p])*(envpars$Param.4[p]-envpars$Param.1[p])+envpars$Param.1[p])
        pars2[p] = ifelse(envpars$Param.1[p]==0 & envpars$Param.2[p]==0,0,0.5*(1-env_vec[p])*(envpars$Param.3[p]-envpars$Param.2[p])+envpars$Param.2[p])
        pars3[p] = envpars$Param.3[p]-0.5*(1-env_vec[p])*(envpars$Param.3[p]-envpars$Param.2[p])
        pars4[p] = envpars$Param.4[p]-0.5*(1-env_vec[p])*(envpars$Param.4[p]-envpars$Param.1[p])
      }
    }
    env_pars = paste(respfxn_type,sprintf("%.2f",pars1), sprintf("%.2f",pars2), sprintf("%.2f",pars3), sprintf("%.2f",pars4))
    tags.env <- paste0("<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(", respfxn_num,"),", env_pars,", Indexed.Single[]")
  }
  
  tags = c(tags.vul,tags.env,tags.disp)
  
  
  #command files----
  cmd_j <- cmd_base
  cmd_j[startsWith(cmd_j, "<ECOSPACE_OUTPUT_DIR>")] <-
    sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", run_path)
  cmd_j <- c(cmd_j, tags)
  cmd_file <- file.path(run_path, "cmd.txt")
  writeLines(cmd_j, cmd_file)
  return(cmd_file)
  
}


safe_runEwE <- function(cmdfile, do.obj) {
  on.exit(gc(), add = TRUE)
  out <- tryCatch({
    # example: quiet system call
    # system2("EwE.exe", args = c(cmdfile), stdout = FALSE, stderr = FALSE, wait = TRUE)
    fn.runEwE(cmdfile = cmdfile, do.obj = do.obj)
  }, error = function(e) NA_real_)
  if (is.numeric(out) && length(out) >= 1 && is.finite(out[1])) out[1] else NA_real_
}#eof


# safe_runEwE <- function(cmdfile, do.obj, timeout_sec = 1) {
#   if (requireNamespace("R.utils", quietly = TRUE)) {
#     R.utils::withTimeout({
#       fn.runEwE(cmdfile = cmdfile, do.obj = do.obj)
#     }, timeout = timeout_sec, onTimeout = "silent")
#   } else {
#     # Fallback without timeout (consider installing R.utils)
#     fn.runEwE(cmdfile = cmdfile, do.obj = do.obj)
#   }
# }#eof


# safe_runEwE <- function(cmdfile, do.obj) {
#   on.exit(gc(), add = TRUE)   # ensure cleanup
#   out <- tryCatch(
#     fn.runEwE(cmdfile = cmdfile, do.obj = do.obj),
#     error = function(e) NA_real_)
#     if (is.numeric(out) && length(out) >= 1 && is.finite(out[1])){
#       out[1]
#     } else {
#       NA_real_
#     }
# }#eof


#run the population of models----
fn.runEwE.gapop <-  function(
    files.cmd, 
    obj.fxn=1, 
    delete.output=T,
    max_tries=5
){
  #source(file.setup)

  #clusterExport(cl,append(cl.export,c("file.console", "fn.runEwE", "safe_runEwE","fn.objfxn1", "fn.objfxn2","styear","enyear","group.names","df.names","fn.ecospace_predB_ts2array","fn.ecospace_predC_ts2array")))

  #pbar <- winProgressBar("Running Ecospace GApop",label=paste0("Simulation 0 of ",length(files.cmd)),max=100)
  #prog <- function(n) setWinProgressBar(pbar,(n/length(files.cmd)*100),label=paste("Simulation Run", n,"of",length(files.cmd),"Completed"))
  #opts <- list(progress=prog)
  
  #print(paste('Running',length(files.cmd),'Ecospace simulations'))
  t1 <- Sys.time()
  runLL.out <- foreach(i = seq_along(files.cmd), .errorhandling = 'pass') %dopar% { #, .options.snow = opts
    #fn.runEwE(cmdfile=files.cmd[i], do.obj=obj.fxn)
    safe_runEwE(cmdfile=files.cmd[i], do.obj = obj.fxn)
  }
  runLL <- unlist(runLL.out)
  #close(pbar)
  #print(paste('Run time',round(as.numeric(Sys.time()-t1),2)))
  
  ##missing runs----------------------------------------------------------------
  #tmp.errs = sample(1:340,10)
  #runLL[11:15] <- NA_real_ 
  erruns <- which(is.na(runLL))
  #filecheck <- sapply(dirname(files.cmd),FUN=function(x)length(list.files(x)))
  #erruns <- which(filecheck<=1)  
  tries <- 0
  while(length(erruns)>=1 && tries<max_tries){
    tries <- tries+1
    #print(paste0('Redo missing runs: n=',length(erruns)))
    
    #pbar <- winProgressBar("Running Ecospace GApop: Missing Runs",label=paste0("Simulation 0 of ",length(erruns)),max=100)
    #prog <- function(n) setWinProgressBar(pbar,(n/length(erruns)*100),label=paste("Simulation Run", n,"of", length(erruns),"Completed"))
    #opts <- list(progress=prog)
    
    runLL.err.out <- foreach(i=1:length(erruns),.errorhandling='pass') %dopar% {#,.options.snow=opts
      #fn.runEwE(cmdfile=files.cmd[erruns[i]], do.obj=obj.fxn)
      safe_runEwE(cmdfile=files.cmd[erruns[i]], do.obj=obj.fxn)
    }
    runLL.err <- unlist(runLL.err.out)
    #close(pbar)
    
    #for(k in 1:length(erruns)) runs[[erruns[k]]] <- unlist(runs.erruns[k])
    for (k in seq_along(erruns)) runLL[erruns[k]] <- runLL.err[k]
    
    #filecheck <- sapply(dirname(files.cmd),FUN=function(x)length(list.files(x)))
    #erruns <- which(filecheck<=1)  #if there are many missing runs, then need to do this in parallel
    #runLL[16] <- NA_real_
    erruns <- which(is.na(runLL))
  }
  
  
  # Optional: warn if some runs still missing after retries
  if (length(erruns) > 0) {
    warning(sprintf("Incomplete: %d runs still missing after %d retries at indices: %s",
                    length(erruns), tries, paste(erruns, collapse = ", ")))
  }
  
  # Clean up any generated output directories, if applicable
  if (delete.output && exists("run_dir")) {
    unlink(list.dirs(run_dir, full.names = TRUE, recursive = FALSE), recursive = TRUE)
  }
  
  # Return fitness vector (total log-likelihoods to minimize)
  fitness <- runLL
  # print(paste('Run time', round(as.numeric(Sys.time() - t1), 2)))
  return(fitness)
} # eof

# === Selection ===
select_parents <- function(gapop, fitness) {
  ranks <- rank(-fitness, ties.method='random')
  probs <- ranks / sum(ranks)
  #cbind(ranks,fitness,probs)
  selected <- gapop[sample(1:nrow(gapop), nrow(gapop), replace = TRUE, prob = probs), ]
  return(selected)
}

# === Crossover ===
crossover <- function(parents) {
  offspring <- parents
  for (i in seq(1, pop_size - 1, by = 2)) {
    if (runif(1) < 0.8) {
      point <- sample(1:(n_pars - 1), 1)
      temp <- offspring[i, (point + 1):n_pars]
      offspring[i, (point + 1):n_pars] <- offspring[i + 1, (point + 1):n_pars]
      offspring[i + 1, (point + 1):n_pars] <- temp
    }
  }
  return(offspring)
}

crossover_uniform <- function(parents, p_cross = 0.8, p_gene = 0.5) {
  offspring <- parents
  pop_size <- nrow(parents)
  n_pars <- ncol(parents)
  for (i in seq(1, pop_size - 1, by = 2)) {
    if (runif(1) < p_cross) {
      mask <- runif(n_pars) < p_gene  #mask is the index of parameters to be crossed
      # swap masked genes
      tmp <- offspring[i, mask]
      offspring[i, mask] <- offspring[i + 1, mask]
      offspring[i + 1, mask] <- tmp
    }
  }
  return(offspring)
} #eof

# === Mutation ===
mutate <- function(population, margin=0.2) {
  #low = L.bounds
  #upp = U.bounds
  #original scale
  low = exp(apply(population,2,min))
  upp = exp(apply(population,2,max))
  low = low-margin*abs(low)
  low[vul.par.idx] <- ifelse(low[vul.par.idx]<=1.0,1.01,low[vul.par.idx])
  upp = upp+margin*abs(upp)
  #back to log scale
  low = log(low)
  upp = log(upp)
  for (i in 1:nrow(population)) {
        mask <- runif(n_pars) < mutation_rate #mask are the parameters to mutate
        population[i, mask] <- runif(sum(mask), low[mask], upp[mask])
  }
  return(population)
} #eof


cluster_is_ok <- function() {
  ok <- tryCatch({
    res <- foreach(i = 1:2, .combine = c) %dopar% { Sys.getpid() }
    length(res) == 2
  }, error = function(e) FALSE)
  ok
}

ensure_cluster <- function() {
  if (!cluster_is_ok()) {
    # rebuild
    try(silent = TRUE, stopCluster(cl))
    rm(list = "cl", envir = .GlobalEnv)
    gc()
    closeAllConnections()
    workers <- max(1L, floor(detectCores() / 4))
    cl <<- parallel::makePSOCKcluster(workers, outfile = "cluster_workers.log")
    doParallel::registerDoParallel(cl)
    clusterExport(cl, c(
      # everything used from worker side:
      "safe_runEwE", "fn.runEwE", "fn.objfxn1", "fn.objfxn2",
      "styear", "enyear", "group.names", "df.names", "obs.ts",
      "cmd_base", "myconfig", "run_dir"
    ))
  }
}


#GA function----
fn.GA <- function(myconfig){
  #unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)
  pop_size <<- myconfig$popSize
  n_generations <<- myconfig$n_gen
  mutation_rate <<- myconfig$pmutation
  elitism <<- myconfig$elitism
  do.penalty <<- myconfig$do.penalty
  pen.wt.mult <<- myconfig$pen.wt.mult
  
  #create results file
  file.ga.results <- file.path(dir.main, paste0('ga_results_',timestamp,'.csv'))
  myconfig.string <- paste(paste(names(myconfig),unlist(myconfig),sep="="),collapse="; ")
  write.table('Ecospace GA calibration', file.ga.results,row.names=F, col.names=F,append=F)
  write.table(Sys.time(), file.ga.results,row.names=F, col.names=F,append=T)
  write.table(cmd_base[which(startsWith(cmd_base,"<EWE_MODEL_FILE>"))],file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(cmd_base[which(startsWith(cmd_base,"<ECOSIM_SCENARIO_INDEX>"))],file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(cmd_base[which(startsWith(cmd_base,"<ECOSPACE_SCENARIO_INDEX>"))],file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(myconfig.string,file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(t(c("gen_num","min_fitness","mean_fitness",names(log_par_vec))), file.ga.results, sep=",", row.names = F, col.names = F, append=T)
  
  #base run
  message('Running the base model')
  files.cmd <- fn.parvec2cmd(log_par_vec=log_par_vec,g=999,idx=0)
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1)  
  
  #pipe output
  g=-1
  best_fit = min(fitness, na.rm=T)
  mean_fit = mean(fitness, na.rm=T)
  sd_fit = sd(fitness, na.rm=T)
  max_fit = max(fitness, na.rm=T)
  median_fit = median(fitness, na.rm=T)
  best_pars = round(exp(log_par_vec),4)
  write.table(t(c(g,best_fit,mean_fit,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
  cat(sprintf("%-10s %-15s %-15s %-15s %-15s %-15s | %-10s \n", "Gen", "min_LL", "max_LL", "mean_LL","median_LL","sd_LL","time end"))
  cat(sprintf("%-10d %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f | %s\n", g, best_fit, max_fit, mean_fit, median_fit, sd_fit, Sys.time()))

  #initial population.................................
  #message('Running the initial population')
  gapop <- fn.GApop(usedist='ln')
  gapop[1,] <- log_par_vec  #include base run in initial population
  #files.cmd <- apply(gapop,1,function(x) fn.parvec2cmd(log_par_vec=x, g=0, idx=0)) 
  files.cmd <- lapply(1:nrow(gapop),function(i) fn.parvec2cmd(log_par_vec=gapop[i,], g=0, idx=i)) 
  files.cmd <- unlist(files.cmd, use.names=F)
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1)
  
  #calculate penalty for parameter bounds violations
  if(do.penalty){
  pen.wt = abs(diff(range(fitness)))*pen.wt.mult
  upp.pen = U.bounds*1.5
  low.pen = L.bounds*0.5
  penalties = rep(0,nrow(gapop))
  for(i in 1:nrow(gapop)){
    #i=119
    penalty=rep(0,ncol(gapop))
    penalties[i] <- sum((pmax(0,gapop[i,]-upp.pen)^2 + pmax(0,low.pen-gapop[i,])^2)*pen.wt)
  }
  fitness <- fitness+penalties
  }
  #pipe output
  g=0
  best_fit = min(fitness, na.rm=T)
  mean_fit = mean(fitness, na.rm=T)
  sd_fit = sd(fitness, na.rm=T)
  max_fit = max(fitness, na.rm=T)
  median_fit = median(fitness, na.rm=T)
  best_pars = round(exp(gapop[which.min(fitness),]),4)
  write.table(t(c(g,best_fit,mean_fit,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
  cat(sprintf("%-10d %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f | %s\n", g, best_fit, max_fit, mean_fit, median_fit, sd_fit, Sys.time()))
  #message(sprintf("Generation 0: Lowest LL score  = %.4f\n", min(fitness)))
  #message(sprintf("Generation 0: Avg pop LL score  = %.4f\n", mean(fitness)))
  
  for (gen in 1:n_generations) {
    #gen=1
    if(gen==1) unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)
    
    #reset workers after every 5 generations - this didn't help but we'll keep it jic
    if(gen %in% seq(6,n_generations,5)){
      try(silent = TRUE, stopCluster(cl))
      rm(cl); gc()
      closeAllConnections()
      workers <- floor(detectCores() / 4)
      cl <- parallel::makePSOCKcluster(workers, outfile = "cluster_workers.log")
      doParallel::registerDoParallel(cl)
      clusterExport(cl, c("file.console","cmd_base","myconfig","run_dir",
                          "safe_runEwE","fn.runEwE","fn.objfxn1","fn.objfxn2",
                          "fn.ecospace_predB_ts2array","fn.ecospace_predC_ts2array",
                          "styear","enyear","group.names","df.names","obs.ts"))
    }
    
    # Elitism - keep the top n runs
    elite_idx <- order(fitness, decreasing=F)[1:elitism]
    elite <- gapop[elite_idx, ]
    
    # Selection, Crossover, Mutation
    parents <- select_parents(gapop, fitness) #resamples the population, with replacement, with rank-based probabilities in the sample draws
    #offspring <- crossover(parents) #offspring are when two parents crossover a part of their parameter vector
    offspring <- crossover_uniform(parents, p_cross=0.8, p_gene=0.5)
    offspring <- mutate(offspring, margin=0.2) #randomly draw new parameter values to mutate the individual
    
    # Evaluate new population
    ensure_cluster()  # optional but recommended
    files.cmd <- lapply(1:nrow(offspring),function(i) fn.parvec2cmd(log_par_vec=offspring[i,], g=gen, idx=i)) 
    files.cmd <- unlist(files.cmd, use.names=F)
    new_fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1, delete.output=T)
    
    #calculate penalty for parameter bounds violations
    # pen.wt = abs(diff(range(fitness)))*0.1
    # upp.pen = U.bounds*1.5
    # low.pen = L.bounds*0.5
    if(do.penalty){
    penalties = rep(0,nrow(offspring))
    for(i in 1:nrow(offspring)){
      #i=119
      penalties[i] <- sum((pmax(0,offspring[i,]-upp.pen)^2 + pmax(0,low.pen-offspring[i,])^2)*pen.wt)
    }
    new_fitness <- new_fitness+penalties
    }
    
    # Combine elite + offspring: keep elites and the best (lowest) offspring
    worst_offspring_idx <- order(new_fitness, decreasing = TRUE)[1:elitism]  # drop worst 'elitism'
    gapop   <- rbind(elite, offspring[-worst_offspring_idx, ])
    fitness <- c(fitness[elite_idx], new_fitness[-worst_offspring_idx])
    
    # Combine elite + offspring, drop worst offspring and replace with elites
    # offspring.rank = rank(-new_fitness, ties.method='random')
    # drop.idx = which(offspring.rank<=elitism)
    # gapop <- rbind(elite, offspring[-drop.idx,])
    # fitness <- c(fitness[elite_idx], new_fitness[-drop.idx])
    
    #pipe output
    g=gen
    best_fit = min(fitness, na.rm=T)
    mean_fit = mean(fitness, na.rm=T)
    sd_fit = sd(fitness, na.rm=T)
    max_fit = max(fitness, na.rm=T)
    median_fit = median(fitness, na.rm=T)
    best_pars = round(exp(gapop[which.min(fitness),]),4)
    write.table(t(c(g,best_fit,mean_fit,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
    cat(sprintf("%-10d %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f | %s\n", g, best_fit, max_fit, mean_fit, median_fit, sd_fit, Sys.time()))
    #message(sprintf("Generation %d: Avg pop LL score = %.4f\n", gen, mean(fitness)))
    #unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)
    rm(parents,offspring,files.cmd);gc()
  }
}#eof



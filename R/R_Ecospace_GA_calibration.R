
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
  disp.cv=0.1,
  do.med = TRUE,
  med_pars = med_pars,
  med.xbase.cv = .1){
  
  # do.vuls = TRUE
  # vul_pars = predprey_pars
  # vul.min = 1.01
  # vul.max = 1e6
  # vul.cv = 0.6
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
  # do.med = FALSE
  # med_pars = NULL
  # med.xbase.cv = .1
  
  #validate setup----
  if(do.vuls){
    pdcol = unique(vul_pars$pred[is.na(vul_pars$prey)])
    pdpy = c(0,unique(vul_pars$pred[!is.na(vul_pars$prey)]))
    if(length(which(pdcol %in% pdpy))>0){
      stop("Trying to estimate ki and kij for at least on predator. Check vul_pars input.")
    }
  }
  
  #build parameter vector----
  n_vuls <- n_env <- n_disp <- n_med <- 0
  if(do.vuls) n_vuls <- nrow(vul_pars)
  if(do.env) n_env <- nrow(env_pars)  #need to melt the env resp fxn dataframe
  if(do.disp) n_disp <- nrow(disp_pars)
  if(do.med) n_med <- nrow(med_pars)
  #medations functions
  
  n_pars <- n_vuls + n_env + n_disp + n_med
  if (n_pars == 0) stop("n_pars==0: No predator-prey pairs or environmental responses parameters to estimate.  Check inputs.")
  
  
  #because med_pars and likely env_pars are more constrained and could be negative, we do not log them
  #why log any of them?
  #log_vuln_vec <- env_vec <- log_disp_vec <- med_vec <- numeric()
  vuln_vec <- env_vec <- disp_vec <- med_vec <- numeric()
  par.idx <- par.labels <- par.groups <- character()
  #log.par.idx <- logical()
  
  if(do.vuls) {
    #log_vuln_vec <- log(vul_pars$base.val)
    vuln_vec <- vul_pars$base.val
    par.idx <- c(par.idx,rep('vul',n_vuls))
    par.labels <- c(par.labels,paste('vul',vul_pars$pred,vul_pars$prey,sep="_"))
    par.groups <- c(par.groups,paste('vul',vul_pars$pred,vul_pars$prey,sep="_"))
    #log.par.idx <- c(log.par.idx,rep(TRUE,length(log_vuln_vec)))
  }
  if(do.env) {
    env_vec = rep(0,n_env)
    par.idx = c(par.idx,rep('env',n_env))
    par.labels = c(par.labels,paste('env',env_pars$Function.number,sep="_"))
    par.groups = c(par.groups,paste('env',env_pars$Function.number,sep="_"))
    #log.par.idx <- c(log.par.idx, rep(FALSE,length(env_vec)))
  }
  if(do.disp){
    #log_disp_vec = log(disp_pars$base.val)
    disp_vec = disp_pars$base.val
    par.idx = c(par.idx,rep('disp',n_disp))
    par.labels = c(par.labels,paste('disp',disp_pars$pred,sep="_"))
    par.groups = c(par.groups,paste('disp',disp_pars$pred,sep="_"))
    #log.par.idx <- c(log.par.idx,rep(TRUE,length(log_disp_vec)))
  }
  
  med.cv = numeric()
  if(do.med){
    for(i in 1:n_med){
      #i=2
      #HYPERBOLIC SHAPE
      if(med_pars$shape_type[i]==3){
        #par_vec.i <- c(log(med_pars$x_base[i]), med_pars$par1[i], med_pars$par2[i], med_pars$par3[i])
        par_vec.i <- c(med_pars$x_base[i], med_pars$par1[i], med_pars$par2[i], med_pars$par3[i])
        par.labels.i <- paste(paste0('med',rep(med_pars$med_index[i],length(par_vec.i))), paste0('type',rep(med_pars$shape_type[i],length(par_vec.i))), c('xbase','Yzero','Yend','Ybase'),sep="_")
        par.groups.i <- paste0('med',rep(med_pars$med_index[i],length(par_vec.i)))
      }
      #LINEAR SHAPE
      if(med_pars$shape_type[i]==1){
        slope.i <- med_pars$par2[i] - med_pars$par1[i]
        #par_vec.i <- c(log(med_pars$x_base[i]), slope.i)
        par_vec.i <- c(med_pars$x_base[i], slope.i)
        par.labels.i <- paste(paste0('med',rep(med_pars$med_index[i],length(par_vec.i))), paste0('type',rep(med_pars$shape_type[i],length(par_vec.i))), c('xbase','slope'),sep="_")
        par.groups.i <- paste0('med',rep(med_pars$med_index[i],length(par_vec.i)))
      }
      med_vec <- c(med_vec,par_vec.i)
      par.labels <- c(par.labels,par.labels.i)
      par.groups <- c(par.groups,par.groups.i)
      #log.par.idx <- c(log.par.idx,c(TRUE,rep(FALSE,length(par_vec.i)-1)))
    }
    par.idx <- c(par.idx,rep('med',length(med_vec)))
  }
  
  #n_pars <<- length(log_vuln_vec) + length(env_vec) + length(log_disp_vec) + length(med_vec)
  n_pars <<- length(vuln_vec) + length(env_vec) + length(disp_vec) + length(med_vec)
  if (n_pars == 0) stop("n_pars==0: No parameters to estimate.  Check inputs.")
  
  #est_par_vec <- c(log_vuln_vec, env_vec, log_disp_vec, med_vec)
  est_par_vec <- c(vuln_vec, env_vec, disp_vec, med_vec)
  names(est_par_vec) <- par.labels
  
  # index parameter types
  vul.par.idx <<- which(par.idx=='vul')
  env.par.idx <<- which(par.idx=='env')
  disp.par.idx <<- which(par.idx=='disp')
  med.par.idx <<- which(par.idx=='med')
  #log.par.idx <<- log.par.idx
  
  # super assignment so these are available
  vul_pars <<- vul_pars
  disp_pars <<- disp_pars
  med_pars <<- med_pars
  
  # set parameter bounds--------------------------------------------------------
  ##vulnerabilities - do these on the log scale
  lower.vuls <- upper.vuls <- numeric()
  lower.vuls <- exp(log(vul_pars$base.val) + exp(vul.cv)^2 * qnorm(0.05))
  lower.vuls <- ifelse(lower.vuls<=1.0,1.01,lower.vuls)
  upper.vuls <- exp(log(vul_pars$base.val) + exp(vul.cv)^2 * qnorm(0.95))
  
  ##environmental shapes - TBD
  lower.env <- upper.env <- numeric()
  lower.env <- rep(env.min,n_env)
  upper.env <- rep(env.max,n_env)
  
  ##dispersal
  lower.disp <- upper.disp <- numeric()
  lower.disp <- disp_vec-2*disp.cv*disp_vec
  upper.disp <- disp_vec+2*disp.cv*disp_vec
  
  ##mediation
  lower.med <- upper.med <- numeric()  
  if(n_med>0){
  med.cv = numeric()  
  for(i in 1:length(med_vec)){
    #i=1
    label.i <- par.labels[med.par.idx[i]]
    shape.i <- as.numeric(gsub('med','',unlist(strsplit(label.i,"_"))[1]))
    type.i <- as.numeric(gsub('type','',unlist(strsplit(label.i,"_"))[2]))
    par.i <- unlist(strsplit(label.i,"_"))[3]
    val.i <- med_vec[i]
    if(par.i=="xbase"){
      cv.i <- med.xbase.cv
      lower.i <- floor(val.i-2*val.i*med.xbase.cv)
      upper.i <- ceiling(val.i+2*val.i*med.xbase.cv)
      #cv.i = med.xbase.cv
    } else{
      cv.i <- med_pars$cv[med_pars$med_index==shape.i]
        if(val.i>0) lower.i <- val.i-2*cv.i*val.i
        if(val.i<0) lower.i <- val.i+2*cv.i*val.i
        if(val.i>0) upper.i <- val.i+2*cv.i*val.i
        if(val.i<0) upper.i <- val.i-2*cv.i*val.i
      }
    lower.med <- c(lower.med,lower.i)
    upper.med <- c(upper.med,upper.i)
    med.cv <- c(med.cv,cv.i)
  }
  }
  
  #super assignment
  est_par_vec <<- round(est_par_vec,2)
  L.bounds <<- round(c(lower.vuls,lower.env,lower.disp,lower.med),2)
  U.bounds <<- round(c(upper.vuls, upper.env, upper.disp,upper.med),2)
  n_vuls <<- n_vuls
  n_env <<- n_env
  n_disp <<- n_disp
  n_med <<- length(med_vec)
  par.groups <<- par.groups
  par.labels <<- par.labels
  
  #create vector of cv for rnorm() draws
  par_cv_vec <- est_par_vec
  par_cv_vec[vul.par.idx] <- vul.cv
  par_cv_vec[env.par.idx] <- env.cv
  par_cv_vec[disp.par.idx] <- disp.cv
  par_cv_vec[med.par.idx] <- med.cv
  par_cv_vec <<- par_cv_vec
  
} #end function

#make GA populations----
fn.GApop = function(run_config=myconfig, usedist='unif'){
  #run_config <- myconfig
  if(usedist=='unif'){
    mat <- matrix(runif(run_config$popSize*n_pars, L.bounds, U.bounds), nrow=run_config$popSize, byrow=T)
  }
  if(usedist=='ln'){
    mat <- matrix(NA,nrow=run_config$popSize, ncol=n_pars)
    colnames(mat) <- par.labels
    par_sd <- rep(0,n_pars)
    par_sd[log.par.idx] <- sqrt(log(par_cv_vec[log.par.idx]^2+1))
    par_sd[!log.par.idx] <- sqrt((par_cv_vec[!log.par.idx]+1))
    tmp.log <- matrix(rlnorm(run_config$popSize*sum(log.par.idx), meanlog=est_par_vec[log.par.idx], sdlog=par_sd[log.par.idx]),nrow=run_config$popSize, byrow=T)
    tmp.norm <- matrix(rnorm(run_config$popSize*sum(!log.par.idx), mean=est_par_vec[!log.par.idx], sd=par_sd[!log.par.idx]),nrow=run_config$popSize, byrow=T)
    mat[,log.par.idx] <- tmp.log
    mat[,!log.par.idx] <- tmp.norm
    
    idx <- mat[, vul.par.idx] < 1.0
    mat[,vul.par.idx][idx] <- 1.01
    integer.idx <- grep("xbase",par.labels)
    mat[,integer.idx] <- round(mat[,integer.idx])
    mat[,log.par.idx] <- log(mat[,log.par.idx])
    
    # graphics.off();rm(.SavedPlots);windows(record=T)
    # par(mfrow=c(3,3))
    # for(i in 1:ncol(mat)){
    #   hist(mat[,i],main=par.labels[i])
    # }
  }
  return(mat)
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
fn.parvec2cmd <- function(par_vec=est_par_vec, g=gen, idx=i, out_dir=run_dir){
  
  # par_vec = est_par_vec
  # g = 999
  # idx = 0
  # out_dir=run_dir
  
  #output directory----
  run_path <- safe_make_run_dir_tempfile(base_dir=out_dir,
                                         par_vec=par_vec,
                                         gen = g,
                                         idx = idx,
                                         prefix = "run",
                                         random_len = 8)
  
  #parameter tags----
  tags.vul <- tags.env <- tags.disp <- tags.med <- character()
  if(length(vul.par.idx)>0){
    vuln_vec <- par_vec[vul.par.idx]
    tags.vul <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", vul_pars$pred, " ", ifelse(is.na(vul_pars$prey),"",vul_pars$prey),
                       "), ", sprintf("%.5f", vuln_vec), ", Indexed.Single")
  }

  if(length(disp.par.idx)>0){
    disp_vec <- par_vec[disp.par.idx]
    tags.disp <- paste0("<ECOSPACE_DISPERSAL_RATE_INDEXED>(", disp_pars$pred,"), ",disp_vec, ", Indexed.Single")
  }
  
  if(length(env.par.idx)>0){
    env_pars = paste(respfxn_type,sprintf("%.2f",pars1), sprintf("%.2f",pars2), sprintf("%.2f",pars3), sprintf("%.2f",pars4))
    tags.env <- paste0("<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(", respfxn_num,"),", env_pars,", Indexed.Single[]")
  }
  
  if(length(med.par.idx)>0){
    med_groups <- unique(par.groups[med.par.idx])
    tags.med <- character()
    for(i in 1:length(med_groups)){
      #i=1
      par_set.i <- par_vec[which(par.groups==med_groups[i])]
      par_labels.i <- par.labels[which(par.groups==med_groups[i])]
      type.i <- as.numeric(gsub('type','',strsplit(par_labels.i[1],"_")[[1]][2]))
      shape.i <- as.numeric(gsub('med','',strsplit(par_labels.i[1],"_")[[1]][1]))
      #LINEAR SHAPE
      if(type.i==1){
        par.string.i <- paste(par_set.i[1], type.i, 1-0.5*par_set.i[2], 1+0.5*par_set.i[2], sep=" ")
        tag.i <- paste0("<MEDIATION_FUNCTION_INDEXED>(",shape.i,"), ",par.string.i,", Indexed.Single[]")
      }
      #HYPERBOLIC SHAPE
      if(type.i==3){
        par.string.i <- paste(par_set.i[1], type.i, par_set.i[2], par_set.i[3], par_set.i[4], 1, sep=" ")
        tag.i <- paste0("<MEDIATION_FUNCTION_INDEXED>(",shape.i,"), ",par.string.i,", Indexed.Single[]")
      }
      
      tags.med <- c(tags.med,tag.i)
    }
  }
  
  tags = c(tags.vul,tags.env,tags.disp, tags.med)
  
  
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
    #i=1
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
  gapop.dist <<- myconfig$gapop.dist
  
  #create results file
  file.ga.results <- file.path(dir.main, paste0('ga_results_',timestamp,'.csv'))
  myconfig.string <- paste(paste(names(myconfig),unlist(myconfig),sep="="),collapse="; ")
  write.table('Ecospace GA calibration', file.ga.results,row.names=F, col.names=F,append=F)
  write.table(Sys.time(), file.ga.results,row.names=F, col.names=F,append=T)
  write.table(cmd_base[which(startsWith(cmd_base,"<EWE_MODEL_FILE>"))],file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(cmd_base[which(startsWith(cmd_base,"<ECOSIM_SCENARIO_INDEX>"))],file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(cmd_base[which(startsWith(cmd_base,"<ECOSPACE_SCENARIO_INDEX>"))],file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(myconfig.string,file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(t(c("gen_num","min_fitness","mean_fitness",names(est_par_vec))), file.ga.results, sep=",", row.names = F, col.names = F, append=T)
  
  #base run
  files.cmd <- fn.parvec2cmd(par_vec=est_par_vec,g=999,idx=0)
  message('Running the base model')
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1)  
  
  #pipe output
  g=-1
  best_fit = min(fitness, na.rm=T)
  mean_fit = mean(fitness, na.rm=T)
  sd_fit = sd(fitness, na.rm=T)
  max_fit = max(fitness, na.rm=T)
  median_fit = median(fitness, na.rm=T)
  best_pars = est_par_vec 
  write.table(t(c(g,best_fit,mean_fit,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
  cat(sprintf("%-10s %-15s %-15s %-15s %-15s %-15s | %-10s \n", "Gen", "min_LL", "max_LL", "mean_LL","median_LL","sd_LL","time end"))
  cat(sprintf("%-10d %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f | %s\n", g, best_fit, max_fit, mean_fit, median_fit, sd_fit, Sys.time()))

  #initial population.................................
  #PICKUP HERE - NEED TO CHECK gn.GApop to work with mediation (1/26/2020)
  #message('Running the initial population')
  gapop <- fn.GApop(usedist=gapop.dist)
  gapop[1,] <- est_par_vec  #include base run in initial population
  #files.cmd <- apply(gapop,1,function(x) fn.parvec2cmd(log_par_vec=x, g=0, idx=0)) 
  files.cmd <- lapply(1:nrow(gapop),function(i) fn.parvec2cmd(par_vec=gapop[i,], g=0, idx=i)) 
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



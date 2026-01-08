
#prepare parameter vector----
fn.makeparvec <- function(
  do.vuls = TRUE, 
  vul_pars = predprey_pars,
  vul.min = 1.01,
  vul.max = 1e6,
  vul.cv = 0.3,
  do.env = FALSE, 
  env_pars = NULL,
  env_pars.min=0.25,
  env_pars.max=2,
  do.disp = TRUE,
  disp_pars = disp_pars,
  disp.min=3,
  disp.max=3000,
  disp.cv=0.1){
  
  # do.vuls = TRUE
  # vul_pars = predprey_pars
  # vul.min = 1.01
  # vul.max = 1e6
  # vul.cv = 0.3
  # do.env = FALSE
  # env_pars = NULL
  # env_pars.min=0.25
  # env_pars.max=2
  # do.disp = TRUE
  # disp_pars = disp_pars
  # disp.min=3
  # disp.max=3000
  # disp.cv=0.1
  
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
  lower.vuls <- vul_pars$base.val-1.96*vul_pars$base.val*vul.cv
  lower.vuls <- ifelse(lower.vuls<=1.0,1.01,lower.vuls)
  lower.vuls <- log(lower.vuls)
  upper.vuls <- log(vul_pars$base.val+1.96*vul_pars$base.val*vul.cv)

  lower.env <- rep(-1,n_env)
  upper.env <- rep(1,n_env)
  #lower.disp <- rep(log(disp.min),n_disp)
  #upper.disp <- rep(log(disp.max),n_disp)
  lower.disp <- log(disp_pars$base.val-1.96*disp_pars$base.val*disp.cv)
  upper.disp <- log(disp_pars$base.val+1.96*disp_pars$base.val*disp.cv)
  
  #super assignment
  log_par_vec <<- log_par_vec
  L.bounds <<- c(lower.vuls,lower.env,lower.disp)
  U.bounds <<- c(upper.vuls, upper.env, upper.disp)
  n_vuls <<- n_vuls
  n_env <<- n_env
  n_disp <<- n_disp
} #end function

#make GA populations----
fn.GApop = function(object){
  run_config <- myconfig
  matrix(runif(run_config$popSize*n_pars, L.bounds, U.bounds), nrow=run_config$popSize, byrow=T)
  #matrix(rep(log_par_vec, run_config$popSize), nrow=run_config$popSize, byrow=T)
}


# Alternative minimalistic approach:
safe_make_run_dir_tempfile <- function(run_dir, log_par_vec) {
  h <- digest::digest(log_par_vec, algo = "md5", serialize = TRUE)
  run_path <- tempfile(pattern = paste0("run_", h, "_"), tmpdir = run_dir)
  dir.create(run_path, showWarnings = TRUE, recursive = TRUE)
  run_path
}


#write command files----
fn.parvec2cmd <- function(log_par_vec){
  #par_vec <- log_par_vec
  #log_par_vec=gapop.init[1,]
  
  #output directory----
  #run_id <- paste0("run_", digest(log_par_vec, algo=c("md5")))
  #run_path <- file.path(run_dir, run_id)
  #dir.create(run_path, showWarnings = TRUE)
  run_path <- safe_make_run_dir_tempfile(run_dir, log_par_vec)
  
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

#run the population of models----
fn.runEwE.gapop <-  function(
    files.cmd, 
    obj.fxn=1, 
    cl.export = list('files.cmd','obs.ts')
){
  #source(file.setup)

  clusterExport(cl,append(cl.export,list("file.console", "fn.runEwE", "fn.objfxn1", "fn.objfxn2","styear","enyear","group.names","df.names","fn.ecospace_predB_ts2array","fn.ecospace_predC_ts2array")))

  #runlist=runlist[1:20,]
  pbar <- winProgressBar("Running Ecospace GApop",label=paste0("Simulation 0 of ",length(files.cmd)),max=100)
  prog <- function(n) setWinProgressBar(pbar,(n/length(files.cmd)*100),label=paste("Simulation Run", n,"of", 
                                                                                   length(files.cmd),"Completed"))
  opts <- list(progress=prog)
  
  #print(paste('Running',length(files.cmd),'Ecospace simulations'))
  t1 <- Sys.time()
  runs <- foreach(i = 1:length(files.cmd), .errorhandling = 'pass', .options.snow = opts) %dopar% {
    fn.runEwE(cmdfile=files.cmd[i], do.obj = obj.fxn)
  }
  close(pbar)
  #print(paste('Run time',round(as.numeric(Sys.time()-t1),2)))
  
  ##missing runs----
  filecheck <- sapply(dirname(files.cmd),FUN=function(x)length(list.files(x)))
  erruns <- which(filecheck<=1)  
  
  while(length(erruns)>=1){
    #print(paste0('Redo missing runs: n=',length(erruns)))
    
    pbar <- winProgressBar("Running Ecospace GApop: Missing Runs",label=paste0("Simulation 0 of ",length(erruns)),max=100)
    prog <- function(n) setWinProgressBar(pbar,(n/length(erruns)*100),label=paste("Simulation Run", n,"of", length(erruns),"Completed"))
    opts <- list(progress=prog)
    
    runs.erruns <- foreach(i=1:length(erruns),.errorhandling='pass',.options.snow=opts) %dopar% {
      fn.runEwE(cmdfile=files.cmd[erruns[i]], do.obj=obj.fxn)
    }
    close(pbar)
    
    for(k in 1:length(erruns)) runs[[erruns[k]]] <- unlist(runs.erruns[k])
    
    filecheck <- sapply(dirname(files.cmd),FUN=function(x)length(list.files(x)))
    erruns <- which(filecheck<=1)  #if there are many missing runs, then need to do this in parallel
  }
  
  fitness <- sapply(runs, function(x) x[1])   #cbind(runlist, do.call(rbind, runs))
  #print('All runs completed')
  #stopCluster(cl);
  #closeAllConnections()
  unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)
  #fitness <- -1*fitness
  return(fitness) #fitness is the total loglikelihood that I want to MINIMIZE
} #eof

# === Selection ===
select_parents <- function(gapop, fitness) {
  ranks <- rank(-fitness, ties.method='random')
  probs <- ranks / sum(ranks)
  #cbind(ranks,fitness,probs)
  selected <- gapop[sample(1:nrow(gapop), pop_size, replace = TRUE, prob = probs), ]
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

# === Mutation ===
mutate <- function(population) {
  #population=gapop
  #low = apply(population,2,min)
  #upp = apply(population,2,max)
  low = L.bounds
  upp = U.bounds
  for (i in 1:nrow(population)) {
    for (j in 1:n_pars) {
      if (runif(1) < mutation_rate) {
        population[i, j] <- runif(1, low[j], upp[j])
      }
    }
  }
  return(population)
}


#GA function----
fn.GA <- function(myconfig){
  pop_size <<- myconfig$popSize
  n_generations <<- myconfig$n_gen
  mutation_rate <<- myconfig$pmutation
  elitism <<- myconfig$elitism
  
  #store results
  file.ga.results <- file.path(run_dir, 'ga_results.csv')
  write.table(t(c("gen_num","min_fitness","mean_fitness",names(log_par_vec))), file.ga.results, sep=",", row.names = F, col.names = F)
  
  #base run
  message('Running the base model')
  files.cmd <- fn.parvec2cmd(log_par_vec=log_par_vec)
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1, cl.export=list('obs.ts'))  
  
  #pipe output
  g=-1
  best_fit = min(fitness)
  mean_fit = mean(fitness)
  best_pars = round(exp(log_par_vec),4)
  write.table(t(c(g,best_fit,mean_fit,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
  message(sprintf("Base LL score  = %.4f\n", min(fitness)))
  
  #initial population
  gapop <- fn.GApop()
  gapop[1,] <- log_par_vec  #include base run in initial population
  files.cmd <- apply(gapop,1,function(x) fn.parvec2cmd(log_par_vec=x)) 
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1, cl.export = list("obs.ts"))
  
  #pipe output
  g=0
  best_fit = min(fitness)
  mean_fit = mean(fitness)
  best_pars = round(exp(gapop[which.min(fitness),]),4)
  write.table(t(c(g,best_fit,mean_fit,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
  cat(sprintf("%-10s %-15s %-15s\n", "Generation", "Min_LL", "Avg_LL"))
  cat(sprintf("%-10d %-15.6f %-15.6f\n", 0, best_fit, mean_fit))
  #message(sprintf("Generation 0: Lowest LL score  = %.4f\n", min(fitness)))
  #message(sprintf("Generation 0: Avg pop LL score  = %.4f\n", mean(fitness)))
  
  for (gen in 1:n_generations) {
    #gen=1
    # Elitism - keep the top n runs
    elite_idx <- order(fitness, decreasing=F)[1:elitism]
    elite <- gapop[elite_idx, ]
    
    # Selection, Crossover, Mutation
    parents <- select_parents(gapop, fitness) #resamples the population, with replacement, with rank-based probabilities in the sample draws
    offspring <- crossover(parents) #offspring are when two parents crossover a part of their parameter vector
    offspring <- mutate(offspring) #randomly draw new parameter values to mutate the individual
    
    # Evaluate new population
    files.cmd <- apply(offspring,1,function(x) fn.parvec2cmd(x)) 
    new_fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1, cl.export = list("obs.ts"))
    
    # Combine elite + offspring, drop worst offspring and replace with elites
    offspring.rank = rank(-new_fitness, ties.method='random')
    drop.idx = which(offspring.rank<=elitism)
    gapop <- rbind(elite, offspring[-drop.idx,])
    fitness <- c(fitness[elite_idx], new_fitness[-drop.idx])
    
    #pipe output
    g=gen
    best_fit = min(fitness)
    mean_fit = mean(fitness)
    best_pars = round(exp(gapop[which.min(fitness),]),4)
    write.table(t(c(g,best_fit,mean_fit,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
    cat(sprintf("%-10d %-15.6f %-15.6f\n", gen, best_fit, mean_fit))
    #message(sprintf("Generation %d: Lowest LL score = %.4f\n", gen, min(fitness)))
    #message(sprintf("Generation %d: Avg pop LL score = %.4f\n", gen, mean(fitness)))
  }
}

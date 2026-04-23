#GA function----
#' @title Run genetic algorithm calibration for Ecospace
#' @description
#' Sets GA controls from `myconfig`, writes a results CSV, runs a base model,
#' evaluates an initial population, then iterates generations with elitism,
#' selection, crossover, mutation, and optional parameter-bound penalties.
#' Progress and summary stats are printed each generation.
#'
#' @param myconfig List of GA settings, expected fields:
#'   - popSize: population size
#'   - n_gen: number of generations
#'   - pmutation: per-gene mutation probability
#'   - elitism: number of elites kept each generation
#'   - do.penalty: logical, apply penalties for bound violations
#'   - pen.wt.mult: scalar multiplier for penalty weight
#'   - gapop.dist: initial population distribution ("unif" or "ln")
#'
#' @details
#' Uses helper functions: fn.GApop, fn.parvec2cmd, fn.runEwE.gapop,
#' select_parents, crossover_uniform, mutate, and optionally ensure_cluster.
#' Writes a CSV to dir.main using a timestamped filename, and cleans output
#' directories as needed. Several values are read from the environment, such as
#' est_par_vec, dir.main, timestamp, cmd_base, run_dir, U.bounds, L.bounds.
#' Global variables pop_size, n_generations, mutation_rate, elitism, do.penalty,
#' pen.wt.mult, and gapop.dist are set for use by downstream helpers.
#'
#' @return Invisibly returns NULL. Results (fitness summaries and best parameters)
#' are appended to the CSV and printed to the console.

fn.GA <- function(myconfig){
  #unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)
  pop_size <<- myconfig$popSize
  n_generations <<- myconfig$n_gen
  mutation_rate <<- myconfig$pmutation
  elitism <<- myconfig$elitism
  do.penalty <<- myconfig$do.penalty
  pen.wt.mult <<- myconfig$pen.wt.mult
  gapop.pardist <<- myconfig$gapop.pardist
  gapop.vuldist <<- myconfig$gapop.vuldist
  mutate.margin <<- myconfig$mutate.margin
  dir.results <<- myconfig$dir.ga_results
  #small test change here
  
  #create results file
  file.ga.results <- file.path(dir.results.ga, paste0('ga_results_',timestamp,'.csv'))
  myconfig.string <- paste(paste(names(myconfig),unlist(myconfig),sep="="),collapse="; ")
  write.table('Ecospace GA calibration', file.ga.results,row.names=F, col.names=F,append=F)
  write.table(Sys.time(), file.ga.results,row.names=F, col.names=F,append=T)
  write.table(cmd_base[which(startsWith(cmd_base,"<EWE_MODEL_FILE>"))],file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(cmd_base[which(startsWith(cmd_base,"<ECOSIM_SCENARIO_INDEX>"))],file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(cmd_base[which(startsWith(cmd_base,"<ECOSPACE_SCENARIO_INDEX>"))],file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(myconfig.string,file.ga.results,file.ga.results,row.names=F, col.names=F,append=T, quote=F)
  write.table(t(c("gen","min_LL","max_LL","mean_LL","median_LL","sd_LL","NA_count",names(est_par_vec))), file.ga.results, sep=",", row.names = F, col.names = F, append=T)
  
  #base run
  files.cmd <- fn.parvec2cmd(par_vec=est_par_vec,g=999,idx=0)
  message('Running the base model')
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1,delete.output=T)  

  #pipe output
  g=-1
  best_fit = min(fitness[1], na.rm=T)
  mean_fit = mean(fitness[1], na.rm=T)
  sd_fit = sd(fitness[1], na.rm=T)
  max_fit = max(fitness[1], na.rm=T)
  median_fit = median(fitness[1], na.rm=T)
  na_count = sum(is.na(fitness))
  best_pars = est_par_vec
  write.table(t(c(g, best_fit, max_fit, mean_fit, median_fit, sd_fit, na_count,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
  cat(sprintf("%-5s %-10s %-10s %-10s %-10s %-10s %-6s | %-10s \n", "Gen", "min_LL", "max_LL", "mean_LL","median_LL","sd_LL","NA_count","time end"))
  cat(sprintf("%-5d %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f %-6d | %s\n", g, best_fit, max_fit, mean_fit, median_fit, sd_fit, na_count, Sys.time()))
  
  #initial population.................................
  #message('Running the initial population')
  gapop <- fn.GApop(pardist=gapop.pardist,vuldist=gapop.vuldist)
  gapop[1,] <- est_par_vec  #include base run in initial population
  write.csv(gapop,file.path(dir.results,paste0('init_ga_pop',timestamp,'.csv')))
            
  pdf(file.path(dir.results,paste0('init_ga_pop_distributions_',timestamp,'.pdf')),onefile=T)
  graphics.off();rm(.SavedPlots);windows(record=T)
  par(mfrow=c(3,3))
  for(i in 1:ncol(gapop)){
    hist(gapop[,i],breaks=100,main=colnames(gapop)[i],xlab='value')
    abline(v=est_par_vec[i],col='blue',lty=2)
  }
  dev.off()

  files.cmd <- lapply(1:nrow(gapop),function(i) fn.parvec2cmd(par_vec=gapop[i,], g=0, idx=i)) 
  files.cmd <- unlist(files.cmd, use.names=F)
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=1, delete.out=T)

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
  na_count = sum(is.na(fitness))
  best_pars = round(gapop[which.min(fitness),],4)
  write.table(t(c(g,best_fit, max_fit, mean_fit, median_fit, sd_fit, na_count,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
  cat(sprintf("%-5d %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f %-6d | %s\n", g, best_fit, max_fit, mean_fit, median_fit, sd_fit, na_count, Sys.time()))
  #message(sprintf("Generation 0: Lowest LL score  = %.4f\n", min(fitness)))
  #message(sprintf("Generation 0: Avg pop LL score  = %.4f\n", mean(fitness)))
  
  for (gen in 1:n_generations) {
    #gen=1
    if(gen==1) unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)
    
    #reset workers after every 5 generations - this didn't help but we'll keep it jic
    if(n_generations>5){
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
    }
    # Elitism - keep the top n runs
    elite_idx <- order(fitness, decreasing=F)[1:elitism]
    elite <- gapop[elite_idx, ]
    
    # Selection, Crossover, Mutation
    parents <- select_parents(gapop, fitness) #resamples the population, with replacement, with rank-based probabilities in the sample draws
    #offspring <- crossover(parents) #offspring are when two parents crossover a part of their parameter vector
    offspring <- crossover_uniform(parents, group_id=par.groups, p_cross=0.8, p_group=0.5)
    offspring <- mutate(offspring, margin=mutate.margin) #randomly draw new parameter values to mutate the individual
    
    # graphics.off();rm(.SavedPlots);windows(record=T)
    # par(mfrow=c(3,3))
    # for(i in 1:ncol(offspring)){
    #   hist(offspring[,i],breaks=100,main=colnames(offspring)[i],xlab='value')
    #   abline(v=est_par_vec[i],col='blue',lty=2)
    # }

    
    # Evaluate new population
    ensure_cluster()  # optional but recommended
    files.cmd <- lapply(1:nrow(offspring),function(i) fn.parvec2cmd(par_vec=offspring[i,], g=gen, idx=i)) 
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
    
    #plot and save final generation
    if(gen==n_generations){
      pdf(file.path(dir.results,paste0('final_ga_pop_distributions_',timestamp,'.pdf')),onefile=T)
      graphics.off();rm(.SavedPlots);windows(record=T)
      par(mfrow=c(3,3))
      for(i in 1:ncol(gapop)){
        hist(gapop[,i],breaks=100,main=colnames(gapop)[i],xlab='value')
        abline(v=est_par_vec[i],col='blue',lty=2)
      }
      dev.off()
      write.csv(gapop,file.path(dir.results,paste0('final_ga_pop',timestamp,'.csv')))
    }
    
    
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
    na_count = sum(is.na(fitness))
    best_pars = round(gapop[which.min(fitness),],4)
    write.table(t(c(g,best_fit, max_fit, mean_fit, median_fit, sd_fit, na_count,best_pars)), file.ga.results, sep=",", append=T, row.names=F, col.names=F)
    cat(sprintf("%-5d %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f %-6d | %s\n", g, best_fit, max_fit, mean_fit, median_fit, sd_fit, na_count, Sys.time()))
    #cat(sprintf("%-10d %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f | %s\n", g, best_fit, max_fit, mean_fit, median_fit, sd_fit, Sys.time()))
    #message(sprintf("Generation %d: Avg pop LL score = %.4f\n", gen, mean(fitness)))
    #unlink(list.dirs(run_dir, full.names = T, recursive = F), recursive=T)
    rm(parents,offspring,files.cmd);gc()
  }
}#eof


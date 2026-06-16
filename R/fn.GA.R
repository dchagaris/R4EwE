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
  n_generations <<- myconfig$n_gen
  mutation_rate <<- myconfig$pmutation
  elitism <<- myconfig$elitism
  do.penalty <<- myconfig$do.penalty
  pen.wt.mult <<- myconfig$pen.wt.mult
  gapop.pardist <<- myconfig$gapop.pardist
  gapop.vuldist <<- myconfig$gapop.vuldist
  mutate.margin <<- myconfig$mutate.margin
  dir.results <<- myconfig$dir.ga_results

  # objective function selection (defaults: fn.ecospace_objfxn with fit.abs.catch=FALSE)
  obj.fxn       <- if(is.null(myconfig$obj.fxn))       3L    else as.integer(myconfig$obj.fxn)
  fit.abs.catch <- if(is.null(myconfig$fit.abs.catch)) FALSE else as.logical(myconfig$fit.abs.catch)

  # If TRUE, every individual in the final-generation gapop ends up with a preserved
  # Ecospace output directory so downstream per-run aggregation (e.g. fn.M0_loss_rate
  # across the final pop to derive Mrt uncertainty bounds) has full coverage.
  # Implementation:
  #   1) gen == n_generations runs with delete.output=FALSE, so the offspring runs
  #      that survive into gapop stay on disk at zero extra cost.
  #   2) The elite rows in the final gapop were carried over from gen N-1 and have
  #      no run on disk; they get re-evaluated in a small extra batch (~elitism
  #      runs) after the loop. EwE.exe is deterministic for a given param vector,
  #      so this regenerates the deterministic output rather than adding a new run.
  # A mapping CSV (final_ga_runs_<ts>.csv) records each gapop row's role and its
  # cmd file path. Default FALSE for backward compatibility.
  preserve.final.output <- isTRUE(myconfig$preserve.final.output)
  
  #create results file
  file.ga.results <- file.path(dir.results, paste0('ga_results_',timestamp,'.csv'))
  myconfig.string <- paste(paste(names(myconfig),unlist(myconfig),sep="="),collapse="; ")
  write.table('Ecospace GA calibration', file.ga.results,row.names=F, col.names=F,append=F)
  write.table(Sys.time(), file.ga.results,row.names=F, col.names=F,append=T)
  write.table(cmd_base[which(startsWith(cmd_base,"<EWE_MODEL_FILE>"))],         file.ga.results, row.names=F, col.names=F, append=T, quote=F)
  write.table(cmd_base[which(startsWith(cmd_base,"<ECOSIM_SCENARIO_INDEX>"))],  file.ga.results, row.names=F, col.names=F, append=T, quote=F)
  write.table(cmd_base[which(startsWith(cmd_base,"<ECOSPACE_SCENARIO_INDEX>"))],file.ga.results, row.names=F, col.names=F, append=T, quote=F)
  write.table(myconfig.string,                                                  file.ga.results, row.names=F, col.names=F, append=T, quote=F)
  write.table(t(c("gen","min_LL","max_LL","mean_LL","median_LL","sd_LL","NA_count",names(est_par_vec))), file.ga.results, sep=",", row.names = F, col.names = F, append=T)
  
  #base run
  files.cmd <- fn.parvec2cmd(par_vec=est_par_vec,g=999,idx=0)
  message('Running the base model')
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=obj.fxn, fit.abs.catch=fit.abs.catch, delete.output=T, gen=-1)

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
            
  fn.plot_ga_pop_responses(gapop, est_par_vec, par.groups, par.labels,
                           env_pars = env_pars, redtide_pars = redtide_pars,
                           file = file.path(dir.results,
                                  paste0('init_ga_pop_distributions_', timestamp, '.pdf')))

  files.cmd <- lapply(1:nrow(gapop),function(i) fn.parvec2cmd(par_vec=gapop[i,], g=0, idx=i))
  files.cmd <- unlist(files.cmd, use.names=F)
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=obj.fxn, fit.abs.catch=fit.abs.catch, delete.output=T, gen=0)

  #calculate penalty for parameter bounds violations
  if(do.penalty){
    pen.wt = abs(diff(range(fitness)))*pen.wt.mult
    upp.pen = U.bounds*1.5
    low.pen = L.bounds*0.5
    penalties = rep(0,nrow(gapop))
    for(i in 1:nrow(gapop)){
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
      # Cluster lives in .GlobalEnv as a single source of truth (also used by
      # ensure_cluster). Use cl <<- so this reset writes back to the same slot
      # rather than creating a function-local shadow that leaks zombie workers.
      if(exists("cl", envir = .GlobalEnv)){
        try(stopCluster(get("cl", envir = .GlobalEnv)), silent = TRUE)
        rm(list = "cl", envir = .GlobalEnv)
      }
      gc()
      workers <- if(!is.null(myconfig$workers)) myconfig$workers else floor(detectCores() / 2)
      cl <<- parallel::makePSOCKcluster(workers, outfile = "cluster_workers.log")
      doSNOW::registerDoSNOW(cl)
      clusterExport(cl, .r4ewe_worker_exports())
    }
    }
    # Elitism - keep the top n runs
    elite_idx <- order(fitness, decreasing=F)[1:elitism]
    elite <- gapop[elite_idx, , drop = FALSE]
    
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
    ensure_cluster(workers = myconfig$workers)  # optional but recommended
    files.cmd <- lapply(1:nrow(offspring),function(i) fn.parvec2cmd(par_vec=offspring[i,], g=gen, idx=i))
    files.cmd <- unlist(files.cmd, use.names=F)
    # On the last gen, keep the offspring run dirs so downstream per-run aggregation
    # (e.g. Mrt across the final pop) can read each individual's Ecospace output.
    del.this.gen <- !(preserve.final.output && gen == n_generations)
    new_fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=obj.fxn, fit.abs.catch=fit.abs.catch, delete.output=del.this.gen, gen=gen)
    
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
    
    # Combine elite + offspring: keep elites and the best (lowest) offspring.
    # na.last=FALSE pushes NA fitnesses to the front of the decreasing sort so they're
    # included in worst_offspring_idx and dropped (failed runs shouldn't propagate).
    worst_offspring_idx <- order(new_fitness, decreasing = TRUE, na.last = FALSE)[1:elitism]
    gapop   <- rbind(elite, offspring[-worst_offspring_idx, ])
    fitness <- c(fitness[elite_idx], new_fitness[-worst_offspring_idx])
    
    #plot and save final generation
    if(gen==n_generations){
      # Posterior plot: highlight the min-LL individual (best fit) in blue rather than
      # the baseline est_par_vec. Histograms get a vline at best-fit values, env/rt
      # panels overlay the best-fit response curve in blue.
      best_pars_final <- gapop[which.min(fitness), ]
      fn.plot_ga_pop_responses(gapop, est_par_vec, par.groups, par.labels,
                               env_pars = env_pars, redtide_pars = redtide_pars,
                               highlight_pars = best_pars_final,
                               file = file.path(dir.results,
                                      paste0('final_ga_pop_distributions_', timestamp, '.pdf')))
      write.csv(gapop, file.path(dir.results, paste0('final_ga_pop', timestamp, '.csv')))

      if(preserve.final.output){
        # gapop layout after rebuild: rows 1..elitism are elites (carried over from
        # gen N-1, no output dir on disk); rows elitism+1..popSize are this gen's
        # surviving offspring (output dirs kept because del.this.gen=FALSE above).
        offspring_rows_kept <- setdiff(seq_len(nrow(offspring)), worst_offspring_idx)
        offspring_cmd       <- files.cmd[offspring_rows_kept]

        # Re-run the elite slots to materialize their output dirs. EwE.exe is
        # deterministic so this regenerates the same outputs gen N-1 produced.
        elite_cmd <- character(0L)
        if(elitism > 0L){
          elite_cmd <- vapply(seq_len(elitism),
                              function(i) fn.parvec2cmd(par_vec=gapop[i,], g=998, idx=i),
                              character(1L))
          cat(sprintf("\n[final-pop] Re-running %d elite individuals to preserve their Ecospace output\n", elitism))
          .ignored <- fn.runEwE.gapop(elite_cmd, obj.fxn=obj.fxn,
                                      fit.abs.catch=fit.abs.catch,
                                      delete.output=FALSE, gen="elite")
        }

        mapping <- data.frame(
          gapop_row = seq_len(nrow(gapop)),
          role      = c(rep("elite",     elitism),
                        rep("offspring", length(offspring_cmd))),
          cmd_file  = c(elite_cmd, offspring_cmd),
          fitness   = fitness,
          stringsAsFactors = FALSE
        )
        map_path <- file.path(dir.results, paste0("final_ga_runs_", timestamp, ".csv"))
        write.csv(mapping, map_path, row.names=FALSE)
        message(sprintf("Final-pop output preserved: %d/%d individuals on disk; mapping at %s",
                        sum(!is.na(mapping$cmd_file)), nrow(gapop), map_path))
      }
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


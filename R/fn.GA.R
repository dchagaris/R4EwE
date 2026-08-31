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
#' @export
fn.GA <- function(myconfig){
  n_generations <<- myconfig$n_gen
  mutation_rate <<- myconfig$pmutation
  elitism <<- myconfig$elitism
  do.penalty <<- myconfig$do.penalty
  pen.wt.mult <<- myconfig$pen.wt.mult
  gapop.pardist <<- myconfig$gapop.pardist
  gapop.vuldist <<- myconfig$gapop.vuldist
  gapop.rtdist  <<- if(is.null(myconfig$gapop.rtdist)) "log-uniform" else myconfig$gapop.rtdist
  gapop.fltdist <<- if(is.null(myconfig$gapop.fltdist)) "normal"      else myconfig$gapop.fltdist
  selection.method <<- if(is.null(myconfig$selection.method)) "tournament" else myconfig$selection.method
  tournament.k     <<- if(is.null(myconfig$tournament.k))     3L           else as.integer(myconfig$tournament.k)
  dir.results <<- myconfig$dir.ga_results

  # objective function selection (defaults: fn.ecospace_objfxn with fit.abs.catch=FALSE)
  obj.fxn       <- if(is.null(myconfig$obj.fxn))       3L    else as.integer(myconfig$obj.fxn)
  fit.abs.catch <- if(is.null(myconfig$fit.abs.catch)) FALSE else as.logical(myconfig$fit.abs.catch)

  # Phase-3 elite-carry design: every generation's runs live in a per-gen
  # subfolder run_dir/gen<NN>/. At each gen boundary, the top-`elitism` dirs
  # are MOVED into the next gen's folder (folder names unchanged for
  # traceability back to gen-of-origin), and the rest of the previous gen's
  # dirs are deleted. Elite output on disk is the ORIGINAL evaluation that
  # earned the fitness -- no reruns, no stochastic drift. The final-gen
  # folder therefore already holds all `popSize` dirs, and the mapping CSV
  # is written from the tracked cmd.files array (no g=998 rerun step).
  # `preserve.final.output` is now effectively always TRUE and the flag is
  # retained only for backward compatibility (ignored).
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
  
  #base run -- run the model with cmd_base alone (no appended parameter taglines).
  # cmd_base already carries the user-tuned baseline tags; appending another set
  # from est_par_vec would override those with sensitivity-derived base.val
  # values that may not match the cmd_base baselines, distorting gen=-1.
  base_run_dir  <- file.path(run_dir, "run_base")
  if(!dir.exists(base_run_dir)) dir.create(base_run_dir, recursive = TRUE, showWarnings = FALSE)
  base_out_dir  <- normalizePath(base_run_dir, winslash = "\\", mustWork = FALSE)
  cmd_base_run  <- cmd_base
  cmd_base_run[startsWith(cmd_base_run, "<ECOSPACE_OUTPUT_DIR>")] <-
    sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", base_out_dir)
  base_cmd_file <- file.path(base_run_dir, "cmd.txt")
  writeLines(cmd_base_run, base_cmd_file)

  message('Running the base model')
  # Keep base run on disk for reference/comparison. It lives in run_dir/run_base/
  # and is never touched by the gen-loop cleanup (which only operates on
  # run_dir/gen<NN>/ folders).
  fitness <- fn.runEwE.gapop(base_cmd_file, obj.fxn=obj.fxn, fit.abs.catch=fit.abs.catch, delete.output=F, gen=-1)

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
  gapop <- fn.GApop(pardist=gapop.pardist, vuldist=gapop.vuldist,
                    rtdist=gapop.rtdist, fltdist=gapop.fltdist)
  gapop[1,] <- est_par_vec  #include base run in initial population
  write.csv(gapop,file.path(dir.results,paste0('init_ga_pop',timestamp,'.csv')))
            
  fn.plot_ga_pop_responses(gapop, est_par_vec, par.groups, par.labels,
                           env_pars = env_pars, redtide_pars = redtide_pars,
                           file = file.path(dir.results,
                                  paste0('init_ga_pop_distributions_', timestamp, '.pdf')))

  # Initial population goes into run_dir/gen00/. Keep on disk (delete.output=F)
  # so gen-1's elite carry can move the top-`elitism` dirs into run_dir/gen01/.
  files.cmd <- lapply(1:nrow(gapop),function(i) fn.parvec2cmd(par_vec=gapop[i,], g=0, idx=i))
  files.cmd <- unlist(files.cmd, use.names=F)
  fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=obj.fxn, fit.abs.catch=fit.abs.catch, delete.output=F, gen=0)
  # Sync files.cmd + gapop with any retry substitutions (fresh random parvecs
  # replacing failed originals). Keeps gapop[i,] <-> files.cmd[i] alignment.
  if(!is.null(attr(fitness, "cmd.files")))
    files.cmd <- attr(fitness, "cmd.files")
  subs <- attr(fitness, "substituted")
  if(!is.null(subs) && length(subs) > 0L){
    for(idx_str in names(subs))
      gapop[as.integer(idx_str), ] <- subs[[idx_str]]
  }
  # cmd.files tracks each gapop row's cmd path across gens. Rebuilt after every
  # gen boundary as c(elite.cmd.moved, this-gen's surviving offspring cmds).
  cmd.files <- files.cmd

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
    # NOTE: no blanket unlink here. Gen dirs are managed explicitly below
    # (elites moved forward, dropped offspring deleted after fitness known).

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
    elite     <- gapop[elite_idx, , drop = FALSE]
    elite.cmd <- cmd.files[elite_idx]   # cmd paths in prior gen's folder

    # Move elite dirs into this gen's folder (folder names unchanged for
    # traceability back to the elite's gen-of-origin + offspring idx).
    # Cmd files' <ECOSPACE_OUTPUT_DIR> is rewritten to the new path.
    dst_gen_dir     <- .gen_dir(run_dir, gen)
    elite.cmd.moved <- .move_elite_dirs(elite.cmd, dst_gen_dir)

    # Now delete the previous gen's folder (all surviving offspring from that
    # gen that were NOT promoted to elite are no longer needed).
    prev_gen_dir <- .gen_dir(run_dir, gen - 1L)
    if(dir.exists(prev_gen_dir))
      unlink(prev_gen_dir, recursive = TRUE)

    # Selection, Crossover, Mutation
    parents <- if(identical(selection.method, "rank")){
      select_parents(gapop, fitness)
    } else {
      tournament_select(gapop, fitness, k = tournament.k)
    }
    offspring <- crossover_uniform(parents, group_id=par.groups, p_cross=0.8, p_group=0.3)
    offspring <- mutate(offspring) #randomly draw new parameter values to mutate the individual
    
    # graphics.off();rm(.SavedPlots);windows(record=T)
    # par(mfrow=c(3,3))
    # for(i in 1:ncol(offspring)){
    #   hist(offspring[,i],breaks=100,main=colnames(offspring)[i],xlab='value')
    #   abline(v=est_par_vec[i],col='blue',lty=2)
    # }

    
    # Evaluate new population -- always keep offspring dirs on disk. Selective
    # cleanup of the dropped-worst offspring happens after fitness is known
    # (see below), and elite dirs persist across gens via .move_elite_dirs.
    ensure_cluster(workers = myconfig$workers)  # optional but recommended
    files.cmd   <- lapply(1:nrow(offspring),function(i) fn.parvec2cmd(par_vec=offspring[i,], g=gen, idx=i))
    files.cmd   <- unlist(files.cmd, use.names=F)
    new_fitness <- fn.runEwE.gapop(files.cmd, obj.fxn=obj.fxn, fit.abs.catch=fit.abs.catch, delete.output=F, gen=gen)
    # Sync files.cmd + offspring with any retry substitutions (fresh random
    # parvecs that replaced failed originals). Updating `offspring` here means
    # the substituted parvec flows into the new gapop via the rbind below --
    # so the population matrix stays consistent with what actually ran.
    if(!is.null(attr(new_fitness, "cmd.files")))
      files.cmd <- attr(new_fitness, "cmd.files")
    subs <- attr(new_fitness, "substituted")
    if(!is.null(subs) && length(subs) > 0L){
      for(idx_str in names(subs))
        offspring[as.integer(idx_str), ] <- subs[[idx_str]]
    }
    
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
    
    # Combine elite + offspring: drop ALL NA-fitness offspring (failed runs should
    # never propagate, even if there are more of them than elite slots), then drop
    # enough additional worst non-NA offspring to free up `elitism` slots. If NAs
    # alone already exceed elitism, the resulting population is short of popSize
    # and we backfill with elite clones (with a warning); this is extremely rare
    # given fn.runEwE.gapop's retry loop, but it keeps gapop dimensions stable.
    na_idx     <- which(is.na(new_fitness))
    non_na_idx <- which(!is.na(new_fitness))
    n_drop     <- max(elitism, length(na_idx))
    worst_non_na <- non_na_idx[order(new_fitness[non_na_idx], decreasing = TRUE)]
    drop_extra <- worst_non_na[seq_len(max(0L, n_drop - length(na_idx)))]
    worst_offspring_idx <- c(na_idx, drop_extra)

    # Delete dropped-offspring dirs from this gen's folder before rebuilding.
    dropped_cmd <- files.cmd[worst_offspring_idx]
    if(length(dropped_cmd) > 0L)
      unlink(dirname(dropped_cmd), recursive = TRUE)

    gapop     <- rbind(elite, offspring[-worst_offspring_idx, , drop=FALSE])
    fitness   <- c(fitness[elite_idx], new_fitness[-worst_offspring_idx])
    cmd.files <- c(elite.cmd.moved, files.cmd[-worst_offspring_idx])

    deficit <- myconfig$popSize - nrow(gapop)
    if(deficit > 0L){
      warning(sprintf("Gen %d: %d offspring failed beyond what elites can fill; backfilling with %d elite clones.",
                      gen, length(na_idx), deficit))
      fill_idx  <- rep(seq_len(elitism), length.out = deficit)
      gapop     <- rbind(gapop,   elite[fill_idx, , drop=FALSE])
      fitness   <- c(fitness, fitness[seq_len(elitism)][fill_idx])
      # Backfilled slots point at the same elite dir on disk (harmless: elite
      # is a genome duplicate). Downstream ensemble should dedup on cmd path.
      cmd.files <- c(cmd.files, elite.cmd.moved[fill_idx])
    }
    
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

      # Elite dirs have been carried forward via .move_elite_dirs on every gen
      # boundary, so run_dir/gen<n_generations>/ already holds all `popSize`
      # dirs. The mapping CSV is written directly from the tracked cmd.files
      # array -- no g=998 rerun step needed.
      mapping <- data.frame(
        gapop_row = seq_len(nrow(gapop)),
        role      = c(rep("elite",     elitism),
                      rep("offspring", nrow(gapop) - elitism)),
        cmd_file  = cmd.files,
        fitness   = fitness,
        stringsAsFactors = FALSE
      )
      map_path <- file.path(dir.results, paste0("final_ga_runs_", timestamp, ".csv"))
      write.csv(mapping, map_path, row.names=FALSE)
      message(sprintf("Final-pop mapping written: %d individuals on disk in %s\n  %s",
                      nrow(mapping), .gen_dir(run_dir, gen), map_path))
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


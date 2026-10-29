#GA function----
#' @title Run genetic algorithm calibration for Ecospace
#' @description
#' Sets GA controls from `myconfig`, writes a results CSV, runs a base model,
#' evaluates an initial population, then iterates generations with elitism,
#' selection, crossover, mutation, and optional parameter-bound penalties.
#' Progress and summary stats are printed each generation.
#'
#' Stopping rule: from generation \code{stop.min_gen} on, the GA stops when either
#' (1) the best LL has improved by less than the tolerance over the last
#' \code{stop.patience} generations, or (2) the SD of the population LL is below the
#' tolerance (the population has collapsed onto near-clones). Otherwise it stops at
#' \code{n_gen}. Either way the final-population files are written for the
#' generation it stopped at, and the reason is written to
#' \code{ga_stop_<timestamp>.csv}.
#'
#' @param myconfig List of GA settings, expected fields:
#'   - popSize: population size
#'   - n_gen: maximum number of generations, a hard cap (default 30)
#'   - stop.rtol: relative stopping tolerance (default 0.005). The tolerance in LL
#'     units is \code{stop.rtol * |best LL|}, recomputed each generation. 0 turns the
#'     early stop off, so the GA runs all \code{n_gen} generations.
#'   - stop.patience: generations over which the best LL must improve by at least
#'     the tolerance (default 3)
#'   - stop.min_gen: earliest generation at which the GA may stop (default 5; raised
#'     to \code{stop.patience} if smaller)
#'   - pmutation: per-gene mutation probability
#'   - elitism: number of elites kept each generation
#'   - do.penalty: logical, apply penalties for bound violations
#'   - pen.wt.mult: scalar multiplier for penalty weight
#'   - gapop.dist: initial population distribution ("unif" or "ln")
#'   - resume.gen: NULL (default) for a fresh run. To continue a run that died, the
#'     last generation it completed (or TRUE for the last one in \code{ga_evals}).
#'     The base and gen-0 runs are skipped; see Resume below.
#'
#' Resume: the population, fitness and best-LL history at the end of
#' \code{resume.gen} are rebuilt from \code{ga_evals_<timestamp>.csv} (gen 0, then
#' each generation's top \code{elitism} plus its kept offspring) and checked against
#' the min / median / SD rows in \code{ga_results_<timestamp>.csv}. Only the elites
#' need run dirs; they are found in \code{run_dir/gen<resume.gen+1>/} or
#' \code{gen<resume.gen>/}. Run dirs left in \code{gen<resume.gen+1>/} by the crashed
#' attempt are deleted (their parameter vectors were never logged). The GA then
#' continues from generation \code{resume.gen + 1}, appending to both CSVs. The
#' caller must rebuild the same \code{est_par_vec} and bounds (same sensitivity, same
#' script), and use the same \code{timestamp}, \code{run_dir}, \code{popSize} and
#' \code{elitism}.
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
#' @return Invisibly returns a list with \code{final_gen} and \code{stop_reason}. Results (fitness summaries and best parameters)
#' are appended to the CSV and printed to the console. Every evaluated individual
#' is appended to \code{ga_evals_<timestamp>.csv} in \code{dir.ga_results}: one row
#' per Ecospace run with gen, idx, role (base / base-in-pop / initial / offspring),
#' kept (FALSE for offspring dropped at the generation boundary), fitness, run
#' folder and every parameter value.
#' @export
fn.GA <- function(myconfig){
  n_generations <<- if(is.null(myconfig$n_gen)) 30L else as.integer(myconfig$n_gen)
  stop.rtol     <- if(is.null(myconfig$stop.rtol))     0.005 else as.numeric(myconfig$stop.rtol)
  stop.patience <- if(is.null(myconfig$stop.patience)) 3L    else as.integer(myconfig$stop.patience)
  stop.min_gen  <- if(is.null(myconfig$stop.min_gen))  5L    else as.integer(myconfig$stop.min_gen)
  stop.min_gen  <- max(stop.min_gen, stop.patience)   # need stop.patience gens of history
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

  file.ga.results <- file.path(dir.results, paste0('ga_results_',timestamp,'.csv'))
  file.ga.evals   <- file.path(dir.results, paste0('ga_evals_', timestamp, '.csv'))
  resume.gen <- myconfig$resume.gen
  if(isTRUE(resume.gen)) resume.gen <- max(utils::read.csv(file.ga.evals, check.names = FALSE)$gen)

  if(!is.null(resume.gen)){
    # Resume: rebuild the population at the end of resume.gen (see Resume in the docs)
    if(do.penalty) stop("fn.GA: resume.gen is not supported with do.penalty = TRUE")
    st <- .ga_resume_state(file.ga.evals, file.ga.results, as.integer(resume.gen), elitism,
                           myconfig$popSize, run_dir, names(est_par_vec))
    gapop       <- st$pop
    fitness     <- st$fit
    cmd.files   <- st$cmd
    best_hist   <- st$best_hist
    stop.reason <- NULL
    first.gen   <- as.integer(resume.gen) + 1L
    message(sprintf("Resuming after gen %d: best LL %.2f, %d elite dirs found; continuing at gen %d",
                    resume.gen, min(fitness, na.rm = TRUE), sum(!is.na(cmd.files)), first.gen))
  } else {
    first.gen <- 1L

    #create results file
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
    # every evaluated individual, one row per Ecospace run (see .ga_log_evals)
    file.ga.evals <- file.path(dir.results, paste0('ga_evals_', timestamp, '.csv'))
    if(file.exists(file.ga.evals)) unlink(file.ga.evals)
    .ga_log_evals(file.ga.evals, gen = -1, pop = matrix(est_par_vec, nrow = 1, dimnames = list(NULL, names(est_par_vec))),
                  fit = fitness[1], cmd = base_cmd_file, kept = TRUE, role = "base")
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
    .ga_log_evals(file.ga.evals, gen = 0, pop = gapop, fit = fitness, cmd = unlist(cmd.files),
                  kept = TRUE, role = ifelse(seq_len(nrow(gapop)) == 1, "base-in-pop", "initial"))
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
    best_hist   <- best_fit   # best LL by gen, gen 0 first (elites keep it non-increasing)
    stop.reason <- NULL
    #message(sprintf("Generation 0: Avg pop LL score  = %.4f\n", mean(fitness)))
  }

  for (gen in first.gen:n_generations) {
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
      fn.cluster_export(cl)
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
    .ga_log_evals(file.ga.evals, gen = gen, pop = offspring, fit = new_fitness, cmd = files.cmd,
                  kept = !(seq_len(nrow(offspring)) %in% worst_offspring_idx), role = "offspring")

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
    
    # Stopping rule (see Details): tolerance = stop.rtol x |current best LL|
    best_now  <- min(fitness, na.rm = TRUE)
    best_hist <- c(best_hist, best_now)
    stop.tol  <- stop.rtol * abs(best_now)
    stop.gain <- best_hist[max(1L, gen + 1L - stop.patience)] - best_now   # best_hist[1] is gen 0
    stop.sd   <- sd(fitness, na.rm = TRUE)
    if(stop.rtol > 0 && gen >= stop.min_gen){
      if(stop.gain < stop.tol){
        stop.reason <- sprintf("best LL improved %.1f over the last %d gens (< tol %.1f)",
                               stop.gain, stop.patience, stop.tol)
      } else if(stop.sd < stop.tol){
        stop.reason <- sprintf("population LL SD %.1f < tol %.1f", stop.sd, stop.tol)
      }
    }
    if(is.null(stop.reason) && gen == n_generations)
      stop.reason <- sprintf("reached n_gen = %d", n_generations)
    is.final <- !is.null(stop.reason)

    #plot and save final generation
    if(is.final){
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
    if(is.final){
      message(sprintf("GA stopped at gen %d: %s", gen, stop.reason))
      utils::write.csv(data.frame(gen = gen, reason = stop.reason, best_LL = best_now,
                                  gain = stop.gain, sd_LL = stop.sd, tol = stop.tol,
                                  stop.rtol = stop.rtol, stop.patience = stop.patience,
                                  stop.min_gen = stop.min_gen, n_gen = n_generations),
                       file.path(dir.results, paste0("ga_stop_", timestamp, ".csv")), row.names = FALSE)
      break
    }
  }
  invisible(list(final_gen = gen, stop_reason = stop.reason))
}#eof


# Internal: append one row per evaluated individual to the GA evaluation log
# (ga_evals_<timestamp>.csv): gen, idx, role, kept, fitness, run folder, then
# every parameter value. Offspring dropped at the generation boundary are kept
# in the log (kept = FALSE) -- they are the ones that show which parameter
# ranges produce poor or catastrophic fits. Elites are not re-run, so each
# individual appears once, at the generation it was evaluated.
#' @keywords internal
#' @noRd
.ga_log_evals <- function(file, gen, pop, fit, cmd, kept, role){
  pop <- as.matrix(pop)
  n   <- nrow(pop)
  out <- data.frame(gen = gen, idx = seq_len(n), role = rep_len(role, n), kept = rep_len(kept, n),
                    fitness = as.numeric(fit)[seq_len(n)],
                    run = basename(dirname(as.character(unlist(cmd))))[seq_len(n)],
                    pop, check.names = FALSE, stringsAsFactors = FALSE)
  first <- !file.exists(file)
  utils::write.table(out, file, sep = ",", row.names = FALSE, col.names = first, append = !first)
  invisible(out)
}


# Internal: rebuild the GA state at the end of generation `resume.gen` from the
# evaluation log, replaying fn.GA's generation boundary: the new population is the
# top `elitism` of the previous one (order(), as in fn.GA) followed by that gen's
# kept offspring, backfilled with elite clones when short. Checks each rebuilt
# generation against the min / median / SD rows of ga_results, finds the elite
# dirs on disk, and deletes run dirs that a crashed attempt left in
# gen<resume.gen+1>/ (their parameters were never logged).
#' @keywords internal
#' @noRd
.ga_resume_state <- function(file.evals, file.results, resume.gen, elitism, popSize, run_dir, par.names){
  ev <- utils::read.csv(file.evals, check.names = FALSE, stringsAsFactors = FALSE)
  miss <- setdiff(par.names, names(ev))
  if(length(miss) || length(setdiff(setdiff(names(ev), c("gen", "idx", "role", "kept", "fitness", "run")), par.names)))
    stop("resume: est_par_vec does not match the parameters in ", basename(file.evals),
         if(length(miss)) paste0(" (missing: ", paste(head(miss, 5), collapse = ", "), ")"))
  if(!resume.gen %in% ev$gen) stop("resume: gen ", resume.gen, " not found in ", basename(file.evals))
  l  <- readLines(file.results)
  rs <- utils::read.csv(text = l[grep("^\"?gen", l)[1]:length(l)], check.names = FALSE)

  g0  <- ev[ev$gen == 0, ]; g0 <- g0[order(g0$idx), ]
  pop <- as.matrix(g0[, par.names]); fit <- g0$fitness; run <- g0$run
  best_hist <- min(fit, na.rm = TRUE)
  chk <- function(g){
    r <- rs[rs$gen == g, ]
    if(nrow(r) == 0) stop("resume: gen ", g, " has no row in ", basename(file.results))
    d <- abs(c(min(fit, na.rm = TRUE), stats::median(fit, na.rm = TRUE), stats::sd(fit, na.rm = TRUE)) -
             c(r$min_LL[1], r$median_LL[1], r$sd_LL[1]))
    if(any(d > 0.01)) stop(sprintf("resume: rebuilt gen %d does not match %s (min/median/sd off by %s)",
                                   g, basename(file.results), paste(signif(d, 3), collapse = "/")))
  }
  chk(0)
  for(g in seq_len(resume.gen)){
    o <- ev[ev$gen == g & ev$role == "offspring", ]; o <- o[order(o$idx), ]
    k <- as.logical(o$kept)
    ei  <- order(fit, decreasing = FALSE)[seq_len(elitism)]
    e.pop <- pop[ei, , drop = FALSE]; e.fit <- fit[ei]; e.run <- run[ei]
    pop <- rbind(e.pop, as.matrix(o[k, par.names])); fit <- c(e.fit, o$fitness[k]); run <- c(e.run, o$run[k])
    deficit <- popSize - nrow(pop)
    if(deficit > 0L){
      f <- rep(seq_len(elitism), length.out = deficit)
      pop <- rbind(pop, e.pop[f, , drop = FALSE]); fit <- c(fit, e.fit[f]); run <- c(run, e.run[f])
    }
    chk(g)
    best_hist <- c(best_hist, min(fit, na.rm = TRUE))
  }
  rownames(pop) <- NULL

  # elite dirs: moved into gen<resume.gen+1>/ if that gen had started, else still in gen<resume.gen>/
  nxt <- .gen_dir(run_dir, resume.gen + 1L); cur <- .gen_dir(run_dir, resume.gen)
  ei  <- order(fit, decreasing = FALSE)[seq_len(elitism)]
  cmd <- rep(NA_character_, length(run))
  for(i in ei){
    f <- c(file.path(nxt, run[i], "cmd.txt"), file.path(cur, run[i], "cmd.txt"))
    f <- f[file.exists(f)]
    if(!length(f)) stop("resume: elite run dir not found in ", basename(nxt), " or ", basename(cur), ": ", run[i])
    cmd[i] <- f[1]
  }
  if(dir.exists(nxt)){
    stale <- setdiff(list.dirs(nxt, full.names = FALSE, recursive = FALSE), run[ei])
    if(length(stale)){
      message(sprintf("resume: deleting %d run dirs left in %s by the crashed attempt", length(stale), nxt))
      unlink(file.path(nxt, stale), recursive = TRUE)
    }
  }
  list(pop = pop, fit = fit, cmd = cmd, best_hist = best_hist)
}

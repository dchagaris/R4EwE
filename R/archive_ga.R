#' @title Archive a completed GA calibration run
#' @description Minimum-viable phase-agnostic archive for a completed
#'   \code{fn.GA} run whose final-pop hook (\code{preserve.final.output=TRUE})
#'   has already written \code{final_ga_pop<ts>.csv} and
#'   \code{final_ga_runs_<ts>.csv} to \code{dir.results.ga}. Does four things:
#'   \enumerate{
#'     \item reads the final-pop CSVs to identify the best-fit individual,
#'     \item copies calibration script snapshots (and the R4EwE git rev) to
#'           \code{phase<N>_<ts>/code/},
#'     \item extracts absolute-value GUI-ready CSVs (vuls, disp, fleetdyn,
#'           env+redtide) from the best-fit parvec / cmd file and writes them
#'           to \code{phase<N>_<ts>/gui_ready/},
#'     \item writes a README.md summarising config, best LL, and known caveats.
#'   }
#'   Creates the phase-1 subdir structure (\code{code/}, \code{results/},
#'   \code{plots/}, \code{gui_ready/}, \code{logs/}) for consistency; only
#'   \code{code/} and \code{gui_ready/} are populated. Bundle the rest by hand
#'   or extend with \code{bundle_results}/\code{bundle_plots}/\code{bundle_logs}
#'   in a future revision.
#' @param phase Integer phase number (1, 2, 3, ...). Becomes the archive-dir
#'   prefix (\code{phase<N>_<timestamp>/}) and the README title.
#' @param timestamp GA-run timestamp string (matches the \code{<ts>} suffix on
#'   \code{final_ga_pop<ts>.csv}). Defaults to \code{timestamp} in
#'   \code{.GlobalEnv}.
#' @param dir.results.ga Directory holding \code{final_ga_pop<ts>.csv},
#'   \code{final_ga_runs_<ts>.csv}, and \code{ga_results_<ts>.csv}. Defaults to
#'   \code{dir.results.ga} in \code{.GlobalEnv}.
#' @param code_files Character vector of script paths to copy into
#'   \code{code/}. Relative paths are resolved against
#'   \code{dirname(dir.results.ga)}. Missing files are skipped silently.
#'   Default = the standard four calibration scripts.
#' @param r4ewe_dir Path to the R4EwE git working tree for commit-hash
#'   capture. Set to \code{NULL} to skip.
#' @param notes Optional character vector of extra bullet lines appended to a
#'   trailing \code{## Additional notes} section in the README. Use for
#'   phase-specific caveats.
#' @param overwrite If \code{FALSE} (default) and the target archive dir already
#'   exists, error out. Set \code{TRUE} to reuse and overwrite.
#' @param move_results If \code{TRUE} (default), moves every timestamp-matching
#'   file out of \code{dir.results.ga} into the archive after reads complete:
#'   CSVs to \code{results/}, PDFs/PNGs to \code{plots/}. Set \code{FALSE} to
#'   leave the source files in \code{dir.results.ga} (copy semantics were the
#'   old phase-1 behavior).
#' @return Invisibly returns the archive directory path.
#' @export
archive_ga <- function(phase,
                       timestamp      = NULL,
                       dir.results.ga = NULL,
                       code_files     = NULL,
                       r4ewe_dir      = "C:/dchagaris/GitHub/R4EwE",
                       notes          = NULL,
                       overwrite      = FALSE,
                       move_results   = TRUE){

  # ---- resolve inputs from .GlobalEnv if not supplied ----
  .g <- function(nm)
    if(exists(nm, envir = .GlobalEnv)) get(nm, envir = .GlobalEnv) else NULL

  if(is.null(timestamp))      timestamp      <- .g("timestamp")
  if(is.null(dir.results.ga)) dir.results.ga <- .g("dir.results.ga")
  if(is.null(code_files))     code_files     <- c(
    "setup cal WFS MICE.R",
    "run ecospace GA WFS MICE.R",
    "run par sensitivity WFS MICE.R"
  )
  stopifnot(!is.null(timestamp), !is.null(dir.results.ga),
            length(phase) == 1L, phase == as.integer(phase))
  phase <- as.integer(phase)

  # calibration-side globals needed for GUI-ready CSV construction
  est_par_vec   <- .g("est_par_vec")
  vuls.base     <- .g("vuls.base")
  disp.base     <- .g("disp.base")
  fleetdyn.base <- .g("fleetdyn.base")
  fit.vuls      <- isTRUE(.g("fit.vuls"))
  fit.disp      <- isTRUE(.g("fit.disp"))
  fit.fltdyn    <- isTRUE(.g("fit.fltdyn"))
  fit.env       <- isTRUE(.g("fit.env"))
  fit.redtide   <- isTRUE(.g("fit.redtide"))
  myconfig      <- .g("myconfig")

  # ---- Section 0. Locate inputs from the final-pop hook ----
  gapop_csv <- file.path(dir.results.ga, paste0("final_ga_pop",    timestamp, ".csv"))
  runs_csv  <- file.path(dir.results.ga, paste0("final_ga_runs_",  timestamp, ".csv"))
  ga_csv    <- file.path(dir.results.ga, paste0("ga_results_",     timestamp, ".csv"))
  if(!file.exists(gapop_csv)) stop("Missing: ", gapop_csv)
  if(!file.exists(runs_csv))
    stop("Missing: ", runs_csv,
         "\n(fn.GA's final-pop hook writes it when preserve.final.output=TRUE.)")

  gapop.final   <- read.csv(gapop_csv, row.names = 1, check.names = FALSE)
  runs_map      <- read.csv(runs_csv,  stringsAsFactors = FALSE)
  fitness.final <- runs_map$fitness
  best.idx      <- which.min(fitness.final)
  best.cmd      <- runs_map$cmd_file[best.idx]
  best.role     <- runs_map$role[best.idx]
  if(!file.exists(best.cmd))
    warning("Best-fit cmd file not on disk (path in mapping CSV): ", best.cmd)

  cat(sprintf("Phase-%d archive for GA run %s\n", phase, timestamp))
  cat(sprintf("  best-fit row : %d (%s)\n", best.idx, best.role))
  cat(sprintf("  best LL      : %.2f\n", fitness.final[best.idx]))

  # ---- Section 1. Build archive tree ----
  archive_dir <- file.path(dir.results.ga, sprintf("phase%d_%s", phase, timestamp))
  if(dir.exists(archive_dir) && !overwrite)
    stop("Archive already exists: ", archive_dir,
         "\nDelete/rename it or pass overwrite=TRUE.")
  for(sd in c("code", "results", "plots", "gui_ready", "logs"))
    dir.create(file.path(archive_dir, sd), recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("Archive tree: %s\n", archive_dir))

  # ---- Section 2. Code snapshots + R4EwE git rev ----
  dir.main <- dirname(dir.results.ga)
  copied   <- character(0)
  for(f in code_files){
    src <- if(file.exists(f)) f else file.path(dir.main, f)
    if(file.exists(src)){
      file.copy(src, file.path(archive_dir, "code", basename(src)), overwrite = TRUE)
      copied <- c(copied, basename(src))
    }
  }
  cat(sprintf("Copied %d code files.\n", length(copied)))

  if(!is.null(r4ewe_dir) && dir.exists(file.path(r4ewe_dir, ".git"))){
    tc <- function(cmd, args) tryCatch(
      system2("git", c("-C", shQuote(r4ewe_dir), cmd, args),
              stdout = TRUE, stderr = TRUE),
      error = function(e) sprintf("[git failed: %s]", conditionMessage(e)))
    commit   <- tc("rev-parse", "HEAD")
    branch   <- tc("rev-parse", c("--abbrev-ref", "HEAD"))
    status   <- tc("status",    "--porcelain")
    diff_out <- tc("diff",      "HEAD")
    writeLines(
      c(sprintf("R4EwE repo   : %s", r4ewe_dir),
        sprintf("branch       : %s", paste(branch, collapse = "")),
        sprintf("commit       : %s", paste(commit, collapse = "")),
        sprintf("archived at  : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "",
        "working-tree status (git status --porcelain):",
        if(length(status)) status else "(clean)"),
      file.path(archive_dir, "code", "R4EwE_commit.txt"))
    if(length(diff_out) && any(nzchar(diff_out)))
      writeLines(diff_out, file.path(archive_dir, "code", "R4EwE_diff.patch"))
    cat("Recorded R4EwE git commit + diff.\n")
  }

  # ---- Section 3. GUI-ready CSVs (absolute values) ----
  gui_dir    <- file.path(archive_dir, "gui_ready")
  ga.bestfit <- gapop.final[best.idx, names(est_par_vec)]

  .first_tok <- sapply(strsplit(names(est_par_vec), "_"), `[`, 1)
  .tok2      <- sapply(strsplit(names(est_par_vec), "_"), `[`, 2)
  .tok3      <- sapply(strsplit(names(est_par_vec), "_"), `[`, 3)

  bestfit_pars <- data.frame(
    par_label  = names(est_par_vec),
    pred       = ifelse(.first_tok %in% c('vul','disp'), as.numeric(.tok2), NA),
    prey       = ifelse(.first_tok == 'vul',              as.numeric(.tok3), NA),
    envresp    = ifelse(grepl('^env[0-9]+$', .first_tok), as.numeric(sub('env','', .first_tok)),
                  ifelse(grepl('^rt[0-9]+$',  .first_tok), as.numeric(sub('rt','',  .first_tok)), NA)),
    fleet      = ifelse(.first_tok == 'flt', as.numeric(.tok2), NA),
    base_pars  = est_par_vec,
    final_pars = as.numeric(unlist(ga.bestfit)),
    stringsAsFactors = FALSE
  )

  # vuls
  if(fit.vuls && !is.null(vuls.base)){
    vuls.final <- vuls.base
    vul.rows   <- which(!is.na(bestfit_pars$pred) & !is.na(bestfit_pars$prey) &
                        grepl("^vul_", bestfit_pars$par_label))
    for(i in vul.rows){
      pd <- bestfit_pars$pred[i]; py <- bestfit_pars$prey[i]
      vuls.final[py, pd + 1] <- round(bestfit_pars$final_pars[i], 2)
    }
    write.csv(vuls.final, file.path(gui_dir, paste0("vuls_ga_final_", timestamp, ".csv")))
  }

  # dispersal
  if(fit.disp && !is.null(disp.base)){
    disp.final <- disp.base
    disp.final$disp_final <- disp.final$disp
    disp_rows <- substr(bestfit_pars$par_label, 1, 4) == 'disp'
    disp.final$disp_final[bestfit_pars$pred[disp_rows]] <- bestfit_pars$final_pars[disp_rows]
    write.csv(disp.final, file.path(gui_dir, paste0("disp_ga_final_", timestamp, ".csv")))
  }

  # fleetdyn
  if(fit.fltdyn && !is.null(fleetdyn.base)){
    fleetdyn.final <- fleetdyn.base
    fleetdyn.final$effective.power_final <- fleetdyn.final$effective.power
    fleetdyn.final$effort.mult_final     <- fleetdyn.final$effort.mult
    pow_rows  <- grep("pow",  bestfit_pars$par_label)
    mult_rows <- grep("mult", bestfit_pars$par_label)
    fleetdyn.final$effective.power_final[bestfit_pars$fleet[pow_rows]]  <- bestfit_pars$final_pars[pow_rows]
    fleetdyn.final$effort.mult_final    [bestfit_pars$fleet[mult_rows]] <- bestfit_pars$final_pars[mult_rows]
    write.csv(fleetdyn.final, file.path(gui_dir, paste0("fleetdyn_ga_final_", timestamp, ".csv")))
  }

  # env + redtide (parsed from best.cmd; already absolute)
  if((fit.env || fit.redtide) && file.exists(best.cmd)){
    cmd_lines <- as.character(read.delim(best.cmd, header = FALSE,
                                         blank.lines.skip = TRUE)[, 1])
    tags_env  <- cmd_lines[grep("<ECOSPACE_ENVIRONMENTAL_RESPONSE", cmd_lines)]
    tags_env  <- tags_env[substr(tags_env, 1, 1) == "<"]
    pars_env_final <- do.call(rbind, lapply(tags_env, function(x){
      paren_num   <- as.numeric(gsub(".*\\((\\d+)\\).*", "\\1", x))
      middle_part <- strsplit(x, ",")[[1]][2]
      numbers     <- as.numeric(strsplit(trimws(middle_part), "\\s+")[[1]])
      c(paren_num, numbers)
    }))
    pars_env_final <- t(pars_env_final)
    rownames(pars_env_final) <- c('resp_num', 'shape_type',
                                  paste0('par', 1:(nrow(pars_env_final) - 2)))
    write.csv(pars_env_final,
              file.path(gui_dir, paste0("env_ga_final_", timestamp, ".csv")))
  }
  cat(sprintf("Wrote %d GUI-ready CSVs.\n", length(list.files(gui_dir))))

  # ---- Section 4. README ----
  base_ll <- if(file.exists(ga_csv)){
    tryCatch({r <- read.csv(ga_csv, header = FALSE, skip = 6)
              sprintf("%.2f (gen=-1)", r[which(r[, 1] == -1)[1], 2])},
             error = function(e) "unavailable")
  } else "ga_results CSV missing"

  cfg <- function(nm, fmt = "%s")
    if(!is.null(myconfig) && !is.null(myconfig[[nm]]))
      sprintf(fmt, myconfig[[nm]]) else "?"

  readme <- c(
    sprintf("# Phase %d GA calibration archive", phase),
    "",
    sprintf("**GA timestamp:** `%s`", timestamp),
    sprintf("**Archived:**     `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Config",
    sprintf("- Population size:   %s", cfg("popSize")),
    sprintf("- Generations:       %s", cfg("n_gen")),
    sprintf("- Elitism:           %s", cfg("elitism")),
    sprintf("- Selection:         %s (k=%s)", cfg("selection.method"), cfg("tournament.k")),
    sprintf("- Objective fxn:     %s (fn.ecospace_objfxn)", cfg("obj.fxn")),
    sprintf("- Spatial LL:        %s (weight=%s)", cfg("fit.spatial"), cfg("spatial.weight")),
    sprintf("- fit.abs.catch:     %s", cfg("fit.abs.catch")),
    "",
    "## Results",
    sprintf("- Best LL:           %.2f (gapop row %d, role: %s)",
            fitness.final[best.idx], best.idx, best.role),
    sprintf("- Base LL:           %s", base_ll),
    sprintf("- Best-fit cmd file: %s", best.cmd),
    "",
    "## Layout",
    "- `code/`      -- calibration script snapshots + R4EwE git rev",
    if(isTRUE(move_results))
      "- `results/`   -- moved from dir.results.ga: final_ga_pop, final_ga_runs, ga_results, init_ga_pop, and any timestamp-matched CSVs (e.g. mrt_summary_*)"
    else
      "- `results/`   -- (empty; source CSVs left in `dir.results.ga`)",
    if(isTRUE(move_results))
      "- `plots/`     -- moved from dir.results.ga: timestamp-matched PDFs/PNGs (e.g. final_ga_pop_distributions)"
    else
      "- `plots/`     -- (empty; source PDFs left in `dir.results.ga`)",
    "- `gui_ready/` -- best-fit vuls/disp/fleetdyn/env CSVs (absolute values)",
    "- `logs/`      -- (empty)",
    "",
    "## Reproducing the best-fit in the GUI",
    "Use the CSVs in `gui_ready/`. Env/RT taglines were extracted from the",
    "best-fit cmd.txt, so values are already `base x multiplier` (absolute).",
    "Do NOT apply another transformation.",
    "",
    "## Notes",
    "- Ecospace has stochastic movement (individual-based model). A single",
    "  evaluation of a parvec is one draw from a noisy fitness distribution;",
    "  observed variance for a given parvec can exceed several hundred LL",
    "  units. Elite-carry preserves the ORIGINAL fitness score, which will not",
    "  reproduce exactly on a fresh rerun of the same cmd file.",
    "- For ensemble analysis (Mrt, biomass trajectories), prefer the surviving",
    "  final-gen offspring dirs (tagged `run_g<N>_iNNNN_*`, where N = final",
    "  generation) over the elite rerun dirs (tagged `run_g998_iNNNN_*`).",
    "  Offspring output is what actually earned each individual its fitness;",
    "  elite reruns are fresh stochastic draws off the same parvec and will",
    "  not reproduce the original score."
  )
  if(!is.null(notes)){
    readme <- c(readme, "", "## Additional notes",
                paste0("- ", notes))
  }
  writeLines(readme, file.path(archive_dir, "README.md"))
  cat("README.md written.\n")

  # ---- Section 5. Move timestamp-matched results out of dir.results.ga ----
  # Skip anything already inside the archive dir (defensive: overwrite=TRUE reruns).
  if(isTRUE(move_results)){
    all_ts <- list.files(dir.results.ga, pattern = timestamp,
                         full.names = TRUE, recursive = FALSE)
    all_ts <- all_ts[!startsWith(normalizePath(all_ts, mustWork = FALSE),
                                 normalizePath(archive_dir, mustWork = FALSE))]
    csv_src <- all_ts[grepl("\\.csv$",         all_ts, ignore.case = TRUE)]
    img_src <- all_ts[grepl("\\.(pdf|png)$",   all_ts, ignore.case = TRUE)]

    move_one <- function(src, dst_dir){
      dst <- file.path(dst_dir, basename(src))
      if(!file.rename(src, dst)){
        # Cross-device or locked: copy + remove
        if(file.copy(src, dst, overwrite = TRUE)) file.remove(src) else FALSE
      } else TRUE
    }
    n_csv <- sum(vapply(csv_src, move_one, logical(1),
                        dst_dir = file.path(archive_dir, "results")))
    n_img <- sum(vapply(img_src, move_one, logical(1),
                        dst_dir = file.path(archive_dir, "plots")))
    cat(sprintf("Moved %d CSV(s) to results/ and %d plot(s) to plots/.\n",
                n_csv, n_img))
  }

  cat("\n== archive complete ==\n")
  cat(sprintf("Archive: %s\n", archive_dir))
  invisible(archive_dir)
}

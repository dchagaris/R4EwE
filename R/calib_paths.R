#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calibration run paths ------------------------------------------------------------------------------
# Every calibration run (sensitivity or GA) is identified by a phase and a
# user-chosen label, never by a date. Three roots hold its files:
#   scratch  <root.scratch>/<phase>/<label>   Ecospace run dirs (big, fast disk)
#   archive  <root.archive>/<phase>/<label>   same tree after the manual move
#   results  <root.results>/<phase>/<label>   LL tables, plots, run_info.csv
# Dates and provenance live in results/run_info.csv (field,value rows).

#fn.calib_paths -------------------------------------------------------------------------------------
#' @title Paths for a calibration run
#' @description Build the scratch, archive and results paths for a calibration
#'   run identified by \code{phase} and \code{label}. Pure: creates and writes
#'   nothing, so it is safe to call when only post-processing.
#' @param phase Phase name, e.g. \code{"phase5"}.
#' @param label Run label, e.g. \code{"sens_init"} or \code{"ga01"}. Letters,
#'   digits, \code{_}, \code{-} and \code{.} only.
#' @param root.scratch Root for the Ecospace run dirs.
#' @param root.archive Root the scratch tree is moved to when the run is done.
#'   \code{NULL} if there is no archive.
#' @param root.results Root for results (LL tables, plots, \code{run_info.csv}).
#' @return A list with \code{phase}, \code{label}, \code{dir.scratch},
#'   \code{dir.archive}, \code{dir.results}, \code{file.info}, and
#'   \code{dir.run}: the scratch dir if it exists, else the archive dir if it
#'   exists, else the scratch dir.
#' @examples
#' \dontrun{
#' p <- fn.calib_paths("phase5", "ga01")
#' p$dir.run
#' }
#' @export
fn.calib_paths <- function(phase, label,
                           root.scratch = "C:/WFS MICE",
                           root.archive = "E:/WFS MICE",
                           root.results = file.path(getwd(), "ga_results")){
  for(x in list(phase = phase, label = label)){
    if(length(x) != 1 || !nzchar(x) || !grepl("^[A-Za-z0-9_.-]+$", x))
      stop("phase and label must be a single name of letters, digits, '_', '-' or '.': ", x)
  }
  dir.scratch <- file.path(root.scratch, phase, label)
  dir.archive <- if(is.null(root.archive)) NA_character_ else file.path(root.archive, phase, label)
  dir.results <- file.path(root.results, phase, label)

  dir.run <- if(dir.exists(dir.scratch)) dir.scratch
             else if(!is.na(dir.archive) && dir.exists(dir.archive)) dir.archive
             else dir.scratch

  list(phase       = phase,
       label       = label,
       dir.scratch = dir.scratch,
       dir.archive = dir.archive,
       dir.results = dir.results,
       dir.run     = dir.run,
       file.info   = file.path(dir.results, "run_info.csv"))
}

#fn.calib_info -------------------------------------------------------------------------------------
#' @title Read a calibration run's run_info.csv
#' @param paths A list from \code{\link{fn.calib_paths}}.
#' @return A named character vector (field = value), or \code{NULL} if the run
#'   has never been started.
#' @export
fn.calib_info <- function(paths){
  if(!file.exists(paths$file.info)) return(NULL)
  info <- utils::read.csv(paths$file.info, colClasses = "character")
  stats::setNames(info$value, info$field)
}

# Internal: merge fields into run_info.csv, keeping existing ones.
.calib_info_write <- function(paths, fields){
  info <- fn.calib_info(paths)
  if(is.null(info)) info <- character()
  fields <- vapply(fields, function(v) paste(as.character(v), collapse = " "), "")
  info[names(fields)] <- fields
  dir.create(paths$dir.results, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(field = names(info), value = unname(info)),
                   paths$file.info, row.names = FALSE)
  invisible(info)
}

# Internal: R4EwE commit -- the GitHub SHA when installed with remotes, else
# git HEAD of the source tree when loaded with devtools::load_all().
.r4ewe_commit <- function(){
  sha <- utils::packageDescription("R4EwE")$RemoteSha
  if(!is.null(sha)) return(substr(sha, 1, 7))
  src <- getNamespaceInfo(asNamespace("R4EwE"), "path")
  out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(src), "rev-parse", "--short", "HEAD"),
                                           stdout = TRUE, stderr = FALSE)),
                  error = function(e) character())
  if(length(out) == 1) out else NA_character_
}

# Internal: value of a <TAG> line in an EwE console command file.
.cmd_value <- function(cmd, tag){
  hit <- grep(paste0("^<", tag, ">"), cmd, value = TRUE)
  if(!length(hit)) return(NA_character_)
  trimws(strsplit(hit[1], ",")[[1]][2])
}

#fn.calib_start ------------------------------------------------------------------------------------
#' @title Start a calibration run
#' @description Call once, immediately before launching a run. Creates the
#'   scratch and results dirs and records the start in \code{run_info.csv}
#'   (status \code{running}). Refuses a label that is already \code{done}
#'   unless \code{overwrite = TRUE}, which starts a fresh \code{run_info.csv};
#'   warns on one still marked \code{running}.
#'   Does not delete anything: clearing old run dirs is up to the caller.
#' @param paths A list from \code{\link{fn.calib_paths}}.
#' @param cmd_file Optional EwE console command file; its model file and
#'   Ecosim/Ecospace scenario indices are recorded.
#' @param overwrite Allow restarting a label whose status is \code{done}.
#' @param ... Extra fields to record, e.g. \code{note = "spin-up 1 yr"}.
#' @return Invisibly, the run_info fields.
#' @export
fn.calib_start <- function(paths, cmd_file = NULL, overwrite = FALSE, ...){
  info <- fn.calib_info(paths)
  status <- if(is.null(info)) NA else info[["status"]]
  if(identical(status, "done")){
    if(!overwrite)
      stop("Run '", paths$phase, "/", paths$label, "' is already done (finished ",
           info[["finished"]], "). Pick a new label or set overwrite = TRUE.")
    unlink(paths$file.info)   # a restarted run must not inherit the old run's fields
  }
  if(identical(status, "running"))
    warning("Run '", paths$phase, "/", paths$label, "' is marked running (started ",
            info[["started"]], "); starting it again.")
  if(!is.na(paths$dir.archive) && dir.exists(paths$dir.archive))
    warning("Archive dir already exists and may hold an earlier run: ", paths$dir.archive)

  dir.create(paths$dir.scratch, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$dir.results, recursive = TRUE, showWarnings = FALSE)

  fields <- list(phase    = paths$phase,
                 label    = paths$label,
                 status   = "running",
                 started  = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 finished = "",
                 dir.scratch = paths$dir.scratch,
                 dir.archive = paths$dir.archive,
                 r4ewe_commit = .r4ewe_commit())
  if(!is.null(cmd_file)){
    cmd <- readLines(cmd_file, warn = FALSE)
    fields$cmd_file         <- normalizePath(cmd_file, winslash = "/")
    fields$model_file       <- basename(.cmd_value(cmd, "EWE_MODEL_FILE"))
    fields$ecosim_scenario  <- .cmd_value(cmd, "ECOSIM_SCENARIO_INDEX")
    fields$ecospace_scenario <- .cmd_value(cmd, "ECOSPACE_SCENARIO_INDEX")
  }
  fields <- c(fields, list(...))
  invisible(.calib_info_write(paths, fields))
}

#fn.calib_finish -----------------------------------------------------------------------------------
#' @title Finish a calibration run
#' @description Record the end of a run in \code{run_info.csv}.
#' @param paths A list from \code{\link{fn.calib_paths}}.
#' @param status Final status, \code{"done"} (default) or e.g. \code{"stopped"}.
#' @param ... Extra fields to record, e.g. \code{best_LL = 5124.98}.
#' @return Invisibly, the run_info fields.
#' @export
fn.calib_finish <- function(paths, status = "done", ...){
  if(is.null(fn.calib_info(paths)))
    stop("No run_info.csv for '", paths$phase, "/", paths$label, "'; call fn.calib_start first.")
  invisible(.calib_info_write(paths, c(list(status   = status,
                                            finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                                       list(...))))
}

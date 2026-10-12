#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calibration phases and runs ------------------------------------------------------------------------
# A phase owns its input files and settings; each run within it is identified
# by a user-chosen label, never by a date.
#   <root.results>/<phase>/inputs/           the phase's input files (edited / exported)
#   <root.results>/<phase>/<label>/          LL tables, plots, run_info.csv
#   <root.results>/<phase>/<label>/inputs/   frozen copy taken by fn.calib_start
#   <root.scratch>/<phase>/<label>           Ecospace run dirs (big, fast disk)
#   <root.archive>/<phase>/<label>           same tree after the manual move
# Dates, settings, input hashes and provenance live in run_info.csv (field,value rows).

#fn.calib_phase -------------------------------------------------------------------------------------
#' @title Define a calibration phase
#' @description Bundle a phase's roots, input files, external references and
#'   settings into one object that \code{\link{fn.calib_paths}} builds runs
#'   from. Pure apart from checking that the input files exist.
#' @param phase Phase name, e.g. \code{"phase5"}.
#' @param inputs Named character vector of input file names, resolved against
#'   \code{<root.results>/<phase>/inputs/}. Copied into each run at launch.
#' @param refs Named character vector of shared external files or dirs (obs
#'   maps, console exe, ...). Paths are recorded, nothing is copied.
#' @param settings Named list of run settings (years, objective-function
#'   options, ...). Recorded in each run's \code{run_info.csv}.
#' @param root.scratch Root for the Ecospace run dirs.
#' @param root.archive Root the scratch tree is moved to when a run is done,
#'   or \code{NULL}.
#' @param root.results Root for the phase folder.
#' @return A \code{calib_phase} list with \code{phase}, the three roots,
#'   \code{dir.phase}, \code{dir.inputs}, \code{files} (full paths),
#'   \code{refs} and \code{settings}.
#' @examples
#' \dontrun{
#' ph <- fn.calib_phase("phase5",
#'         inputs   = c(cmdbase = "WFS MICE cmd 5min phase5.txt"),
#'         settings = list(styear = 1985, enyear = 2024))
#' ph$files$cmdbase
#' }
#' @export
fn.calib_phase <- function(phase, inputs, refs = NULL, settings = list(),
                           root.scratch = "C:/WFS MICE",
                           root.archive = "E:/WFS MICE",
                           root.results = file.path(getwd(), "ga_results")){
  .check_name(phase)
  if(is.null(names(inputs)) || any(!nzchar(names(inputs))) || anyDuplicated(names(inputs)))
    stop("inputs must be a named vector with unique names, e.g. c(cmdbase = 'cmd.txt').")
  dir.phase  <- file.path(root.results, phase)
  dir.inputs <- file.path(dir.phase, "inputs")
  files <- as.list(stats::setNames(file.path(dir.inputs, inputs), names(inputs)))
  missing <- !file.exists(unlist(files))
  if(any(missing))
    stop("Input files missing from ", dir.inputs, ":\n  ",
         paste(inputs[missing], collapse = "\n  "))
  if(!is.null(refs)){
    refs <- as.list(refs)
    if(any(!file.exists(unlist(refs))))
      warning("Reference paths not found:\n  ",
              paste(unlist(refs)[!file.exists(unlist(refs))], collapse = "\n  "))
  }
  structure(list(phase        = phase,
                 root.scratch = root.scratch,
                 root.archive = root.archive,
                 root.results = root.results,
                 dir.phase    = dir.phase,
                 dir.inputs   = dir.inputs,
                 files        = files,
                 refs         = refs,
                 settings     = settings),
            class = "calib_phase")
}

#fn.calib_paths -------------------------------------------------------------------------------------
#' @title Paths for a calibration run
#' @description Build the scratch, archive and results paths for run
#'   \code{label} of a phase. Pure: creates and writes nothing, so it is safe
#'   to call when only post-processing.
#' @param phase A \code{calib_phase} from \code{\link{fn.calib_phase}}, or a
#'   phase name (then the \code{root.*} arguments are used and the run has no
#'   input files or settings).
#' @param label Run label, e.g. \code{"sens_init"} or \code{"ga01"}. Letters,
#'   digits, \code{_}, \code{-} and \code{.} only.
#' @param root.scratch,root.archive,root.results Roots, used only when
#'   \code{phase} is a name.
#' @return A list with \code{phase}, \code{label}, \code{dir.scratch},
#'   \code{dir.archive}, \code{dir.results}, \code{file.info},
#'   \code{dir.inputs} (the run's frozen copy), \code{files}, \code{refs},
#'   \code{settings}, and \code{dir.run}: the scratch dir if it holds
#'   anything, else the archive dir if it exists, else the scratch dir.
#'   \code{files} points at the run's frozen copies once the run has been
#'   started, else at the phase inputs.
#' @examples
#' \dontrun{
#' p <- fn.calib_paths(ph, "ga01")
#' p$dir.run
#' }
#' @export
fn.calib_paths <- function(phase, label,
                           root.scratch = "C:/WFS MICE",
                           root.archive = "E:/WFS MICE",
                           root.results = file.path(getwd(), "ga_results")){
  ph <- if(inherits(phase, "calib_phase")) phase
        else list(phase = phase, root.scratch = root.scratch, root.archive = root.archive,
                  root.results = root.results, files = list(), refs = NULL, settings = list())
  .check_name(ph$phase); .check_name(label)

  dir.scratch <- file.path(ph$root.scratch, ph$phase, label)
  dir.archive <- if(is.null(ph$root.archive)) NA_character_ else file.path(ph$root.archive, ph$phase, label)
  dir.results <- file.path(ph$root.results, ph$phase, label)
  dir.inputs  <- file.path(dir.results, "inputs")

  # an empty scratch dir (left behind by robocopy /MOVE) does not count
  has.files <- function(d) dir.exists(d) && length(list.files(d, all.files = TRUE, no.. = TRUE)) > 0
  dir.run <- if(has.files(dir.scratch)) dir.scratch
             else if(!is.na(dir.archive) && dir.exists(dir.archive)) dir.archive
             else dir.scratch

  files.phase <- ph$files
  files.run   <- lapply(files.phase, function(f) file.path(dir.inputs, basename(f)))
  started     <- length(files.run) && all(file.exists(unlist(files.run)))

  list(phase       = ph$phase,
       label       = label,
       dir.scratch = dir.scratch,
       dir.archive = dir.archive,
       dir.results = dir.results,
       dir.run     = dir.run,
       dir.inputs  = dir.inputs,
       file.info   = file.path(dir.results, "run_info.csv"),
       files       = if(started) files.run else files.phase,
       files.phase = files.phase,
       refs        = ph$refs,
       settings    = ph$settings)
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

# Internal: phase / label names must be safe as single folder names.
.check_name <- function(x){
  if(length(x) != 1 || !nzchar(x) || !grepl("^[A-Za-z0-9_.-]+$", x))
    stop("phase and label must be a single name of letters, digits, '_', '-' or '.': ", x)
  invisible(x)
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

# Internal: run_info fields describing an external file (path, size, mtime, md5).
.file_fields <- function(prefix, path){
  if(is.na(path) || !file.exists(path))
    return(stats::setNames(list(path), prefix))
  fi <- file.info(path)
  stats::setNames(list(normalizePath(path, winslash = "/"), fi$size,
                       format(fi$mtime, "%Y-%m-%d %H:%M:%S"), unname(tools::md5sum(path))),
                  paste0(prefix, c("", ".size", ".mtime", ".md5")))
}

#fn.calib_start ------------------------------------------------------------------------------------
#' @title Start a calibration run
#' @description Call once, immediately before launching a run. Creates the
#'   scratch and results dirs, copies the phase input files into the run's
#'   \code{inputs/} folder, and records the start in \code{run_info.csv}
#'   (status \code{running}) with the settings, input md5s, R4EwE commit and,
#'   from the cmd file, the model, spatial-temporal config and time-series
#'   files (path, size, mtime, md5). Refuses a label that is already
#'   \code{done} unless \code{overwrite = TRUE}, which starts a fresh
#'   \code{run_info.csv} and fresh input copies. On a label still marked
#'   \code{running} it warns and keeps the original input copies, warning if
#'   the phase inputs have changed since. Does not delete run dirs: clearing
#'   old output is up to the caller.
#' @param paths A list from \code{\link{fn.calib_paths}}.
#' @param cmd_file EwE console command file. Defaults to \code{files$cmdbase}.
#' @param overwrite Allow restarting a label whose status is \code{done}.
#' @param ... Extra fields to record, e.g. \code{note = "spin-up 1 yr"}.
#' @return Invisibly, \code{paths} with \code{files} pointing at the run's
#'   frozen input copies.
#' @export
fn.calib_start <- function(paths, cmd_file = paths$files.phase$cmdbase, overwrite = FALSE, ...){
  info <- fn.calib_info(paths)
  status <- if(is.null(info)) NA else info[["status"]]
  run <- paste0("'", paths$phase, "/", paths$label, "'")
  fresh <- TRUE
  if(identical(status, "done")){
    if(!overwrite)
      stop("Run ", run, " is already done (finished ", info[["finished"]],
           "). Pick a new label or set overwrite = TRUE.")
    unlink(c(paths$file.info, paths$dir.inputs), recursive = TRUE)  # no inheritance from the old run
  } else if(identical(status, "running")){
    warning("Run ", run, " is marked running (started ", info[["started"]], "); starting it again.")
    fresh <- FALSE
  }
  if(!is.na(paths$dir.archive) && dir.exists(paths$dir.archive))
    warning("Archive dir already exists and may hold an earlier run: ", paths$dir.archive)

  dir.create(paths$dir.scratch, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$dir.inputs,  recursive = TRUE, showWarnings = FALSE)

  # freeze the phase inputs into the run
  src <- unlist(paths$files.phase)
  dst <- file.path(paths$dir.inputs, basename(src))
  if(length(src)){
    keep <- !fresh & file.exists(dst)
    if(any(keep)){
      changed <- unname(tools::md5sum(src[keep])) != unname(tools::md5sum(dst[keep]))
      if(any(changed))
        warning("Phase inputs changed since this run started; keeping the run's original copies of:\n  ",
                paste(basename(src[keep][changed]), collapse = "\n  "))
    }
    ok <- file.copy(src[!keep], dst[!keep], overwrite = TRUE, copy.date = TRUE)
    if(!all(ok)) stop("Failed to copy inputs: ", paste(basename(src[!keep][!ok]), collapse = ", "))
  }
  paths$files <- as.list(stats::setNames(dst, names(paths$files.phase)))

  fields <- list(phase    = paths$phase,
                 label    = paths$label,
                 status   = "running",
                 started  = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 finished = "",
                 dir.scratch  = paths$dir.scratch,
                 dir.archive  = paths$dir.archive,
                 r4ewe_commit = .r4ewe_commit())
  if(length(paths$settings))
    fields <- c(fields, stats::setNames(paths$settings, paste0("setting.", names(paths$settings))))
  for(nm in names(paths$files)){
    fields[[paste0("input.", nm)]]        <- basename(paths$files[[nm]])
    fields[[paste0("input.", nm, ".md5")]] <- unname(tools::md5sum(paths$files[[nm]]))
  }
  if(length(paths$refs))
    fields <- c(fields, stats::setNames(paths$refs, paste0("ref.", names(paths$refs))))
  if(!is.null(cmd_file)){
    cmd <- readLines(cmd_file, warn = FALSE)
    fields$cmd_file          <- normalizePath(cmd_file, winslash = "/")
    fields$ecosim_scenario   <- .cmd_value(cmd, "ECOSIM_SCENARIO_INDEX")
    fields$ecospace_scenario <- .cmd_value(cmd, "ECOSPACE_SCENARIO_INDEX")
    fields <- c(fields,
                .file_fields("model_file",     .cmd_value(cmd, "EWE_MODEL_FILE")),
                .file_fields("st_config_file", .cmd_value(cmd, "SPATIAL_CONFIG_FILE")),
                .file_fields("ts_file",        .cmd_value(cmd, "TIME_SERIES_FILE")))
  }
  fields <- c(fields, list(...))
  .calib_info_write(paths, fields)
  invisible(paths)
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

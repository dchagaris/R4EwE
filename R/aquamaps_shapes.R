#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AquaMaps environmental-envelope response shapes ------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check that a compiler toolchain (Rtools on Windows) is available
#'
#' Installing \pkg{aquamapsdata} from GitHub builds it from source, which needs a
#' working compiler toolchain (Rtools on Windows). Returns invisibly on success,
#' otherwise stops with platform-specific guidance.
#' @return Invisibly \code{TRUE} if build tools are found; otherwise an error.
#' @keywords internal
#' @noRd
fn.check_build_tools <- function() {
  if (!requireNamespace("pkgbuild", quietly = TRUE)) utils::install.packages("pkgbuild")
  if (pkgbuild::has_build_tools(debug = FALSE)) return(invisible(TRUE))

  if (.Platform$OS.type == "windows") {
    stop("Rtools was not found, but it is required to build aquamapsdata from source.\n",
         "Install the version matching your R (", getRversion(), ") from:\n",
         "  https://cran.r-project.org/bin/windows/Rtools/\n",
         "then restart R and try again.", call. = FALSE)
  }
  stop("A C/C++ compiler toolchain was not found, but it is required to build ",
       "aquamapsdata from source. Install the build tools for your OS and retry.",
       call. = FALSE)
}

#' Ensure the packages needed for the AquaMaps workflow are installed
#'
#' Checks for \pkg{aquamapsdata} (GitHub-only) and \pkg{dplyr} (CRAN). If either
#' is missing it prompts the user to install them in an interactive session, or
#' stops with copy-paste install instructions in a non-interactive session.
#' Because \pkg{aquamapsdata} is built from source, Rtools is verified (via
#' \code{fn.check_build_tools()}) before that install runs.
#' @return Invisibly \code{TRUE} when all dependencies are available.
#' @keywords internal
#' @noRd
fn.ensure_aquamaps_deps <- function() {
  deps <- c("aquamapsdata", "dplyr")
  missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))

  msg <- paste0("Required package(s) not installed: ", paste(missing, collapse = ", "), ".")
  instructions <- paste(
    "  # Windows: first install Rtools from https://cran.r-project.org/bin/windows/Rtools/",
    "  install.packages('remotes')",
    "  remotes::install_github('raquamaps/aquamapsdata', dependencies = TRUE)",
    "  install.packages('dplyr')",
    sep = "\n")

  if (!interactive()) {
    stop(msg, "\nInstall them with:\n", instructions, call. = FALSE)
  }

  message(msg)
  ans <- utils::menu(c("Yes, install now", "No, abort"),
                     title = "Install the required AquaMaps package(s) now?")
  if (ans != 1L) stop("Required packages not installed; aborting.", call. = FALSE)

  if (!requireNamespace("remotes", quietly = TRUE)) utils::install.packages("remotes")
  for (pkg in missing) {
    if (pkg == "aquamapsdata") {
      fn.check_build_tools()  # aquamapsdata builds from source -> needs Rtools
      remotes::install_github("raquamaps/aquamapsdata", dependencies = TRUE)
    } else {
      utils::install.packages(pkg)
    }
  }

  # confirm the install actually succeeded before continuing
  still_missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing)) {
    stop("Installation did not succeed for: ", paste(still_missing, collapse = ", "),
         ". Please install manually:\n", instructions, call. = FALSE)
  }
  invisible(TRUE)
}

#' @title Connect to the AquaMaps SQLite database.
#' @description Points the \pkg{aquamapsdata} connection cache at an AquaMaps
#'   SQLite database so that downstream queries (\code{am_search_exact()},
#'   \code{am_hspen()}) read from the chosen file. The \pkg{aquamapsdata} package
#'   hard-codes its database location, so a custom \code{db_path} is supported by
#'   connecting directly and injecting the connection into the package's internal
#'   cache.
#' @param db_path Optional path to an \code{am.db} SQLite file. If \code{NULL}
#'   (default) the package's standard location is used (and an informative error
#'   is raised, pointing to \code{aquamapsdata::download_db()}, if it is missing).
#' @return Invisibly, the path of the database that was connected.
#' @details The database is ~2 GB. Use \code{aquamapsdata::download_db()} once to
#'   fetch it to the default location, or pass \code{db_path} to use a copy kept
#'   elsewhere (e.g. a shared drive) and avoid re-downloading.
#' @export
fn.connect_aquamaps_db <- function(db_path = NULL) {
  fn.ensure_aquamaps_deps()

  am_db_sqlite <- utils::getFromNamespace("am_db_sqlite", "aquamapsdata")
  default_path <- am_db_sqlite()

  if (is.null(db_path)) {
    db_path <- default_path
    if (!file.exists(db_path)) {
      stop("No AquaMaps database found at the default location:\n  ", db_path,
           "\nDownload it (~2 GB, one time) with:\n",
           "  aquamapsdata::download_db()\n",
           "or, if you already have an am.db file elsewhere, pass its path via db_path=.",
           call. = FALSE)
    }
    aquamapsdata::default_db("sqlite")
  } else {
    if (!file.exists(db_path)) stop("No file found at db_path: ", db_path, call. = FALSE)
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path,
                          synchronous = "normal", flags = RSQLite::SQLITE_RO)
    db_cache <- utils::getFromNamespace("db_cache", "aquamapsdata")
    assign("local_db", con, envir = db_cache)
    assign("source", "sqlite", envir = db_cache)
  }

  message("Using AquaMaps database: ", db_path)
  invisible(db_path)
}

#' @title Build AquaMaps trapezoid response shapes from a species list.
#' @description Reads a "species list for AquaMaps query" CSV, queries the
#'   AquaMaps HSPEN (environmental envelope preference) data, and returns a
#'   \code{shapes.out} data frame of trapezoid (EwE function type 9) response
#'   shapes for depth, temperature, salinity, distance-to-shore, and dissolved
#'   oxygen for every species found in AquaMaps. Each shape's four trapezoid
#'   parameters come from the species' \code{Min}, \code{PrefMin}, \code{PrefMax},
#'   and \code{Max} envelope values for the corresponding variable.
#' @param species_csv Path to the species list CSV. Must contain columns
#'   \code{number} (functional group number), \code{name} (functional group
#'   name), and \code{species} (scientific name). Any other columns are ignored.
#' @param db_path Optional path to an \code{am.db} SQLite file. If \code{NULL}
#'   (default) the \pkg{aquamapsdata} package's standard database location is
#'   used; pass a path to use a copy kept elsewhere (see
#'   \code{\link{fn.connect_aquamaps_db}}).
#' @return A data frame with 11 columns: \code{Function name}, \code{Function
#'   type}, \code{Param 1}-\code{Param 6}, \code{var}, \code{fg num}, and
#'   \code{source}. Parameter columns (\code{Function type}, \code{Param 1-6},
#'   \code{fg num}) are numeric and rounded to 3 decimal places.
#' @details Scientific names in the species list that are not matched in AquaMaps
#'   are reported via \code{message()} and simply omitted from the output. The
#'   first call may prompt to install \pkg{aquamapsdata}/\pkg{dplyr} (and verify
#'   Rtools) if they are not already present.
#' @seealso \code{\link{fn.connect_aquamaps_db}}
#' @examples
#' \dontrun{
#' # use the package's default database location:
#' shapes.out <- fn.make_aquamaps_shapes("species list for aquamaps query.csv")
#'
#' # use an am.db stored on a shared drive:
#' shapes.out <- fn.make_aquamaps_shapes("species list for aquamaps query.csv",
#'                                       db_path = "D:/data/aquamaps/am.db")
#' }
#' @export
fn.make_aquamaps_shapes <- function(species_csv, db_path = NULL) {

  # make sure aquamapsdata + dplyr are available (prompts to install if not)
  fn.ensure_aquamaps_deps()

  # full species list
  aqm.spp <- utils::read.csv(species_csv)
  aqm.spp$species[aqm.spp$species == ""] <- NA
  do.spp <- sort(unique(aqm.spp$species))

  # connect to aquamaps data and resolve scientific names to AquaMaps records
  fn.connect_aquamaps_db(db_path)
  spplist <- aquamapsdata::am_search_exact()
  spplist$sciname <- paste(spplist$Genus, spplist$Species)
  do.spplist <- spplist[spplist$sciname %in% do.spp, ]

  missing <- do.spp[which(!do.spp %in% do.spplist$sciname)]
  if (length(missing)) {
    message("Species not found in AquaMaps: ", paste(missing, collapse = ", "))
  }

  # get environmental envelope preference data
  hspen <- as.data.frame(aquamapsdata::am_hspen())
  hspen <- merge(do.spplist, hspen, all.x = TRUE)

  shapes.out <- data.frame(matrix(NA, nrow = 0, ncol = 11))
  names(shapes.out) <- c('Function name', 'Function type',
                         'Param 1', 'Param 2', 'Param 3',
                         'Param 4', 'Param 5', 'Param 6',
                         'var', 'fg num', 'source')

  # AquaMaps shapes
  for (i in 1:nrow(aqm.spp)) {
    fg.num  <- aqm.spp$number[i]
    fg.name <- aqm.spp$name[i]
    if (aqm.spp$species[i] %in% hspen$sciname) {
      sub <- hspen[hspen$sciname == aqm.spp$species[i], ]

      fxn.names <- paste0(fg.num, "_", fg.name, "_",
                          c('depth', 'temp', 'salinity', 'dist2shore', 'DO'), "_aqm")
      fxn.types <- rep(9, length(fxn.names))

      # Find matching columns (Min/PrefMin/PrefMax/Max for each variable)
      matching_cols <- which(
        grepl(paste0("(", paste(c("Min", "PrefMin", "PrefMax", "Max"), collapse = "|"), ")$"), names(sub)) &
          grepl(paste(c("Depth", "Temp", "Salinity", "LandDist", "Oxy"), collapse = "|"), names(sub))
      )
      pars   <- as.numeric(sub[, matching_cols])
      parmat <- cbind(matrix(pars, ncol = 4, byrow = TRUE),
                      matrix(NA, nrow = length(fxn.types), ncol = 2))

      dat.i <- as.data.frame(cbind(fxn.names, fxn.types, parmat))
      dat.i$var    <- sapply(strsplit(dat.i$fxn.names, "_"), `[`, 3)
      dat.i$fg.num <- fg.num
      dat.i$source <- 'aqm'
      names(dat.i) <- names(shapes.out)
      shapes.out <- rbind(shapes.out, dat.i)
    }
  }

  # coerce parameter columns (Function type, Param 1-6, fg num) to numeric and round
  shapes.out[, c(2:8, 10)] <- lapply(shapes.out[, c(2:8, 10)], as.numeric)
  shapes.out[, c(2:8, 10)] <- round(shapes.out[, c(2:8, 10)], 3)

  shapes.out
}

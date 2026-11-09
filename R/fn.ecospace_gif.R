#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecospace output animations-----------------------------------------------------
# fn.ecospace_gif() animates Ecospace spatial output. It reads either the per-group
# map csv files (EcospaceMap<Variable>-<group>.csv) or the per-timestep ascii grids
# in an asc/ folder, and writes one gif per variable/group series.

#' @title Parse the metadata header of an Ecospace map csv.
#' @description Reads the key/value block above \code{<HEADER end/>} in an
#'   \code{EcospaceMap<Variable>-<group>.csv} file and returns it as a named list,
#'   converting the map geometry and timing fields to numeric.
#' @param path Path to an EcospaceMap csv file.
#' @return Named list of header fields.
#' @keywords internal
#' @noRd
.ecospace_map_header <- function(path){
  ln      <- readLines(path, n = 200L, warn = FALSE)
  end_idx <- grep("<HEADER end/>", ln, fixed = TRUE)[1]
  if(is.na(end_idx)) stop("Could not locate <HEADER end/> in ", path, call. = FALSE)
  ln <- ln[seq_len(end_idx - 1L)]
  ln <- ln[grepl(",", ln, fixed = TRUE) & !grepl("^\"?<HEADER", ln)]

  strip <- function(x) gsub('^"|"$', "", trimws(x))
  key <- strip(sub(",.*$", "", ln))
  val <- strip(sub("^[^,]*,", "", ln))
  out <- as.list(val)
  names(out) <- key

  for(k in c("MapRows", "MapCols", "MapCellSize", "MapCellLength",
             "MapTopLeftLat", "MapTopLeftLon", "StartYear",
             "EcoSpaceTimeStepLength")){
    if(!is.null(out[[k]])) out[[k]] <- suppressWarnings(as.numeric(out[[k]]))
  }
  out
}

#' @title Read an Ecospace map csv into a georeferenced raster stack.
#' @description Wraps \code{.parse_ecospace_map_csv()}, taking the grid dimensions and
#'   georeferencing from the file's own header rather than requiring them as arguments.
#'   NODATA cells (\code{-9999}) become \code{NA}. Layer order follows the Step blocks
#'   in the file.
#' @param path Path to an EcospaceMap csv file.
#' @param model_styear Calendar year corresponding to model year 0. \code{NULL} (default)
#'   reads \code{StartYear} from the header.
#' @return A \code{terra} SpatRaster with one layer per timestep, carrying \code{step},
#'   \code{year} (model-relative), \code{time} (calendar) and \code{header} attributes.
#' @keywords internal
#' @noRd
.ecospace_map_stack <- function(path, model_styear = NULL){
  hdr <- .ecospace_map_header(path)
  nr  <- hdr$MapRows
  nc  <- hdr$MapCols
  if(is.null(nr) || is.null(nc) || is.na(nr) || is.na(nc))
    stop("MapRows / MapCols missing from header of ", path, call. = FALSE)

  arr <- .parse_ecospace_map_csv(path, nrows = as.integer(nr), ncols = as.integer(nc))
  step <- attr(arr, "step")
  yr   <- attr(arr, "year")
  arr[arr == -9999] <- NA_real_

  r <- terra::rast(arr)

  # extent from the top-left corner + cell size, when the header provides them
  cs <- hdr$MapCellSize
  if(!is.null(cs) && !is.na(cs) && !is.null(hdr$MapTopLeftLon) && !is.null(hdr$MapTopLeftLat)){
    xmin <- hdr$MapTopLeftLon
    ymax <- hdr$MapTopLeftLat
    terra::ext(r) <- c(xmin, xmin + nc * cs, ymax - nr * cs, ymax)
    # EwE writes the WKT with single quotes; terra needs double
    wkt <- hdr$CoordinateSystemWKT
    if(!is.null(wkt) && nzchar(wkt)){
      ok <- tryCatch({terra::crs(r) <- gsub("'", '"', wkt, fixed = TRUE); TRUE},
                     error = function(e) FALSE, warning = function(w) FALSE)
      if(!ok) try(terra::crs(r) <- "EPSG:4326", silent = TRUE)
    }
  }

  if(is.null(model_styear)) model_styear <- hdr$StartYear
  time <- if(!is.null(model_styear) && !is.na(model_styear)) model_styear + yr else yr
  names(r) <- sprintf("step_%s", step)
  attr(r, "step")   <- step
  attr(r, "year")   <- yr
  attr(r, "time")   <- time
  attr(r, "header") <- hdr
  r
}

#' @title Index the Ecospace map csv files in a run folder.
#' @description GUI runs put the per-group map csvs in a \code{csv/} subfolder; console
#'   runs write them at the top level of the run folder. Both are handled.
#' @param dir.out Ecospace run output folder.
#' @return Data frame with columns \code{file}, \code{variable}, \code{group}, or
#'   \code{NULL} if no map csvs are present.
#' @keywords internal
#' @noRd
.ecospace_map_files <- function(dir.out){
  csv_dir <- if(dir.exists(file.path(dir.out, "csv"))) file.path(dir.out, "csv") else dir.out
  f <- list.files(csv_dir, pattern = "^EcospaceMap.*\\.csv$", full.names = TRUE)
  if(length(f) == 0) return(NULL)
  b <- basename(f)
  data.frame(file     = f,
             variable = sub("^EcospaceMap([^-]*)-.*$", "\\1", b),
             group    = sub("\\.csv$", "", sub("^EcospaceMap[^-]*-", "", b)),
             stringsAsFactors = FALSE)
}

#' @title Format a calendar time vector as YYYY-MM frame labels.
#' @param time Numeric vector of fractional calendar years.
#' @return Character vector of \code{YYYY-MM} labels.
#' @keywords internal
#' @noRd
.ecospace_datestamp <- function(time){
  yr <- floor(time + 1e-8)
  mo <- round((time - yr) * 12) + 1
  yr[mo > 12] <- yr[mo > 12] + 1
  mo[mo > 12] <- 1
  sprintf("%d-%02d", yr, mo)
}

#' @title Read a spatial-temporal driver folder into a dated raster stack.
#' @description Reads the monthly ascii grids that make up one ST driver variable --
#'   Ecospace \emph{input}, not output. Files are dated from their own name suffix:
#'   \code{<var>_YYYYMM.asc} for real months, and \code{<var>_YYYY_MM.asc} for the
#'   12-month climatology written with a placeholder year (e.g. \code{_0000_01}).
#'   Only the folder itself is scanned, so a \code{static/} subfolder of long-term
#'   means is left alone.
#' @param dir.stdriver Folder holding one driver variable's ascii grids.
#' @param datestamps Optional character vector of labels, one per file, replacing the
#'   dates parsed from the file names.
#' @return List with \code{rast}, \code{time} (fractional calendar year), \code{var}
#'   and \code{stamps}, or \code{NULL} when no dated ascii grids are present.
#' @keywords internal
#' @noRd
.ecospace_stdriver_stack <- function(dir.stdriver, datestamps = NULL){
  f <- list.files(dir.stdriver, pattern = "\\.asc$", full.names = TRUE,
                  recursive = FALSE)
  if(length(f) == 0) return(NULL)
  b <- basename(f)

  # climatology first: its YYYY_MM form would otherwise be misread as YYYYMM
  yr <- mo <- rep(NA_real_, length(b))
  m1 <- regmatches(b, regexec("_([0-9]{4})_([0-9]{2})\\.asc$", b))
  m2 <- regmatches(b, regexec("_([0-9]{4})([0-9]{2})\\.asc$", b))
  for(i in seq_along(b)){
    hit <- if(length(m1[[i]]) == 3L) m1[[i]] else if(length(m2[[i]]) == 3L) m2[[i]] else NULL
    if(!is.null(hit)){ yr[i] <- as.numeric(hit[2]); mo[i] <- as.numeric(hit[3]) }
  }
  keep <- !is.na(yr) & mo >= 1 & mo <= 12
  if(!any(keep)) return(NULL)
  f <- f[keep]; b <- b[keep]; yr <- yr[keep]; mo <- mo[keep]

  # a real-dated series and a placeholder-year climatology usually share a folder;
  # animate the real one unless it is all there is (do.years can select the other)
  if(any(yr >= 1000) && any(yr < 1000)){
    sel <- yr >= 1000
    f <- f[sel]; b <- b[sel]; yr <- yr[sel]; mo <- mo[sel]
  }

  ord <- order(yr, mo)
  f <- f[ord]; b <- b[ord]; yr <- yr[ord]; mo <- mo[ord]

  r <- terra::rast(f)
  time   <- yr + (mo - 1) / 12
  stamps <- if(is.null(datestamps)) sprintf("%04d-%02d", yr, mo) else
            rep_len(datestamps, length(f))
  names(r) <- stamps
  list(rast = r, time = time, var = basename(dir.stdriver), stamps = stamps)
}

#' @title Read a per-region Ecospace timeseries for the animation inset.
#' @description Locates \code{Ecospace_Average_Region_<r>_<Var>.csv} (monthly) or
#'   \code{Ecospace_Annual_Average_Region_<r>_<Var>.csv} (annual) and returns it with a
#'   calendar time vector, so a trajectory can be drawn alongside the map frames.
#'   Region 0 is the whole model domain.
#' @param dir.out Ecospace run output folder.
#' @param variable Variable name as it appears in the file name, e.g. \code{'Biomass'}.
#' @param region Region number; 0 is the full domain.
#' @param timestep \code{'monthly'}, \code{'annual'}, or \code{'auto'} to prefer monthly.
#' @param model_styear Calendar year for model year 0. \code{NULL} reads the header.
#' @return List with \code{time}, \code{data}, \code{timestep} and \code{file}, or
#'   \code{NULL} when no matching file exists.
#' @keywords internal
#' @noRd
.ecospace_region_ts <- function(dir.out, variable, region = 0,
                                timestep = "auto", model_styear = NULL){
  find1 <- function(ts){
    prefix <- if(ts == "annual") "Ecospace_Annual_Average_Region_" else
                                 "Ecospace_Average_Region_"
    f <- list.files(dir.out,
                    pattern = paste0("^", prefix, region, "_", variable, "\\.csv$"),
                    recursive = TRUE, full.names = TRUE)
    if(length(f)) f[1] else NULL
  }
  # "Ecospace_Average_Region_" also prefixes the annual files, so try monthly first
  # and confirm it really carries a TimeStep header before accepting it.
  ts.use <- timestep
  file   <- NULL
  if(timestep %in% c("auto", "monthly")){
    file <- find1("monthly")
    if(!is.null(file) && is.na(.find_year_skip(file, "monthly"))) file <- NULL
    if(!is.null(file)) ts.use <- "monthly"
  }
  if(is.null(file) && timestep %in% c("auto", "annual")){
    file <- find1("annual")
    if(!is.null(file)) ts.use <- "annual"
  }
  if(is.null(file)) return(NULL)

  nskip <- .find_year_skip(file, ts.use)
  if(length(nskip) == 0 || is.na(nskip)) return(NULL)
  dt <- utils::read.csv(file, skip = nskip, check.names = FALSE)
  if(nrow(dt) == 0 || ncol(dt) < 2) return(NULL)

  idx <- suppressWarnings(as.numeric(dt[[1]]))
  dt  <- dt[, -1, drop = FALSE]
  if(is.null(model_styear)){
    hdr <- tryCatch(.ecospace_map_header(file), error = function(e) NULL)
    model_styear <- if(!is.null(hdr)) hdr$StartYear else NULL
  }
  # monthly files index by Ecospace step (12 per model year); annual by model year
  yr <- if(ts.use == "monthly") idx / 12 else idx
  time <- if(!is.null(model_styear) && !is.na(model_styear)) model_styear + yr else yr

  list(time = time, data = dt, timestep = ts.use, file = file)
}

#' @title Pull one group's column out of a per-region timeseries.
#' @description Matches on the group name, falling back to summing \code{fleet|group}
#'   columns whose trailing group part matches (as in the Catch / Landings files).
#' @param dt Data frame of timeseries columns.
#' @param group Group name to extract.
#' @return Numeric vector, or \code{NULL} when the group is absent.
#' @keywords internal
#' @noRd
.ecospace_ts_column <- function(dt, group){
  nm  <- trimws(names(dt))
  hit <- which(tolower(nm) == tolower(group))
  if(length(hit)) return(suppressWarnings(as.numeric(dt[[hit[1]]])))
  tail.nm <- vapply(strsplit(nm, "|", fixed = TRUE),
                    function(x) trimws(x[length(x)]), character(1))
  hit <- which(tolower(tail.nm) == tolower(group))
  if(length(hit))
    return(rowSums(sapply(dt[, hit, drop = FALSE],
                          function(z) suppressWarnings(as.numeric(z))), na.rm = TRUE))
  NULL
}

#' @title Device coordinates for the trajectory inset.
#' @param pos One of 'bottomleft', 'bottomright', 'topleft', 'topright'.
#' @param size Length-2 width and height as fractions of the device.
#' @param map.x Length-2 device x range occupied by the map panel.
#' @return Length-4 vector suitable for \code{par(fig=)}.
#' @keywords internal
#' @noRd
.ecospace_inset_fig <- function(pos, size, map.x = c(0, 0.82)){
  # keep clear of the map's own axis annotation: more room at the bottom, where the
  # longitude labels sit, than at the top, where only the title does
  pad.x <- 0.03; pad.bot <- 0.07; pad.top <- 0.03
  w <- min(size[1], diff(map.x) - 2 * pad.x)
  h <- min(size[2], 1 - pad.bot - pad.top)
  left <- grepl("left", pos)
  bot  <- grepl("bottom", pos)
  x <- if(left) c(map.x[1] + pad.x, map.x[1] + pad.x + w)
       else     c(map.x[2] - pad.x - w, map.x[2] - pad.x)
  y <- if(bot)  c(pad.bot, pad.bot + h) else c(1 - pad.top - h, 1 - pad.top)
  c(x, y)
}

#' @title Resolve a color palette given by name or by value.
#' @description Lets \code{plt.col} be either an explicit vector of colors or the name
#'   of a palette: \code{'matlab'} / \code{'jet'} (via \code{.matlab_palette()}, which
#'   falls back to an equivalent jet ramp when \code{colorRamps} is unavailable), any
#'   \code{viridis} option, or any ramp exported by \code{colorRamps}.
#' @param x \code{NULL}, a palette name, or a vector of colors.
#' @param n Number of colors to generate when \code{x} names a palette.
#' @return A character vector of colors.
#' @keywords internal
#' @noRd
.ecospace_palette <- function(x, n = 99){
  if(is.null(x)) return(viridis::viridis(n))

  # a ramp function: colorRampPalette() style takes n and returns colors, while
  # colorRamp() style takes values in [0, 1] and returns an RGB 0-255 matrix
  if(is.function(x)){
    out <- tryCatch(x(n), error = function(e) NULL)
    if(is.character(out) && length(out) == n) return(out)
    m <- x(seq(0, 1, length.out = n))
    if(is.matrix(m) && ncol(m) >= 3 && nrow(m) == n)
      return(grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255))
    stop("'plt.col' is a function but returned neither ", n, " colors nor an ",
         n, "x3 RGB matrix. Pass colorRampPalette(...) or colorRamp(...).",
         call. = FALSE)
  }

  if(length(x) > 1 || !is.character(x)) return(x)   # already a color vector

  key <- tolower(x)
  if(key %in% c("matlab", "matlab.like", "jet")) return(.matlab_palette(n))
  if(key %in% c("viridis", "magma", "plasma", "inferno", "cividis",
                "rocket", "mako", "turbo"))
    return(viridis::viridis(n, option = key))
  if(requireNamespace("colorRamps", quietly = TRUE)){
    ex  <- getNamespaceExports("colorRamps")
    hit <- ex[tolower(ex) == key]
    if(length(hit)) return(get(hit[1], envir = asNamespace("colorRamps"))(n))
  }
  # a single valid color name is a (degenerate but legal) color vector
  if(isTRUE(tryCatch({grDevices::col2rgb(x); TRUE}, error = function(e) FALSE)))
    return(x)
  stop("Unrecognised 'plt.col': ", x, ". Use a vector of colors, or one of ",
       "'matlab' / 'jet', a viridis option ('viridis', 'magma', 'turbo', ...), ",
       "or a colorRamps ramp ('matlab.like2', 'blue2red', ...).", call. = FALSE)
}

#' @title Animate Ecospace spatial output as a gif.
#' @description Builds one movie (gif) per variable/group series from an Ecospace run's
#'   spatial output. Two output formats are supported and detected automatically: the
#'   per-group map csv files (\code{EcospaceMap<Variable>-<group>.csv}, written to a
#'   \code{csv/} subfolder by the GUI and to the top level of the run folder by the
#'   console), and the per-timestep ascii grids in an \code{asc/} folder, read via
#'   \code{\link{fn.ecospace_ascii2stack}}. All frames in a series share one color
#'   scale, so change across the animation reflects real change rather than rescaling.
#' @details For \code{source = 'stdriver'}, \code{dir.out} is a single ST driver
#'   variable's folder -- Ecospace \emph{input}, a flat set of monthly ascii grids --
#'   rather than a run output folder. There is one series, so \code{do.vars},
#'   \code{do.grps}, \code{group.names} and \code{fleet.names} do not apply and are
#'   ignored with a warning; everything else (color scale, zero bin, coastlines,
#'   \code{do.years} / \code{do.steps}, inset) works the same. Frames are dated from
#'   the file names, \code{<var>_YYYYMM.asc}, and a 12-month climatology written with a
#'   placeholder year (\code{<var>_YYYY_MM.asc}, e.g. \code{_0000_01}) is set aside
#'   when a real-dated series is present in the same folder -- select it with
#'   \code{do.years} if you want to animate it instead. Only the folder itself is
#'   scanned, so a \code{static/} subfolder of long-term means is ignored. The gif is
#'   written beside the driver folder rather than inside it. With \code{inset = TRUE}
#'   the trajectory is the driver's own domain mean, since drivers have no per-region
#'   timeseries file.
#' @param dir.out Ecospace run output folder, or a single ST driver variable's folder
#'   when \code{source} resolves to \code{'stdriver'}. Defaults to \code{dir.pred}.
#' @param do.vars Character vector of output variables to animate, matching the EwE map
#'   file naming: \code{'Biomass'}, \code{'Catch'}, \code{'Effort'}. Matched
#'   case-insensitively, so the lowercase names used by
#'   \code{\link{fn.ecospace_ascii2stack}} also work.
#' @param do.grps Groups to animate, or fleets when \code{do.vars} is \code{'Effort'}
#'   (effort maps are written per fleet, not per functional group). Either a character
#'   vector of names as they appear in the map file names, or a numeric vector of
#'   indices -- resolved against \code{fleet.names} for \code{'Effort'} and against
#'   \code{group.names} for every other variable. \code{NULL} (default) animates every
#'   group or fleet found for each requested variable, which can be a lot of gifs --
#'   subset it.
#' @param do.years Calendar years to keep as frames. \code{NULL} (default) uses every
#'   timestep in the output.
#' @param do.steps Optional numeric vector indexing which of the retained timesteps to
#'   plot, applied after \code{do.years}. Useful for thinning a monthly run to a
#'   manageable frame count, e.g. \code{do.steps = seq(1, 241, by = 12)}.
#' @param group.names Character vector of group names indexed by pool code, used to
#'   resolve a numeric \code{do.grps} and to label \code{asc} series. \code{NULL}
#'   (default) falls back to \code{df.names$group.names} if that object is visible.
#' @param fleet.names Character vector of fleet names indexed by fleet number, used the
#'   same way for \code{'Effort'}. \code{NULL} (default) falls back to
#'   \code{df.names$fleet.names} if that object is visible.
#' @param source One of \code{'auto'} (default), \code{'csv'}, \code{'asc'} or
#'   \code{'stdriver'}, selecting which format to read. \code{'stdriver'} animates an
#'   Ecospace \emph{input} driver folder rather than model output; see Details.
#' @param transform Optional function taking and returning a \code{terra} SpatRaster,
#'   applied to the whole stack before plotting. Use it for unit conversions, e.g.
#'   \code{transform = function(x) x * (365*24*60*60) * 14.007} to take NPP from
#'   molN/m2/s to gN/m2/yr. Pair it with \code{legend.lab} to label the new units.
#'   Currently applied to the \code{'stdriver'} source only.
#' @param model_styear Calendar year corresponding to model year 0, used to date the
#'   frames. \code{NULL} (default) reads \code{StartYear} from the map csv header.
#' @param zlim Numeric length-2 color scale limits shared by every frame. \code{NULL}
#'   (default) uses the range of the series being animated.
#' @param zlim.probs Optional length-2 probability vector, e.g. \code{c(0, 0.99)}, used
#'   in place of the full range so a few extreme cells do not flatten the map. Ignored
#'   when \code{zlim} is supplied. Cells outside the resulting scale are clamped to the
#'   end colors rather than dropped, so they stay visually distinct from land.
#' @param plt.col Colors for the graded scale, given as any of: an explicit vector of
#'   colors; a ramp function, either \code{colorRampPalette()} style (takes \code{n},
#'   returns colors) or \code{colorRamp()} style (takes \code{[0,1]}, returns an RGB
#'   matrix); or the name of a palette -- \code{'matlab'} (equivalently \code{'jet'})
#'   for the MATLAB jet look, any \code{viridis} option (\code{'viridis'},
#'   \code{'magma'}, \code{'turbo'}, ...), or any ramp exported by \code{colorRamps}
#'   (\code{'matlab.like2'}, \code{'green2red'}, ...). \code{NULL} (default) uses
#'   \code{viridis}. A custom ramp is the way to get an arbitrary sequence, e.g.
#'   \code{colorRampPalette(c("darkgreen", "orange", "red"))} for a severity scale.
#' @param n.col Number of colors to generate when \code{plt.col} names a palette
#'   (default 99). Ignored when \code{plt.col} is given as an explicit color vector.
#' @param zero.col Optional single color for cells at or below \code{zero.thresh},
#'   drawn as a separate bin below the graded \code{plt.col} ramp. Use it to knock
#'   near-zero cells out to a flat background color, e.g. \code{zero.col = "white"}, so
#'   occupied cells read against empty ones. \code{NULL} (default) uses one continuous
#'   ramp over the whole range. Note that the graded ramp is then stretched over
#'   \code{(zero.thresh, max]} rather than the full range, so the colors of non-zero
#'   cells change too.
#' @param zero.thresh Upper edge of the \code{zero.col} bin, in data units. Default
#'   \code{0} colors only exactly-zero cells; raise it (e.g. \code{1e-4}) to catch
#'   near-zero cells. Ignored when \code{zero.col} is \code{NULL}.
#' @param log.scale Logical; space the graded color bins logarithmically instead of
#'   linearly. Default \code{FALSE}. Use it for heavily skewed fields such as red tide,
#'   where a few extreme cells otherwise flatten everything else. Only the \emph{colors}
#'   are logged -- the data are not transformed and the legend stays in native units,
#'   with ticks at the decades. Combines with \code{zero.col}: zeros keep their own flat
#'   bin below the log ramp, which is what you want for a field that is exactly zero
#'   much of the time.
#' @param log.floor Lower bound of the log ramp, in data units. \code{NULL} (default)
#'   uses \code{zero.thresh} when that is positive, otherwise the smallest positive
#'   value in the series. Set it to stop a single tiny value from stretching the ramp
#'   over more decades than are worth showing. Ignored when \code{log.scale} is
#'   \code{FALSE}.
#' @param na.col Color for NA cells -- land and inactive water. Default \code{NA}
#'   leaves them transparent, which shows as white. Set it (e.g. \code{'gray95'}) when
#'   using a white \code{zero.col}, so empty model cells stay distinguishable from
#'   cells that are outside the model domain entirely.
#' @param legend.lab Legend label. \code{NULL} (default) uses the variable name.
#' @param add.map Logical; overlay world and state coastlines. \code{NULL} (default)
#'   adds them only when the raster extent looks like unprojected lon/lat.
#' @param inset Logical; draw a small timeseries panel inside the map showing the
#'   group's trajectory for \code{inset.region}, with a marker that advances in step
#'   with the frames. Default \code{FALSE}. The marker is placed at the frame's own
#'   calendar time, so map and trajectory stay in sync by construction even when
#'   \code{do.years} / \code{do.steps} thin the frames.
#' @param inset.region Region number for the inset timeseries. Default \code{0}, the
#'   whole model domain.
#' @param inset.var Variable for the inset timeseries. \code{NULL} (default) uses the
#'   mapped variable. Effort has no per-region file, so set this (e.g. to
#'   \code{'Catch'}) when animating effort, or the inset is skipped with a warning.
#' @param inset.timestep \code{'auto'} (default, prefers the monthly file),
#'   \code{'monthly'} or \code{'annual'}.
#' @param inset.pos Corner of the map to place the inset in: \code{'bottomleft'}
#'   (default), \code{'bottomright'}, \code{'topleft'} or \code{'topright'}.
#' @param inset.size Length-2 width and height of the inset as fractions of the frame.
#' @param inset.col Color of the trajectory line.
#' @param inset.mark Color of the marker showing the current frame's position.
#' @param dir.gif Folder to write gifs into. Defaults to a \code{gifs} folder inside
#'   \code{dir.out}.
#' @param width,height,res Pixel dimensions and resolution of the frames.
#' @param delay Seconds between frames.
#' @param loop Logical; loop the gif forever (default) or play it once.
#' @param keep.png Logical; keep the intermediate png frames. Default \code{FALSE}
#'   deletes the temporary png folder.
#' @param datestamps Optional character vector of frame labels, one per retained
#'   timestep, overriding the dates derived from the output.
#' @param verbose Logical; report progress.
#' @return Invisibly, a character vector of the gif files written.
#' @examples
#' # example code:
#' \dontrun{
#' # every gag stanza, thinned to one frame per year
#' fn.ecospace_gif(dir.out  = dir.pred,
#'                 do.vars  = "Biomass",
#'                 do.grps  = c("gag 0","gag 1","gag 2","gag 3","gag 4","gag 5+"),
#'                 do.steps = seq(1, 241, by = 12))
#'
#' # one group, robust color scale, 2010-2020 only
#' fn.ecospace_gif(dir.out = dir.pred, do.vars = "Biomass", do.grps = "snappers",
#'                 do.years = 2010:2020, zlim.probs = c(0, 0.99))
#'
#' # an ST driver (Ecospace input), with its domain-mean trajectory as an inset
#' fn.ecospace_gif(dir.out = file.path(dir.stdriver, "btm_temp"),
#'                 plt.col = "matlab", inset = TRUE)
#'
#' # red tide: log color scale, exact zeros knocked out to white
#' fn.ecospace_gif(dir.out = file.path(dir.stdriver, "red tide/sdmTMB"),
#'                 plt.col = "matlab", log.scale = TRUE, log.floor = 1e3,
#'                 zero.col = "white", na.col = "gray95", inset = TRUE)
#' }
#' @seealso \code{\link{fn.ecospace_ascii2stack}}, \code{\link{fn.STdriver_gif}}
#' @export
fn.ecospace_gif <- function(dir.out      = dir.pred,
                            do.vars      = "Biomass",
                            do.grps      = NULL,
                            do.years     = NULL,
                            do.steps     = NULL,
                            group.names  = NULL,
                            fleet.names  = NULL,
                            source       = c("auto", "csv", "asc", "stdriver"),
                            transform    = NULL,
                            model_styear = NULL,
                            zlim         = NULL,
                            zlim.probs   = NULL,
                            plt.col      = NULL,
                            n.col        = 99,
                            zero.col     = NULL,
                            zero.thresh  = 0,
                            na.col       = NA,
                            log.scale    = FALSE,
                            log.floor    = NULL,
                            legend.lab   = NULL,
                            add.map      = NULL,
                            inset        = FALSE,
                            inset.region = 0,
                            inset.var    = NULL,
                            inset.timestep = c("auto", "monthly", "annual"),
                            inset.pos    = c("bottomleft", "bottomright",
                                             "topleft", "topright"),
                            inset.size   = c(0.34, 0.26),
                            inset.col    = "gray30",
                            inset.mark   = "red",
                            dir.gif      = NULL,
                            width        = 1400,
                            height       = 1000,
                            res          = 150,
                            delay        = 0.15,
                            loop         = TRUE,
                            keep.png     = FALSE,
                            datestamps   = NULL,
                            verbose      = TRUE){

  # ---- Input validation ----
  if(missing(dir.out) && !exists("dir.pred"))
    stop("Provide 'dir.out', the Ecospace run output folder.", call. = FALSE)
  if(is.null(dir.out) || !dir.exists(dir.out))
    stop("Provide a valid 'dir.out' folder. Got: ", dir.out, call. = FALSE)
  if(!is.null(zlim) && (length(zlim) != 2 || any(!is.finite(zlim))))
    stop("'zlim' must be a finite numeric vector of length 2.", call. = FALSE)
  if(!is.null(zlim.probs) && (length(zlim.probs) != 2 || any(zlim.probs < 0 | zlim.probs > 1)))
    stop("'zlim.probs' must be two probabilities in [0, 1].", call. = FALSE)

  # ---- Which output format? ----
  source  <- match.arg(source)
  map.df  <- .ecospace_map_files(dir.out)
  has.asc <- dir.exists(file.path(dir.out, "asc")) &&
             length(list.files(file.path(dir.out, "asc"), pattern = "\\.asc$",
                               recursive = TRUE)) > 0
  # an ST driver folder holds its dated ascii grids directly, not under asc/
  has.std <- length(list.files(dir.out, pattern = "\\.asc$", recursive = FALSE)) > 0
  if(source == "auto"){
    source <- if(!is.null(map.df)) "csv" else if(has.asc) "asc" else
              if(has.std) "stdriver" else NA_character_
    if(is.na(source))
      stop("No spatial data found in ", dir.out,
           ". Expected Ecospace output (EcospaceMap*.csv files, or an asc/ folder of ",
           "ascii grids) or an ST driver folder of dated monthly .asc grids.",
           call. = FALSE)
  }
  if(source == "stdriver" && !has.std)
    stop("No .asc grids directly in ", dir.out,
         ". For an ST driver, point 'dir.out' at a single variable's folder.",
         call. = FALSE)
  if(source == "stdriver"){
    for(a in c("do.vars", "do.grps", "group.names", "fleet.names"))
      if(!eval(call("missing", as.name(a))))
        warning("'", a, "' does not apply to an ST driver folder and is ignored.",
                call. = FALSE, immediate. = TRUE)
    # one variable per folder, so the gif belongs beside it, not inside it
    if(is.null(dir.gif)) dir.gif <- file.path(dirname(dir.out), "gifs")
  }
  if(source == "csv" && is.null(map.df))
    stop("No EcospaceMap*.csv files found in ", dir.out, " or its csv/ subfolder.",
         call. = FALSE)
  if(source == "asc" && !has.asc)
    stop("No .asc files found in ", file.path(dir.out, "asc"), call. = FALSE)

  # Numeric selections index different name vectors depending on the variable:
  # effort maps are per fleet, everything else is per functional group. Resolve
  # against both when both kinds of variable were requested.
  want.effort <- any(tolower(do.vars) == "effort")
  want.group  <- any(tolower(do.vars) != "effort")
  if(is.null(group.names) && exists("df.names"))
    group.names <- get("df.names")$group.names
  if(is.null(fleet.names) && exists("df.names"))
    fleet.names <- get("df.names")$fleet.names
  if(is.numeric(do.grps)){
    if(want.group && is.null(group.names))
      stop("Numeric 'do.grps' needs 'group.names' (or a 'df.names' object) to name ",
           "the functional groups.", call. = FALSE)
    if(want.effort && is.null(fleet.names))
      stop("Numeric 'do.grps' selects fleets when do.vars includes 'Effort', so it ",
           "needs 'fleet.names' (or a 'df.names' object). Naming the fleets directly, ",
           "e.g. do.grps = 'shrimp fishery', avoids this.", call. = FALSE)
    grp.idx <- do.grps
    nm <- character(0)
    if(want.group)  nm <- c(nm, group.names[do.grps])
    if(want.effort) nm <- c(nm, fleet.names[do.grps])
    do.grps <- gsub("_", " ", unique(nm[!is.na(nm)]))
    if(length(do.grps) == 0)
      stop("Numeric 'do.grps' did not index any name in ",
           if(want.effort) "'fleet.names'" else "'group.names'", ".", call. = FALSE)
  } else {
    grp.idx <- NULL
  }

  if(is.null(dir.gif)) dir.gif <- file.path(dir.out, "gifs")
  if(!dir.exists(dir.gif)) dir.create(dir.gif, recursive = TRUE)
  dir.png <- file.path(dir.gif, "png_tmp")
  if(!dir.exists(dir.png)) dir.create(dir.png, recursive = TRUE)
  if(!keep.png) on.exit(unlink(dir.png, recursive = TRUE, force = TRUE), add = TRUE)

  # ---- Assemble the series to animate ----
  series <- list()

  if(source == "stdriver"){
    if(verbose) message("Reading ST driver grids in ", basename(dir.out))
    sd <- .ecospace_stdriver_stack(dir.out)   # 'datestamps' is applied per frame below
    if(is.null(sd))
      stop("No dated .asc grids in ", dir.out,
           ". Expected names ending '_YYYYMM.asc' or '_YYYY_MM.asc'.", call. = FALSE)
    if(!is.null(transform)){
      if(!is.function(transform))
        stop("'transform' must be a function taking and returning a SpatRaster.",
             call. = FALSE)
      sd$rast <- transform(sd$rast)
    }
    series[[1]] <- list(var = sd$var, group = "", rast = sd$rast, time = sd$time,
                        stamps = sd$stamps)
  } else if(source == "csv"){
    keep <- tolower(map.df$variable) %in% tolower(do.vars)
    if(!any(keep))
      stop("No map files match do.vars = ", paste(do.vars, collapse = ", "),
           ". Available: ", paste(unique(map.df$variable), collapse = ", "),
           call. = FALSE)
    map.df <- map.df[keep, , drop = FALSE]
    if(!is.null(do.grps)){
      keep <- tolower(map.df$group) %in% tolower(do.grps)
      missed <- setdiff(tolower(do.grps), tolower(map.df$group))
      if(!any(keep))
        stop("No map files match do.grps. Available groups for the requested ",
             "variable(s): ", paste(unique(map.df$group), collapse = ", "), call. = FALSE)
      if(length(missed))
        warning("No map file for group(s): ", paste(missed, collapse = ", "),
                call. = FALSE, immediate. = TRUE)
      map.df <- map.df[keep, , drop = FALSE]
    }
    for(i in seq_len(nrow(map.df))){
      if(verbose) message("Reading ", basename(map.df$file[i]))
      r <- .ecospace_map_stack(map.df$file[i], model_styear = model_styear)
      series[[length(series) + 1L]] <- list(var   = map.df$variable[i],
                                            group = map.df$group[i],
                                            rast  = r,
                                            time  = attr(r, "time"))
    }
  } else {
    # ascii grids: reuse the existing reader, then split its stack into series
    all.idx <- 1:9999
    stk <- fn.ecospace_ascii2stack(dir.out    = dir.out,
                                   do.vars    = tolower(do.vars),
                                   do.grps    = if(!want.group) integer(0) else
                                                if(is.null(grp.idx)) all.idx else grp.idx,
                                   do.fleets  = if(!want.effort) integer(0) else
                                                if(is.null(grp.idx)) all.idx else grp.idx,
                                   do.months  = 1:12,
                                   do.years   = if(is.null(do.years)) 1:9999 else do.years,
                                   annualmaps = FALSE)
    lab   <- names(stk)                       # type_NN_YYYY_MM
    parts <- regmatches(lab, regexec("^(.*)_([0-9]+)_([0-9]{4})_([0-9]{2})$", lab))
    ok    <- lengths(parts) == 5L
    if(!any(ok))
      stop("Could not parse layer names returned by fn.ecospace_ascii2stack().",
           call. = FALSE)
    idx <- which(ok)
    typ <- vapply(parts[ok], function(x) x[2], character(1))
    num <- vapply(parts[ok], function(x) x[3], character(1))
    yr  <- as.numeric(vapply(parts[ok], function(x) x[4], character(1)))
    mo  <- as.numeric(vapply(parts[ok], function(x) x[5], character(1)))
    tm  <- yr + (mo - 1) / 12
    key <- paste(typ, num, sep = "_")
    for(k in unique(key)){
      sel <- which(key == k)
      ord <- sel[order(tm[sel])]
      n1  <- as.integer(num[ord[1]])
      # effort series are numbered by fleet, everything else by functional group
      nms   <- if(tolower(typ[ord[1]]) == "effort") fleet.names else group.names
      gname <- if(!is.null(nms) && !is.na(nms[n1])) gsub("_", " ", nms[n1]) else num[ord[1]]
      # a character do.grps is not seen by fn.ecospace_ascii2stack(), so apply it here
      if(!is.null(do.grps) && is.null(grp.idx) && !(tolower(gname) %in% tolower(do.grps)))
        next
      series[[length(series) + 1L]] <- list(var   = typ[ord[1]],
                                            group = gname,
                                            rast  = stk[[idx[ord]]],
                                            time  = tm[ord])
    }
    if(length(series) == 0)
      stop("No .asc series matched do.grps. Name the group/fleet as it appears in ",
           "'group.names' / 'fleet.names', or pass numeric indices.", call. = FALSE)
    do.years <- NULL   # already applied inside fn.ecospace_ascii2stack()
  }

  # ---- Render ----
  plt.col        <- .ecospace_palette(plt.col, n.col)
  inset.pos      <- match.arg(inset.pos)
  inset.timestep <- match.arg(inset.timestep)
  ts.cache       <- list()   # one read per variable, reused across groups
  gif.files      <- character(0)

  for(s in series){
    r     <- s$rast
    tm    <- s$time
    label <- trimws(paste(s$var, s$group))   # ST driver series carry no group

    # frame subsetting
    keep <- if(is.null(do.years)) seq_along(tm) else which(floor(tm + 1e-8) %in% do.years)
    if(length(keep) == 0){
      warning("No timesteps in do.years for ", label, "; skipped.",
              call. = FALSE, immediate. = TRUE)
      next
    }
    if(!is.null(do.steps)){
      sub.steps <- do.steps[do.steps >= 1 & do.steps <= length(keep)]
      if(length(sub.steps) == 0){
        warning("'do.steps' selected no frames for ", label, "; skipped.",
                call. = FALSE, immediate. = TRUE)
        next
      }
      keep <- keep[sub.steps]
    }
    r  <- r[[keep]]
    tm <- tm[keep]
    if(terra::nlyr(r) < 2)
      warning("Only ", terra::nlyr(r), " frame for ", label,
              "; the gif will be a still.", call. = FALSE, immediate. = TRUE)

    stamps <- if(!is.null(datestamps)) rep_len(datestamps, terra::nlyr(r)) else
              if(!is.null(s$stamps)) s$stamps[keep] else .ecospace_datestamp(tm)

    # one color scale for the whole series
    if(!is.null(zlim)){
      plt.zlims <- zlim
    } else if(!is.null(zlim.probs)){
      plt.zlims <- as.numeric(stats::quantile(terra::values(r), probs = zlim.probs,
                                              na.rm = TRUE))
    } else {
      plt.zlims <- c(min(terra::global(r, "min", na.rm = TRUE)[, 1], na.rm = TRUE),
                     max(terra::global(r, "max", na.rm = TRUE)[, 1], na.rm = TRUE))
    }
    if(!all(is.finite(plt.zlims))){
      warning("No finite values for ", label, "; skipped.", call. = FALSE,
              immediate. = TRUE)
      next
    }
    if(diff(plt.zlims) <= 0)
      plt.zlims <- plt.zlims + c(-0.5, 0.5) * max(1e-9, abs(plt.zlims[1]))
    # optional flat bin at the bottom, so near-zero cells read as background
    use.zero <- !is.null(zero.col)
    if(use.zero && (!is.finite(zero.thresh) || zero.thresh >= plt.zlims[2])){
      warning("'zero.thresh' is not below the top of the color scale for ", label,
              "; ignoring 'zero.col'.", call. = FALSE, immediate. = TRUE)
      use.zero <- FALSE
    }

    # A log scale needs a positive lower bound for the graded ramp: zero and negative
    # cells cannot sit on it. Take it from 'log.floor', else the zero bin's own upper
    # edge, else the smallest positive value in the series.
    use.log <- isTRUE(log.scale)
    lo <- NA_real_
    if(use.log){
      lo <- log.floor
      if(is.null(lo)){
        lo <- if(use.zero && zero.thresh > 0) zero.thresh else {
          v <- terra::values(r); v <- v[is.finite(v) & v > 0]
          if(length(v)) min(v) else NA_real_
        }
      }
      if(!is.finite(lo) || lo <= 0 || lo >= plt.zlims[2]){
        warning("No positive values below the maximum for ", label,
                "; drawing it on a linear scale. Set 'log.floor' to force one.",
                call. = FALSE, immediate. = TRUE)
        use.log <- FALSE
      }
    }

    # plt.brks are in data units (the map); leg.brks/leg.zlim drive the color bar,
    # which is uniform in log space when logging so each decade gets equal width
    if(use.log){
      lg   <- seq(log10(lo), log10(plt.zlims[2]), length.out = length(plt.col) + 1L)
      ramp <- 10^lg
      if(use.zero){
        clamp.lo <- min(plt.zlims[1], zero.thresh)
        cols.use <- c(zero.col, plt.col)
        plt.brks <- c(clamp.lo, ramp)
        leg.brks <- c(lg[1] - diff(lg)[1], lg)
      } else {
        clamp.lo <- lo          # sub-floor cells clamp onto the bottom color
        cols.use <- plt.col
        plt.brks <- ramp
        leg.brks <- lg
      }
      leg.zlim <- range(leg.brks)
      # tick the decades, but label them in native units
      d  <- seq(floor(log10(lo)), ceiling(log10(plt.zlims[2])))
      tk <- 10^d
      if(sum(tk >= lo & tk <= plt.zlims[2]) < 3) tk <- sort(c(tk, 3 * 10^d))
      tk <- tk[tk >= lo & tk <= plt.zlims[2]]
      tick_at  <- log10(tk)
      tick_lab <- formatC(tk, format = "g", digits = 3)
    } else {
      if(use.zero){
        clamp.lo <- min(plt.zlims[1], zero.thresh)
        cols.use <- c(zero.col, plt.col)
        plt.brks <- c(clamp.lo,
                      seq(zero.thresh, plt.zlims[2], length.out = length(plt.col) + 1))
      } else {
        clamp.lo <- plt.zlims[1]
        cols.use <- plt.col
        plt.brks <- seq(plt.zlims[1], plt.zlims[2], length.out = length(plt.col) + 1)
      }
      leg.brks <- NULL         # set below, after the breaks are made strictly increasing
      leg.zlim <- plt.zlims
      tick_at  <- pretty(plt.zlims, n = 6)
      tick_at  <- tick_at[tick_at >= plt.zlims[1] & tick_at <= plt.zlims[2]]
      tick_lab <- format(tick_at, digits = 3, scientific = FALSE)
    }
    # widen the outer breaks so cells sitting exactly on a limit still get a color
    # (a zero bin whose edge equals the data minimum leaves a duplicate break here,
    # so the legend must take its copy after this, not before)
    eps <- diff(range(plt.brks)) * 1e-6
    plt.brks[1] <- plt.brks[1] - eps
    plt.brks[length(plt.brks)] <- plt.brks[length(plt.brks)] + eps
    if(is.null(leg.brks)) leg.brks <- plt.brks

    # coastlines only make sense on unprojected lon/lat
    e <- terra::ext(r)
    do.map <- if(is.null(add.map))
      (e[1] >= -180 && e[2] <= 180 && e[3] >= -90 && e[4] <= 90) else isTRUE(add.map)

    # trajectory for the inset, if requested and available for this series
    ts.x <- NULL; ts.y <- NULL; ts.lab <- NULL
    if(isTRUE(inset) && source == "stdriver"){
      # no per-region file for a driver, so the trajectory is its own domain mean
      yv <- as.numeric(terra::global(s$rast, "mean", na.rm = TRUE)[, 1])
      if(all(is.na(yv))){
        warning("ST driver ", label, " is all NA; drawing without an inset.",
                call. = FALSE, immediate. = TRUE)
      } else {
        ts.x   <- s$time
        ts.y   <- yv
        ts.lab <- paste(s$var, "mean")
      }
    } else if(isTRUE(inset)){
      iv <- if(is.null(inset.var)) s$var else inset.var
      ck <- paste(iv, inset.region, inset.timestep, sep = "|")
      if(!ck %in% names(ts.cache))
        ts.cache[[ck]] <- list(v = .ecospace_region_ts(dir.out, iv, inset.region,
                                                       inset.timestep, model_styear))
      ts <- ts.cache[[ck]]$v
      if(is.null(ts)){
        warning("No Region ", inset.region, " ", iv, " timeseries found in ", dir.out,
                "; drawing ", label, " without an inset.",
                call. = FALSE, immediate. = TRUE)
      } else {
        yv <- .ecospace_ts_column(ts$data, s$group)
        if(is.null(yv)){
          warning("Group '", s$group, "' is not a column of ", basename(ts$file),
                  "; drawing without an inset.", call. = FALSE, immediate. = TRUE)
        } else if(all(is.na(yv))){
          warning("Region ", inset.region, " ", iv, " is all NA for '", s$group,
                  "'; drawing without an inset.", call. = FALSE, immediate. = TRUE)
        } else {
          ts.x   <- ts$time
          ts.y   <- yv
          ts.lab <- paste0(iv, " R", inset.region)
        }
      }
    }

    slug      <- gsub("_+$", "", gsub("[^A-Za-z0-9._+-]+", "_",
                                      trimws(paste(s$var, s$group))))
    png.files <- character(terra::nlyr(r))
    if(verbose) message("Rendering ", terra::nlyr(r), " frames for ", label)

    for(i in seq_len(terra::nlyr(r))){
      png.files[i] <- file.path(dir.png, sprintf("%s_%04d.png", slug, i))
      # clamp to the shared scale so out-of-range cells keep the end color instead
      # of dropping out as transparent, which would look like land / NODATA
      ri <- terra::clamp(r[[i]], clamp.lo, plt.zlims[2], values = TRUE)
      grDevices::png(png.files[i], width = width, height = height, res = res)
      # two-column layout: big map + skinny legend
      graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.22))
      graphics::par(mar = c(2, 0.1, 2, 0.1))
      terra::plot(ri, col = cols.use, breaks = plt.brks, legend = FALSE,
                  colNA = na.col, main = paste(label, stamps[i]),
                  mar = c(2, 0, 2, 0))
      # terra::plot silently ignores 'breaks' (and 'range') when a layer holds a
      # single value, painting every cell the LAST color -- so an all-zero month
      # would come out as the scale maximum. Repaint the cells through
      # graphics::image(), which honors breaks in every case. NA cells are left
      # untouched, so terra's colNA shows through.
      mm <- terra::as.matrix(ri, wide = TRUE)
      ex <- as.vector(terra::ext(ri))
      graphics::image(x = seq(ex[1], ex[2], length.out = ncol(mm) + 1L),
                      y = seq(ex[3], ex[4], length.out = nrow(mm) + 1L),
                      z = t(mm)[, rev(seq_len(nrow(mm))), drop = FALSE],
                      col = cols.use, breaks = plt.brks, add = TRUE)
      if(do.map){
        maps::map(database = "world", add = TRUE, fill = TRUE, col = "gray")
        maps::map(database = "state", add = TRUE, fill = FALSE)
      }
      graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
      fields::image.plot(legend.only = TRUE,
                         zlim        = leg.zlim,
                         col         = cols.use,
                         breaks      = leg.brks,
                         axis.args   = list(at = tick_at, labels = tick_lab,
                                            cex.axis = 0.85),
                         legend.lab  = if(is.null(legend.lab)) s$var else legend.lab,
                         legend.cex  = 0.9,
                         legend.line = 3.5,
                         legend.mar  = 6,
                         horizontal  = FALSE)

      # trajectory inset, drawn last so it floats over the map
      if(!is.null(ts.x)){
        op <- graphics::par(no.readonly = TRUE)
        fg <- .ecospace_inset_fig(inset.pos, inset.size)
        # opaque card over the whole inset region, margins included, so the map's
        # own axis annotation cannot show through the panel
        graphics::par(fig = fg, new = TRUE, mar = c(0, 0, 0, 0))
        graphics::plot.new()
        graphics::rect(0, 0, 1, 1, col = "white", border = "gray40", lwd = 0.8)
        graphics::par(fig = fg, new = TRUE, mar = c(1.9, 3.4, 0.5, 1.2),
                      mgp = c(1.6, 0.4, 0), tcl = -0.25,
                      cex.axis = 0.6, cex.lab = 0.65)
        graphics::plot(ts.x, ts.y, type = "n", xlab = "", ylab = "", axes = FALSE)
        graphics::lines(ts.x, ts.y, col = inset.col, lwd = 1.4)
        # marker sits at this frame's own calendar time, so it cannot drift
        yi <- stats::approx(ts.x, ts.y, xout = tm[i], rule = 2)$y
        graphics::abline(v = tm[i], col = inset.mark, lwd = 1)
        graphics::points(tm[i], yi, pch = 19, col = inset.mark, cex = 0.9)
        graphics::axis(1, cex.axis = 0.6)
        graphics::axis(2, at = pretty(range(ts.y, na.rm = TRUE), n = 3),
                       cex.axis = 0.6, las = 1)
        graphics::title(ylab = ts.lab, line = 2.4, cex.lab = 0.65)
        graphics::box()
        graphics::par(op)
      }
      grDevices::dev.off()
    }

    gif.out <- file.path(dir.gif, paste0(slug, ".gif"))
    gifski::gifski(png.files, gif_file = gif.out, width = width, height = height,
                   delay = delay, loop = loop, progress = FALSE)
    if(verbose) message("Wrote ", gif.out)
    gif.files <- c(gif.files, gif.out)
    if(!keep.png) unlink(png.files, force = TRUE)
  }

  invisible(gif.files)
} #eof

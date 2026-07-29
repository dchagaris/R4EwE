#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spatial likelihood component for Ecospace GA calibration (phase 2).
#
# Adds a per-species spatial fit term to fn.ecospace_objfxn: normalised
# observed maps (e.g. GFISHER maxN heatmaps) vs. mean-of-years predicted
# biomass maps parsed from Ecospace's EcospaceMap<Var>-<species>.csv files.
#
# Design notes:
#   - Obs maps are per (species / pool code) and carry their own year window
#     (styear, endyear) so future obs collected over different time periods
#     can be compared against year-matched model output.
#   - Prediction reader uses the per-species CSV (denser than the per-step
#     .asc files).
#   - LL formulation: cross-entropy -sum(obs * log(pred)) on unit-sum
#     normalised maps over the shared valid (non-NODATA) mask. This is the
#     multinomial negative log-likelihood kernel (constant term -sum(o log o)
#     absorbed), so the spatial component is additive with the lognormal
#     timeseries NLL. Scaled by a single myconfig$spatial.weight (interpretable
#     as a prior precision on the spatial fit).
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Read an ASCII grid, returning a numeric matrix + NODATA value.
#' @description Skips the 6-line ESRI header, coerces cells to numeric, and
#'   returns the matrix with attribute \code{"nodata"}.
#' @param path Path to `.asc` file.
#' @return Matrix (nrows x ncols) with `nodata` attribute (numeric).
#' @keywords internal
#' @noRd
.read_asc <- function(path){
  con <- file(path, "r"); on.exit(close(con), add = TRUE)
  hdr <- character(6)
  for(i in seq_len(6)) hdr[i] <- readLines(con, n = 1)
  parse_num <- function(key){
    ln <- hdr[grep(paste0("^\\s*", key), hdr, ignore.case = TRUE)][1]
    as.numeric(sub(paste0("(?i)^\\s*", key, "\\s+"), "", ln, perl = TRUE))
  }
  ncols  <- as.integer(parse_num("ncols"))
  nrows  <- as.integer(parse_num("nrows"))
  nodata <- parse_num("NODATA_value")
  if(is.na(nodata)) nodata <- -9999

  # remaining lines: one row per line, whitespace-separated
  body <- readLines(con)
  vals <- suppressWarnings(as.numeric(unlist(strsplit(trimws(body), "\\s+"))))
  if(length(vals) != ncols * nrows)
    stop(sprintf("ASC data length mismatch in %s: got %d, expected %d",
                 path, length(vals), ncols * nrows))
  m <- matrix(vals, nrow = nrows, ncol = ncols, byrow = TRUE)
  attr(m, "nodata") <- nodata
  m
}

#' @title Normalise a spatial matrix to unit sum on non-NODATA cells.
#' @description Sets NODATA cells to NA, divides remaining cells by their sum.
#'   Returns NA for the NODATA cells so the shared valid-mask logic in
#'   \code{fn.spatial_LL} can use \code{is.finite()} uniformly.
#' @keywords internal
#' @noRd
.normalise_map <- function(m, nodata = attr(m, "nodata")){
  if(!is.null(nodata)) m[m == nodata] <- NA_real_
  s <- sum(m, na.rm = TRUE)
  if(!is.finite(s) || s <= 0) return(m * NA_real_)
  m / s
}

#' @title Default GFISHER mod-code -> WFS MICE pool-code crosswalk.
#' @description gag ages 0-5+ = mod 20-25 = pc 4-9; red grouper ages 0-5+ =
#'   mod 26-31 = pc 10-15. Crosswalk is \code{pc = mod - 16}.
#' @return data.frame with columns \code{mod, pc, species}.
#' @export
fn.gfisher_default_crosswalk <- function(){
  data.frame(
    mod     = 20:31,
    pc      = 4:15,
    species = c(paste("gag",         0:4), "gag 5+",
                paste("red grouper", 0:4), "red grouper 5+"),
    stringsAsFactors = FALSE
  )
}

#' @title Load GFISHER maxN heatmaps as an observed spatial obs frame.
#' @description Reads the 12 gag + red-grouper maxN `.asc` files, normalises
#'   each to unit sum over non-NODATA cells, attaches the WFS MICE pool code
#'   via the crosswalk, and records the collection year window (\code{styear},
#'   \code{endyear}) alongside each row for year-matched comparison in the
#'   objective function.
#' @param gfisher_dir Path to the directory containing the maxN `.asc` files
#'   (default: WFS MICE canonical location).
#' @param styear,endyear Year window this obs represents. Defaults 2015-2024.
#'   Applied to every row unless overridden in the crosswalk.
#' @param crosswalk data.frame with columns \code{mod, pc, species}; see
#'   \code{fn.gfisher_default_crosswalk}. Optionally include a \code{path}
#'   column to point at a specific file; otherwise the file is located by
#'   pattern \code{"_mod<mod>_"}.
#' @param variable Character; portion of the filename identifying which
#'   variable (default \code{"maxn"}).
#' @return A data.frame with columns \code{mod, pc, species, path, styear,
#'   endyear} and a list-column \code{map_norm} (each element is a numeric
#'   matrix, NA on NODATA cells, unit-sum-normalised).
#' @export
fn.load_gfisher_maxn <- function(gfisher_dir = NULL,
                                 styear      = 2015,
                                 endyear     = 2024,
                                 crosswalk   = NULL,
                                 variable    = "maxn"){
  if(is.null(gfisher_dir))
    gfisher_dir <- "C:/Users/dchagaris/OneDrive - University of Florida/WFS Fisheries Ecosystem Modeling/WFS EwE/Ecospace/maps/GFISHER/5min/maxn"
  if(!dir.exists(gfisher_dir))
    stop("GFISHER directory not found: ", gfisher_dir)
  if(is.null(crosswalk)) crosswalk <- fn.gfisher_default_crosswalk()
  required <- c("mod", "pc", "species")
  if(!all(required %in% names(crosswalk)))
    stop("crosswalk must have columns: ", paste(required, collapse = ", "))

  # locate files if not provided per-row
  all_asc <- list.files(gfisher_dir, pattern = "\\.asc$", full.names = TRUE)
  if("path" %in% names(crosswalk)){
    paths <- crosswalk$path
  } else {
    paths <- vapply(crosswalk$mod, function(m){
      hits <- grep(sprintf("_mod%d_", m), basename(all_asc))
      if(length(hits) == 0) NA_character_ else all_asc[hits[1]]
    }, character(1))
  }
  missing <- which(is.na(paths))
  if(length(missing) > 0)
    warning("No GFISHER file found for mod code(s): ",
            paste(crosswalk$mod[missing], collapse = ", "))

  # load + normalise
  map_norm <- vector("list", nrow(crosswalk))
  for(i in seq_len(nrow(crosswalk))){
    if(is.na(paths[i])) next
    m <- .read_asc(paths[i])
    map_norm[[i]] <- .normalise_map(m)
  }

  data.frame(
    mod      = crosswalk$mod,
    pc       = crosswalk$pc,
    species  = crosswalk$species,
    path     = paths,
    styear   = if("styear"  %in% names(crosswalk)) crosswalk$styear  else styear,
    endyear  = if("endyear" %in% names(crosswalk)) crosswalk$endyear else endyear,
    variable = variable,
    stringsAsFactors = FALSE
  ) -> out
  out$map_norm <- map_norm
  drop_missing <- vapply(map_norm, is.null, logical(1))
  if(any(drop_missing)){
    message(sprintf("[fn.load_gfisher_maxn] dropping %d row(s) with no matching .asc file.",
                    sum(drop_missing)))
    out <- out[!drop_missing, , drop = FALSE]
  }
  out
}

#' @title Parse an Ecospace per-species map CSV.
#' @description The EcospaceMap<Var>-<species>.csv files stack all timesteps
#'   for one species in one file. Header block (\code{<HEADER software/>} ...
#'   \code{<HEADER end/>}) then repeated blocks of \code{blank / Step,N /
#'   Year,Y / <nrows> grid rows}. Returns a 3D array \code{[nrows, ncols,
#'   n_steps]} plus attributes \code{step, year}.
#' @param path Path to the per-species CSV.
#' @param nrows,ncols Expected grid dimensions. Defaults 66 x 78 (WFS MICE).
#' @return 3D numeric array with attributes \code{step} (integer vector),
#'   \code{year} (numeric vector). NODATA -9999 preserved (mask downstream).
#' @keywords internal
#' @noRd
.parse_ecospace_map_csv <- function(path, nrows = 66L, ncols = 78L){
  ln <- readLines(path)
  end_idx <- grep("<HEADER end/>", ln)[1]
  if(is.na(end_idx))
    stop("Could not locate <HEADER end/> in ", path)

  # skip header + Variable + Contaminant lines
  body <- ln[(end_idx + 1L):length(ln)]

  step_lines <- grep("^Step,", body)
  n_steps    <- length(step_lines)
  if(n_steps == 0L) stop("No Step blocks in ", path)

  steps  <- integer(n_steps)
  years  <- numeric(n_steps)
  arr    <- array(NA_real_, dim = c(nrows, ncols, n_steps))

  for(k in seq_len(n_steps)){
    idx     <- step_lines[k]
    steps[k] <- as.integer(sub("^Step,", "", body[idx]))
    yr_line <- body[idx + 1L]
    if(!grepl("^Year,", yr_line))
      stop(sprintf("Expected 'Year,' after Step block at line %d of %s", idx, path))
    years[k] <- as.numeric(sub("^Year,", "", yr_line))

    # next nrows lines are the grid
    grid_lines <- body[(idx + 2L):(idx + 1L + nrows)]
    vals <- suppressWarnings(as.numeric(unlist(strsplit(grid_lines, ","))))
    if(length(vals) != nrows * ncols)
      stop(sprintf("Grid dim mismatch at step %d of %s (got %d, expected %d)",
                   steps[k], path, length(vals), nrows * ncols))
    arr[, , k] <- matrix(vals, nrow = nrows, ncol = ncols, byrow = TRUE)
  }
  attr(arr, "step") <- steps
  attr(arr, "year") <- years
  arr
}

#' @title Read predicted Ecospace spatial biomass averaged over target years.
#' @description For each requested pool code, locates the matching
#'   \code{EcospaceMap<variable>-<species>.csv} (species inferred via
#'   \code{group.names[pc]}), parses the stacked timesteps, extracts steps
#'   whose \code{Year} attribute (fractional-year of model time; e.g. 30.083
#'   for step 361 of a monthly model starting at year 0) falls in the target
#'   year window, and averages the corresponding grids.
#' @param dir.pred Path to an Ecospace output folder (parent of the
#'   \code{csv/} subdir).
#' @param pool_codes Integer vector of WFS MICE pool codes to read.
#' @param group.names Character vector indexed by pool code.
#' @param styear,endyear Calendar year window (inclusive) to average over.
#'   Converted to model-relative years via \code{model_styear}: model years
#'   run from \code{model_styear} at Ecospace step 1.
#' @param model_styear Calendar year corresponding to Ecospace step 1.
#'   Default 1985 (WFS MICE convention).
#' @param variable Character; part of the filename (default \code{"Biomass"}).
#' @param nrows,ncols Grid dimensions. Defaults 66 x 78.
#' @return Named list of numeric matrices (one per pool code); NA cells are
#'   NODATA. NULL entry for any pool code whose file was not found.
#' @export
fn.ecospace_spatial_pred <- function(dir.pred, pool_codes, group.names,
                                     styear       = 2015,
                                     endyear      = 2024,
                                     model_styear = 1985,
                                     variable     = "Biomass",
                                     nrows        = 66L,
                                     ncols        = 78L){
  # GUI-run outputs put per-species map CSVs in a csv/ subfolder;
  # CLI-run (EcoSpace Console) writes them at the top level of dir.pred.
  # Support both.
  csv_dir <- if(dir.exists(file.path(dir.pred, "csv"))) file.path(dir.pred, "csv") else dir.pred
  if(length(list.files(csv_dir, pattern = "^EcospaceMap.*\\.csv$")) == 0)
    stop("No EcospaceMap<Var>-<species>.csv files found in ", csv_dir,
         " (nor a csv/ subfolder under ", dir.pred, ")")

  # convert target calendar years -> model fractional-year window.
  # Ecospace's CSV Year field is fractional years from model start (step*dt).
  yr_min <- styear   - model_styear
  yr_max <- endyear  - model_styear + 1  # inclusive of endyear months

  out <- vector("list", length(pool_codes))
  names(out) <- as.character(pool_codes)

  all_csv <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

  for(i in seq_along(pool_codes)){
    pc <- pool_codes[i]
    if(is.na(pc) || pc < 1 || pc > length(group.names)) next
    # Ecospace writes group names to the CSV filename with SPACES; R4EwE
    # group.names is stored with UNDERSCORES. Convert on the fly.
    species <- gsub("_", " ", group.names[pc])

    pat <- paste0("^EcospaceMap", variable, "-",
                  gsub("\\+", "\\\\+", species),
                  "\\.csv$")
    hits <- grep(pat, basename(all_csv))
    if(length(hits) == 0){
      warning(sprintf("No spatial CSV for pool code %d (%s) in %s",
                      pc, species, csv_dir))
      next
    }
    arr3 <- .parse_ecospace_map_csv(all_csv[hits[1]],
                                    nrows = nrows, ncols = ncols)
    yrs  <- attr(arr3, "year")
    keep <- which(yrs >= yr_min & yrs <= yr_max)
    if(length(keep) == 0){
      warning(sprintf("No timesteps in year window [%g, %g] for pool code %d",
                      yr_min, yr_max, pc))
      next
    }
    m <- apply(arr3[, , keep, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
    m[m == -9999] <- NA_real_
    out[[i]] <- m
  }
  out
}

#' @title Spatial likelihood component for fn.ecospace_objfxn.
#' @description Per row of \code{spatial.obs}: builds predicted map for that
#'   pool code, averages over the row's own \code{styear}..\code{endyear}
#'   window, unit-sum normalises both obs and pred on the shared valid mask
#'   (non-NA in both), and accumulates cross-entropy
#'   \eqn{-\sum \tilde o \log \tilde p} (multinomial NLL kernel). Returns
#'   weighted total plus a per-species named vector.
#' @param dir.pred Path to an Ecospace output folder (parent of \code{csv/}).
#' @param spatial.obs Output of \code{fn.load_gfisher_maxn} (or similar), with
#'   columns \code{pc, species, styear, endyear} and list-col \code{map_norm}.
#' @param spatial.weight Numeric scalar multiplier applied to the summed cross
#'   entropy (interpretable as a prior precision on the spatial fit).
#' @param group.names Character vector indexed by pool code (needed by
#'   \code{fn.ecospace_spatial_pred}).
#' @param model_styear Calendar year corresponding to Ecospace step 1
#'   (default 1985).
#' @return List with:
#'   \code{total}  = weighted total cross-entropy (scalar),
#'   \code{per_species} = named numeric vector, one weighted cross-entropy per
#'     row of \code{spatial.obs}. NA entries indicate rows that could not be
#'     scored (missing prediction, empty window, all-NODATA overlap).
#' @export
fn.spatial_LL <- function(dir.pred, spatial.obs,
                          spatial.weight = 1,
                          group.names,
                          model_styear   = 1985){
  if(is.null(spatial.obs) || nrow(spatial.obs) == 0)
    return(list(total = 0, per_species = numeric(0)))

  pool_codes <- spatial.obs$pc

  # Group rows by (styear, endyear) so we read each species' CSV only once
  # per distinct window. spatial.obs will typically have a single common
  # window (2015-2024) so this collapses to one read per species.
  windows <- unique(spatial.obs[, c("styear", "endyear")])
  pred_maps <- vector("list", length(pool_codes))
  names(pred_maps) <- as.character(pool_codes)

  for(w in seq_len(nrow(windows))){
    rows <- which(spatial.obs$styear == windows$styear[w] &
                  spatial.obs$endyear == windows$endyear[w])
    pcs_w <- pool_codes[rows]
    got <- fn.ecospace_spatial_pred(dir.pred, pcs_w, group.names,
                                    styear = windows$styear[w],
                                    endyear = windows$endyear[w],
                                    model_styear = model_styear)
    for(k in seq_along(rows)){
      pred_maps[[as.character(pcs_w[k])]] <- got[[as.character(pcs_w[k])]]
    }
  }

  per <- setNames(rep(NA_real_, nrow(spatial.obs)), spatial.obs$species)
  for(i in seq_len(nrow(spatial.obs))){
    obs_m  <- spatial.obs$map_norm[[i]]
    pred_m <- pred_maps[[as.character(spatial.obs$pc[i])]]
    if(is.null(obs_m) || is.null(pred_m)) next
    if(!identical(dim(obs_m), dim(pred_m))){
      warning(sprintf("Dim mismatch obs %s vs pred (%s vs %s) - skipping.",
                      spatial.obs$species[i],
                      paste(dim(obs_m),  collapse = "x"),
                      paste(dim(pred_m), collapse = "x")))
      next
    }
    # shared valid mask (both finite)
    valid <- is.finite(obs_m) & is.finite(pred_m)
    if(!any(valid)) next

    # normalise on the SHARED mask (so both integrate to 1 over identical cells)
    obs_v  <- obs_m[valid];  obs_v  <- obs_v  / sum(obs_v,  na.rm = TRUE)
    pred_v <- pred_m[valid]; sp     <- sum(pred_v, na.rm = TRUE)
    if(!is.finite(sp) || sp <= 0) next
    pred_v <- pred_v / sp

    # Cross-entropy (multinomial NLL kernel). eps floor on pred prevents
    # log(0) blow-up where obs has probability mass in a cell pred zeroed out.
    eps    <- .Machine$double.eps^0.5
    per[i] <- -sum(obs_v * log(pmax(pred_v, eps)), na.rm = TRUE)
  }

  total <- spatial.weight * sum(per, na.rm = TRUE)
  list(total = total, per_species = spatial.weight * per)
}

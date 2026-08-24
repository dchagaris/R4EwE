#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecosim timeseries type-code lookup--------------------------------------------------------------------------------
#' @title Load the Ecosim timeseries type-code lookup table.
#' @description Returns a data frame describing every Ecosim timeseries type code (biomass, catch,
#'   mortality, effort, discards, prices, forcing series, etc.), the pool-code role required by each
#'   type, and flags indicating whether the series is in absolute units and whether it is a reference
#'   (calibration) series. The lookup CSV ships with the package at
#'   \code{inst/extdata/time_series_codes.csv}.
#' @return A data frame with columns: \code{code} (integer type code), \code{name} (canonical type
#'   name, e.g. "BiomassRel"), \code{absolute} (logical), \code{reference} (logical), and
#'   \code{pool_role} (text describing what the pool code(s) refer to).
#' @export
fn.get_ts_type_table <- function(){
  f <- system.file("extdata", "time_series_codes.csv", package = "R4EwE")
  if(!nzchar(f) || !file.exists(f)){
    # devtools::load_all fallback: look under the package source tree
    pkg <- tryCatch(find.package("R4EwE"), error = function(e) NULL)
    if(!is.null(pkg)) f <- file.path(pkg, "inst", "extdata", "time_series_codes.csv")
  }
  if(!file.exists(f)) stop("Could not locate time_series_codes.csv lookup file.")
  raw <- utils::read.csv(f, stringsAsFactors = FALSE)  # default check.names dedupes "Type" -> "Type.1"
  # Required columns in the CSV (after read.csv name-cleaning):
  #   Title, Weight, Pool.code, Pool.code.2, Type (text like "-18 or FixedCostRel"),
  #   Type.code (numeric), Type.1 (canonical name e.g. "FixedCostRel"), absolute, reference
  out <- data.frame(
    code      = suppressWarnings(as.integer(raw$Type.code)),
    name      = trimws(raw$Type.1),
    absolute  = toupper(trimws(as.character(raw$absolute)))  == "TRUE",
    reference = toupper(trimws(as.character(raw$reference))) == "TRUE",
    pool1     = trimws(as.character(raw$Pool.code)),
    pool2     = trimws(as.character(raw$Pool.code.2)),
    stringsAsFactors = FALSE
  )
  out$pool_role <- ifelse(nzchar(out$pool2),
                          paste(out$pool1, out$pool2, sep = " + "),
                          out$pool1)
  out <- out[!is.na(out$code), c("code", "name", "absolute", "reference", "pool_role")]
  out
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecosim observed timeseries---------------------------------------------------------------------------------------
#' @title Read an Ecosim timeseries file.
#' @description Reads an Ecosim timeseries CSV and parses every series type defined in the type-code
#'   lookup (see \code{\link{fn.get_ts_type_table}}). Handles 4-row headers (Title, Weight, Pool
#'   code, Type), 5-row headers (with Pool code 1 / Pool code 2), and the extended 6-row header
#'   used for Ecospace runs (adds a "Region" row identifying the Ecospace region per series). The
#'   "Type" and optional "Region" rows are auto-detected; the data block is taken to start after
#'   the last header row found.
#' @param filename Path to the Ecosim timeseries CSV.
#' @param types Optional subset filter. Either a numeric vector of type codes (e.g. \code{c(0, 1, 6)})
#'   or a character vector of canonical type names (e.g. \code{c("BiomassRel", "Catches")}). Two
#'   shortcuts are also accepted: \code{"reference"} (only calibration/reference series) or
#'   \code{"forcing"} (only forcing series). Default \code{NULL} returns everything.
#' @param long Logical. If \code{TRUE}, the returned list also includes a tidy long-format data frame
#'   \code{ts.long} with columns Year, Series, TypeCode, Type, Poolcode, Poolcode2, Weight, Absolute,
#'   Reference, Value. Defaults to \code{FALSE}.
#' @return A list with:
#'   \itemize{
#'     \item \code{ts.head} -- enriched header table for all (filtered) series, including the
#'       canonical type name and absolute / reference flags.
#'     \item \code{ts} -- wide data frame of values (rows = years, columns = series), with column
#'       names set to the cleaned series titles.
#'     \item \code{by.type} -- named list keyed by canonical type name; each entry is
#'       \code{list(head = ..., data = ...)} containing only the series of that type present in the
#'       file.
#'     \item \code{obsB.head}, \code{obsB}, \code{obsC.head}, \code{obsC} -- backwards-compatible
#'       biomass (codes 0, 1) and catch (codes 6, 61, -6) views. \code{obsB.head}/\code{obsC.head}
#'       retain the legacy 4-column layout (Title, Weight, Poolcode, Type) so existing callers in
#'       \code{plotting.R} and \code{Ecospace_objective_functions.R} keep working.
#'     \item \code{ts.long} -- only if \code{long = TRUE}.
#'   }
#' @examples
#' \dontrun{
#' ts <- fn.read_ecosim_timeseries("ts_mice_v4_discards.csv")
#' names(ts$by.type)
#' ts$by.type$BiomassRel$data
#' ts_ref <- fn.read_ecosim_timeseries("ts.csv", types = "reference", long = TRUE)
#' }
#' @export
fn.read_ecosim_timeseries <- function(filename, types = NULL, long = FALSE){

  # ---- 1. detect number of header rows by finding the "Type" row and (optional) "Region" row ----
  raw_lines <- readLines(filename, n = 20L)
  first_cell <- vapply(strsplit(raw_lines, ","), function(x){
    if(length(x) == 0) "" else gsub('"', '', trimws(x[1]))
  }, character(1))
  type_row <- which(tolower(first_cell) %in% c("type", "type code", "typecode"))
  if(length(type_row) == 0)
    stop("Could not locate the 'Type' header row in: ", filename)
  type_row   <- type_row[1]
  region_row <- which(tolower(first_cell) == "region")
  region_row <- if(length(region_row) > 0) region_row[1] else NA_integer_
  nhead <- max(type_row, region_row, na.rm = TRUE)

  # ---- 2. read header block, transpose so each row = one series ----
  hdr.raw <- utils::read.csv(filename, header = FALSE, nrows = nhead, stringsAsFactors = FALSE)
  row_labels <- as.character(hdr.raw[, 1])
  hdr <- as.data.frame(t(hdr.raw[, -1, drop = FALSE]), stringsAsFactors = FALSE)
  names(hdr) <- row_labels
  rownames(hdr) <- NULL

  # standardize header column names so downstream code can rely on them
  std_names <- vapply(names(hdr), function(x){
    y <- gsub("\\s+", "", x)            # "Pool code 1" -> "Poolcode1", "Pool code" -> "Poolcode"
    y <- sub("code1$", "code", y, ignore.case = TRUE)  # "Poolcode1" -> "Poolcode"
    y <- sub("code2$", "code2", y, ignore.case = TRUE) # "Poolcode2" stays as "Poolcode2"
    y
  }, character(1))
  names(hdr) <- std_names

  # preserve raw character strings for potentially comma-separated cells
  # BEFORE numeric coercion drops them to NA. Enables multi-pool / multi-region
  # type-0/1 series (v5+ format).
  raw_poolcode <- if("Poolcode" %in% names(hdr)) as.character(hdr$Poolcode) else NULL
  raw_region   <- if("Region"   %in% names(hdr)) as.character(hdr$Region)   else NULL

  # coerce non-Title columns to numeric
  for(j in setdiff(names(hdr), "Title")){
    hdr[[j]] <- suppressWarnings(as.numeric(hdr[[j]]))
  }
  if(!"Poolcode2" %in% names(hdr)) hdr$Poolcode2 <- NA_real_
  if(!"Region"    %in% names(hdr)) hdr$Region    <- NA_integer_

  # parse comma-separated pool codes / regions into list-columns. Scalar
  # Poolcode / Region are set from the first element (back-compat for
  # type-12/13/61 paths). Empty cells -> integer(0).
  .parse_int_list <- function(x){
    lapply(x, function(s){
      if(is.null(s) || is.na(s) || !nzchar(trimws(s))) return(integer(0))
      v <- suppressWarnings(as.integer(strsplit(trimws(s), "\\s*,\\s*")[[1]]))
      v[!is.na(v)]
    })
  }
  hdr$Poolcodes <- if(!is.null(raw_poolcode)) .parse_int_list(raw_poolcode)
                   else rep(list(integer(0)), nrow(hdr))
  hdr$Regions   <- if(!is.null(raw_region))   .parse_int_list(raw_region)
                   else rep(list(integer(0)), nrow(hdr))

  # rewrite scalar Poolcode / Region from first parsed element (so legacy
  # code paths still see a single value; multi-cell entries kept only in the
  # list-columns).
  hdr$Poolcode <- vapply(hdr$Poolcodes,
                         function(v) if(length(v)) as.numeric(v[1]) else NA_real_,
                         numeric(1))
  hdr$Region   <- vapply(hdr$Regions,
                         function(v) if(length(v)) as.integer(v[1]) else NA_integer_,
                         integer(1))

  # report multi-value counts (v5+ diagnostic)
  n_multi_pc <- sum(lengths(hdr$Poolcodes) > 1)
  n_multi_r  <- sum(lengths(hdr$Regions)   > 1)
  if(n_multi_pc + n_multi_r > 0)
    message(sprintf("[fn.read_ecosim_timeseries] multi-value entries: %d series with multi Pool code, %d series with multi Region.",
                    n_multi_pc, n_multi_r))

  # ---- 3. read data block ----
  dat <- utils::read.csv(filename, header = FALSE, skip = nhead, stringsAsFactors = FALSE)
  years <- dat[, 1]
  dat <- dat[, -1, drop = FALSE]
  rownames(dat) <- years
  if(nrow(hdr) != ncol(dat))
    warning("HEADER AND TIMESERIES DIMENSION DO NOT MATCH (header rows=", nrow(hdr),
            ", data cols=", ncol(dat), ")")
  colnames(dat) <- gsub(" ", "_", hdr$Title)

  # ---- 4. enrich header with type-code lookup ----
  codes <- fn.get_ts_type_table()
  m <- match(hdr$Type, codes$code)
  hdr$Type.name <- codes$name[m]
  hdr$Absolute  <- codes$absolute[m]
  hdr$Reference <- codes$reference[m]
  hdr$PoolRole  <- codes$pool_role[m]
  unknown <- which(is.na(m))
  if(length(unknown))
    warning("Unknown type code(s) encountered: ",
            paste(unique(hdr$Type[unknown]), collapse = ", "),
            " (series: ", paste(hdr$Title[unknown], collapse = "; "), ")")

  # ---- 5. resolve types filter ----
  keep <- seq_len(nrow(hdr))
  if(!is.null(types)){
    if(is.character(types) && length(types) == 1 && tolower(types) == "reference"){
      keep <- which(isTRUE_vec(hdr$Reference))
    } else if(is.character(types) && length(types) == 1 && tolower(types) == "forcing"){
      keep <- which(!isTRUE_vec(hdr$Reference))
    } else if(is.numeric(types)){
      keep <- which(hdr$Type %in% types)
    } else if(is.character(types)){
      keep <- which(hdr$Type.name %in% types)
    } else {
      stop("'types' must be numeric (type codes) or character (type names or 'reference'/'forcing').")
    }
  }

  # ---- 6. by.type list (one entry per type name present in the filtered selection) ----
  by.type <- list()
  for(tn in unique(stats::na.omit(hdr$Type.name[keep]))){
    j <- intersect(keep, which(hdr$Type.name == tn))
    by.type[[tn]] <- list(head = hdr[j, , drop = FALSE],
                          data = dat[, j, drop = FALSE])
  }

  # ---- 7. backwards-compatible legacy views: obsB / obsC with 4-col head ----
  legacy.head <- data.frame(
    Title    = hdr$Title,
    Weight   = hdr$Weight,
    Poolcode = hdr$Poolcode,
    Type     = hdr$Type,
    stringsAsFactors = FALSE
  )
  bIdx <- which(hdr$Type %in% c(0, 1))
  cIdx <- which(hdr$Type %in% c(6, 61, -6))
  obsB.head <- legacy.head[bIdx, , drop = FALSE]
  obsB      <- dat[, bIdx, drop = FALSE]
  obsC.head <- legacy.head[cIdx, , drop = FALSE]
  obsC      <- dat[, cIdx, drop = FALSE]
  rownames(obsB.head) <- gsub(" ", "_", obsB.head$Title)
  rownames(obsC.head) <- gsub(" ", "_", obsC.head$Title)

  # ---- 8. assemble return ----
  ret <- list(
    ts.head   = hdr[keep, , drop = FALSE],
    ts        = dat[, keep, drop = FALSE],
    by.type   = by.type,
    obsB.head = obsB.head, obsB = obsB,
    obsC.head = obsC.head, obsC = obsC
  )

  if(long){
    sh <- hdr[keep, , drop = FALSE]
    sd <- dat[, keep, drop = FALSE]
    yrs <- suppressWarnings(as.numeric(rownames(sd)))
    ret$ts.long <- do.call(rbind, lapply(seq_len(ncol(sd)), function(j){
      data.frame(
        Year      = yrs,
        Series    = sh$Title[j],
        TypeCode  = sh$Type[j],
        Type      = sh$Type.name[j],
        Poolcode  = sh$Poolcode[j],
        Poolcode2 = sh$Poolcode2[j],
        Region    = sh$Region[j],
        Weight    = sh$Weight[j],
        Absolute  = sh$Absolute[j],
        Reference = sh$Reference[j],
        Value     = sd[, j],
        stringsAsFactors = FALSE
      )
    }))
  }

  return(ret)
}

# small internal helper: TRUE for TRUE, FALSE for NA or FALSE
#' @keywords internal
#' @noRd
isTRUE_vec <- function(x) !is.na(x) & x

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecosim output to arrays---------------------------------------------------------------------------------  
#' @title Read Ecosim predicted biomass timeseries output.
#' @description Reads biomass timeseries output csv into an array.
#' @param dir.out Output directory.  All .csv biomass files created with Ecospace naming conventions 
#' in this directory will be read.  Can be a vector of directories.
#' @param timestep Read 'annual' or 'monthly' data.
#' @return An array of biomass output with dimensions (nyrs,ngrps,nregions,nruns)
#' #examples
#' #example code:
#' predB <- fn.ecospace_predB_ts2array (dir.out="./output/", timestep='annual',n.reg=0)
#' @export
fn.ecosim_predB_ts2array = function(dir.out=dir.pred, timestep='annual'){
  
  if(timestep=='annual') files.bio = list.files(dir.out,pattern="^biomass_annual.csv",recursive=T,full.names = T)
  if(timestep=='monthly') files.bio = list.files(dir.out,pattern="^biomass_monthly.csv",recursive=T,full.names = T)
  
  x = list.files(dir.out,pattern="^biomass_annual.csv",recursive=T,full.names = T)[1]
  nskip = which(tolower(substr(readLines(x),1,4))=='year')-1  
  styear <- as.numeric(read.csv(x, as.is=T, skip=nskip)[1,1])
  
  bio = lapply(files.bio,FUN=function(x){
      if(timestep=='annual') nskip = which(tolower(substr(readLines(x),1,4))=='year')-1
      if(timestep=='monthly') nskip = which(tolower(substr(readLines(x),1,8))=='timeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1)
    })
    
  bio.array = array(dim=c(dim(bio[[1]])[1],dim(bio[[1]])[2],length(files.bio)),
                    dimnames=list(rownames(bio[[1]]),names(bio[[1]]),basename(dirname(files.bio))))
  for(r in 1:length(bio)){
    tmp = as.matrix(bio[[r]])
    bio.array[,,r] <- tmp
  }
  
  if(timestep=='monthly') dimnames(bio.array)[[1]] <- length(seq(styear,(styear+dim(bio.array)[1]-1/12),1/12))
  return(bio.array)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecosim output to arrays---------------------------------------------------------------------------------  
#' @title Read Ecosim predicted catch timeseries output.
#' @description Reads catch timeseries output csv, sums over fleets for each species, and puts into 
#' an array.
#' @param dir.out Output directory.  All .csv catch files created with Ecospace naming conventions 
#' in this directory will be read. Can be a vector of directories.
#' @param timestep Read 'annual' or 'monthly' data.
#' @return An array of catch output with dimensions (nyrs,ngrps,nregions,nruns)
#' #examples
#' #example code:
#' predB <- fn.ecospace_predB_ts2array (dir.out="./output/", timestep='annual',n.reg=0)
#' @export
fn.ecosim_predC_ts2array = function(dir.out=dir.out, timestep='annual'){
    if(timestep=='annual') files.cat = list.files(dir.out,pattern="^catch_annual.csv",recursive=T,full.names = T)
    if(timestep=='monthly') files.cat = list.files(dir.out,pattern="^catch_monthly.csv",recursive=T,full.names = T)
    
    x = list.files(dir.out,pattern="^catch_annual.csv",recursive=T,full.names = T)[1]
    nskip = which(tolower(substr(readLines(x),1,4))=='year')-1  
    styear <- as.numeric(read.csv(x, as.is=T, skip=nskip)[1,1])
    
    cat = lapply(files.cat,FUN=function(x){
      if(timestep=='annual') nskip = which(tolower(substr(readLines(x),1,4))=='year')-1
      if(timestep=='monthly') nskip = which(tolower(substr(readLines(x),1,8))=='timeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1, check.names=F)
    })
    cat.array = array(dim=c(dim(cat[[1]])[1],dim(cat[[1]])[2],length(files.cat)),
                      dimnames=list(rownames(cat[[1]]),names(cat[[1]]),basename(dirname(files.cat))))
    for(r in 1:length(cat)){
      tmp = as.matrix(cat[[r]])
      cat.array[,,r] <- tmp
    }
  
  #if(timestep=='annual') dimnames(cat.array)[[1]] <- styear:(styear+dim(cat.array)[1]-1)
  if(timestep=='monthly') dimnames(cat.array)[[1]] <- length(seq(styear,(styear+dim(cat.array)[1]-1/12),1/12))
  return(cat.array)
}#eof




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecospace output to arrays---------------------------------------------------------------------------------  
#' @title Read Ecospace predicted biomass timeseries output.
#' @description Reads biomass timeseries output csv into an array.
#' @param dir.out Output directory.  All .csv biomass files created with Ecospace naming conventions 
#' in this directory will be read.
#' @param timestep Read 'annual' or 'monthly' data.
#' @param n.reg If output for multiple regions, how many?  Default to 0.
#' @return An array of biomass output with dimensions (nyrs,ngrps,nregions,nruns)
#' #examples
#' #example code:
#' predB <- fn.ecospace_predB_ts2array (dir.out="./output/", timestep='annual',n.reg=0)
#' @export

fn.ecospace_predB_ts2array = function(dir.out=dir.pred, timestep='annual',n.reg=0){
  
  #dir.out = "C:/NWACS MICE/GA output 2026-01-06/GA_Run_20260106_114454/run_00abeb6e03a8897940ee4aa45b593526"
  if(n.reg==0){
    if(timestep=='annual') files.bio = list.files(dir.out,pattern="Ecospace_Annual_Average_Biomass.csv",recursive=T,full.names = T)
    if(timestep=='monthly') files.bio = list.files(dir.out,pattern="Ecospace_Average_Biomass.csv",recursive=T,full.names = T)
    
    file.bio <- files.bio[grep("Region_0",basename(files.bio))]
    
    bio = lapply(files.bio,FUN=function(x){
      if(timestep=='annual') nskip = which(substr(readLines(x),1,4)=='Year')-1
      if(timestep=='monthly') nskip = which(substr(readLines(x),1,8)=='TimeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1)
    })
    
    bio.array = array(dim=c(dim(bio[[1]])[1],dim(bio[[1]])[2],length(files.bio)),
                      dimnames=list(rownames(bio[[1]]),names(bio[[1]]),basename(dirname(files.bio))))
    for(r in 1:length(bio)){
      tmp = as.matrix(bio[[r]])
      bio.array[,,r] <- tmp
    }
  }
  
  if(n.reg>0){
    if(timestep=='annual') files.bio = list.files(dir.out,pattern="^Ecospace_Annual_Average_Region.*\\Biomass.csv$",recursive=T,full.names = T)
    if(timestep=='monthly')  files.bio = list.files(dir.out,pattern="^Ecospace_Average_Region.*\\Biomass.csv$",recursive=T,full.names = T)
    
    bio = lapply(files.bio,FUN=function(x){
      if(timestep=='annual') nskip = which(substr(readLines(x),1,4)=='Year')-1
      if(timestep=='monthly') nskip = which(substr(readLines(x),1,8)=='TimeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1)
    })
    
    bio.array = array(dim=c(dim(bio[[1]])[1],dim(bio[[1]])[2],n.reg+1,length(unique(dirname(files.bio)))),
                      dimnames=list(rownames(bio[[1]]),names(bio[[1]]),paste0('reg',0:n.reg),unique(basename(dirname(files.bio)))))
    for(r in 1:length(bio)){
      if(timestep=='annual') reg.idx=as.numeric(substr(basename(files.bio[r]),32,32))+1
      if(timestep=='monthly') reg.idx=as.numeric(substr(basename(files.bio[r]),25,25))+1
      run.idx=match(basename(dirname(files.bio[r])),dimnames(bio.array)[[4]])
      tmp = as.matrix(bio[[r]])
      bio.array[,,reg.idx,run.idx] <- tmp
    }
  }
  if(timestep=='annual') dimnames(bio.array)[[1]] <- styear:(styear+dim(bio.array)[1]-1)
  if(timestep=='monthly') dimnames(bio.array)[[1]] <- length(seq(styear,(styear+dim(bio.array)[1]-1/12),1/12))
  return(bio.array)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecospace output to arrays---------------------------------------------------------------------------------  
#' @title Read Ecospace predicted catch timeseries output.
#' @description Reads catch timeseries output csv, sums over fleets for each species, and puts into 
#' an array.
#' @param dir.out Output directory.  All .csv catch files created with Ecospace naming conventions 
#' in this directory will be read.
#' @param timestep Read 'annual' or 'monthly' data.
#' @param n.reg If output for multiple regions, how many?  Default to 0.
#' @return An array of catch output with dimensions (nyrs,ngrps,nregions,nruns)
#' #examples
#' #example code:
#' predB <- fn.ecospace_predB_ts2array (dir.out="./output/", timestep='annual',n.reg=0)
#' @export
fn.ecospace_predC_ts2array = function(dir.out=dir.pred, timestep='annual',n.reg=0){
  if(n.reg==0){
    if(timestep=='annual') files.cat = list.files(dir.out,pattern="Ecospace_Annual_Average_Catch.csv",recursive=T,full.names = T)
    if(timestep=='monthly') files.cat = list.files(dir.out,pattern="Ecospace_Average_Catch.csv",recursive=T,full.names = T)
    
    file.cat <- files.cat[grep("Region_0",basename(files.cat))]
    
    cat = lapply(files.cat,FUN=function(x){
      if(timestep=='annual') nskip = which(substr(readLines(x),1,4)=='Year')-1
      if(timestep=='monthly') nskip = which(substr(readLines(x),1,8)=='TimeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1, check.names=F)
    })
    sppnames = sapply(strsplit(names(cat[[1]]),"\\|"),tail,1)
    cat.array = array(dim=c(dim(cat[[1]])[1],dim(cat[[1]])[2],length(files.cat)),
                      dimnames=list(rownames(cat[[1]]),names(cat[[1]]),basename(dirname(files.cat))))
    for(r in 1:length(cat)){
      tmp = as.matrix(cat[[r]])
      cat.array[,,r] <- tmp
    }
  }
  
  if(n.reg>0){
    if(timestep=='annual') files.cat = list.files(dir.out,pattern="^Ecospace_Annual_Average_Region.*\\Catch.csv$",recursive=T,full.names = T)
    if(timestep=='monthly')  files.cat = list.files(dir.out,pattern="^Ecospace_Average_Region.*\\Catch.csv$",recursive=T,full.names = T)
    
    cat = lapply(files.cat,FUN=function(x){
      if(timestep=='annual') nskip = which(substr(readLines(x),1,4)=='Year')-1
      if(timestep=='monthly') nskip = which(substr(readLines(x),1,8)=='TimeStep')-1
      read.csv(x, as.is=T, skip=nskip, row.names=1, check.names=F)
    })
    
    cat.array = array(dim=c(dim(cat[[1]])[1],dim(cat[[1]])[2],n.reg+1,length(unique(dirname(files.cat)))),
                      dimnames=list(rownames(cat[[1]]),names(cat[[1]]),paste0('reg',0:n.reg),unique(basename(dirname(files.cat)))))
    for(r in 1:length(cat)){
      if(timestep=='annual') reg.idx=as.numeric(substr(basename(files.cat[r]),32,32))+1
      if(timestep=='monthly') reg.idx=as.numeric(substr(basename(files.cat[r]),25,25))+1
      run.idx=match(basename(dirname(files.cat[r])),dimnames(cat.array)[[4]])
      tmp = as.matrix(cat[[r]])
      cat.array[,,reg.idx,run.idx] <- tmp
    }
  }
  if(timestep=='annual') dimnames(cat.array)[[1]] <- styear:(styear+dim(cat.array)[1]-1)
  if(timestep=='monthly') dimnames(cat.array)[[1]] <- length(seq(styear,(styear+dim(cat.array)[1]-1/12),1/12))
  return(cat.array)
}#eof

#aggregate catch by group
#' @keywords internal
#' @noRd
fn.agg_catch_by_group <- function(predC){
  print(dimnames(predC)[[2]])
  grps = sapply(strsplit(dimnames(predC)[[2]],"\\|"),tail,1)
  print(grps)
  length(unique(grps))
  
  #now aggregate columns in predC.agg to grps
  predC.agg <- apply(predC, c(1,3), function(v) {tapply(v, grps, sum)})
  predC.agg <- aperm(predC.agg, c(2,1,3))
  dimnames(predC.agg)
  return(predC.agg)
}#eof

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
#' @title Read Ecospace output ascii files into a raster stack.
#' @description This function reads a subset of asc grids produced by an Ecospace run into a raster 
#' stack.  Requires the terra package for faster processing.
#' @param dir.out The run output folder containing.  It should contain a folder named 'asc', which 
#' holds all the output maps.
#' @param do.vars The model output variables to read. Must match the naming convention of the output
#' files, all  lowercase.  Current options are 'biomass','catch', and 'effort'.
#' @param do.grps The model groups for which to read data.
#' @param do.fleets The fleet numbers for which to read effort data.
#' @param do.months The months to read output data, only used if monthly maps are output.
#' @param do.years Years to read output data, will combine with do.months
#' @param annualmaps Import annual (TRUE, default) or monthly (FALSE) maps from Ecospace.
#' @return A SpatRaster stack from the terra package.
#' @examples
#' # example code:
#' \dontrun{predB.maps <- fn.ecospace_ascii2stack(dir.out=dir.pred, do.vars=c('biomass'), do.grps=c(3,5,6,8,10,12), do.fleets=1:6, do.months=1, do.years=2010:2019, annualmaps=TRUE)}
#' @export
fn.ecospace_ascii2stack <- function(dir.out=dir.pred,   
                                    do.vars=c('biomass'), 
                                    do.grps=1:12, 
                                    do.fleets=1:6, 
                                    do.months=1, 
                                    do.years=2010:2019,
                                    annualmaps=TRUE){
  # dir.out=dir.pred
  # do.vars=c('biomass')
  # do.grps=1:12
  # do.fleets=1:6
  # do.months=1
  # do.years=2010:2019
  # annualmaps=TRUE
  files.ascii <- list.files(file.path(dir.out,'asc'),pattern=".asc$", recursive=T, full.names=T)
  grpsplit = strsplit(basename(files.ascii),"-")
  grpnames = character()
  for(g in 1:length(grpsplit)){
    #g=1326
    name.g = grpsplit[[g]]
    name.g = name.g[-c(1,length(name.g))]
    name.g = paste(name.g,collapse="-")
    grpnames = c(grpnames,name.g)
  }
  unique(grpnames)
  timestep = substr(basename(files.ascii),nchar(basename(files.ascii))-8,nchar(basename(files.ascii))-4)
  time = styear + (as.numeric(timestep)-1)/12
  #output map metadata----
  ascii.df = data.frame(type=ifelse(grepl('Biomass',basename(files.ascii)),'biomass',
                                    ifelse(grepl('Catch',basename(files.ascii)),'catch',
                                           ifelse(grepl('Effort',basename(files.ascii)),'effort','other'))),
                        group=grpnames,
                        group.num = match(gsub(" ","_",grpnames),df.names$group.names),
                        fleet.num = match(grpnames,df.names$fleet.names),
                        timestep =timestep,
                        time = time,
                        year = floor(time),
                        month = round(1+((time-styear)-floor(time-styear))*12,0),
                        file=files.ascii)
  #subset----
  if(annualmaps) ascii.df.sub = subset(ascii.df, type%in%do.vars & (group.num%in%do.grps | fleet.num%in%do.fleets) & year%in%do.years)
  if(!annualmaps) ascii.df.sub = subset(ascii.df, type%in%do.vars & (group.num%in%do.grps | fleet.num%in%do.fleets) & year%in%do.years & month%in%do.months)
  ascii.df.sub$label = ifelse(ascii.df.sub$type=='effort',
                              paste(ascii.df.sub$type,formatC(ascii.df.sub$fleet.num,width=2,flag="0"),ascii.df.sub$year,formatC(ascii.df.sub$month,width=2,flag="0"),sep="_"),
                              paste(ascii.df.sub$type,formatC(ascii.df.sub$group.num,width=2,flag="0"),ascii.df.sub$year,formatC(ascii.df.sub$month,width=2,flag="0"),sep="_"))
  
  #stack and return----
  #out.stack <- stack(ascii.df.sub$file)
  out.stack <- rast(ascii.df.sub$file)
  names(out.stack) <- ascii.df.sub$label
  return(out.stack)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecospace per-region prediction readers------------------------------------------
# Shared readers used by fn.ecospace_plot_fits() and fn.ecospace_objfxn(). Reads
# the per-region Ecospace_Annual_Average_Region_<n>_<Var>.csv files into a 4-D
# array [years, groups (or fleet|group), regions, runs].

#' @keywords internal
#' @noRd
.find_year_skip <- function(x, timestep){
  key_fun <- if(timestep == "annual")
               function(L) which(tolower(substr(L, 1, 4)) == "year")
             else
               function(L) which(tolower(substr(L, 1, 8)) == "timestep")
  for(nread in c(30L, 100L, 500L, 5000L)){
    L <- readLines(x, n = nread)
    hit <- key_fun(L)
    if(length(hit)) return(hit[1] - 1)
    if(length(L) < nread) return(NA_integer_)
  }
  NA_integer_
}

#' @keywords internal
#' @noRd
.detect_ecospace_regions <- function(dir.out, timestep){
  prefix <- if(timestep == "annual") "Ecospace_Annual_Average_Region_" else "Ecospace_Average_Region_"
  files <- unique(unlist(lapply(dir.out, function(d){
    list.files(d, pattern = paste0("^", prefix, "\\d+_Biomass\\.csv$"),
               recursive = TRUE)
  })))
  if(length(files) == 0) return(integer(0))
  m <- regmatches(basename(files), regexec("Region_(\\d+)_", basename(files)))
  regs <- as.integer(vapply(m, `[`, character(1), 2))
  sort(unique(regs))
}

#' @title Read region areas for area-weighted multi-region aggregation.
#' @description Ecospace per-region output is a spatially averaged \emph{density}
#'   (\code{"Average regional biomass by group (t/km^2)"}), so a timeseries assigned
#'   to several regions must be combined as an area-weighted mean, not a sum or a
#'   plain mean. Regions differ enormously in size -- in the WFS MICE grid region 12
#'   is a single cell (74 km^2) while region 10 is 734 cells (55,956 km^2) -- so an
#'   unweighted combination lets one cell carry as much of the predicted signal as
#'   several hundred.
#' @details Region 0 is the \strong{whole grid} in Ecospace output
#'   (\code{Ecospace_Annual_Average_Region_0_Biomass.csv} carries numerically
#'   identical data to the whole-grid file; only the \code{Data}/\code{Area} label
#'   rows differ). Its weight is therefore the sum of every non-land region, NOT the
#'   attributes row whose ID happens to be 0 -- in the WFS MICE table that row is
#'   "Water (unsampled)", a genuine sub-region of the grid. Areas are taken in km^2
#'   rather than cell counts so that the latitude-varying size of a 5-min cell is
#'   respected.
#' @param path Path to the region attributes CSV. Defaults to the WFS MICE
#'   \code{combined_regions_5min_attributes.csv}.
#' @param id.col,area.col Column names holding the region ID and its area.
#' @param land.ids Region IDs to drop as land / outside the model domain.
#' @return Named numeric vector of areas, names = region ID as character.
#' @seealso \code{\link{fn.ecospace_objfxn}}, which consumes this as
#'   \code{region.areas}.
#' @examples
#' \dontrun{region.areas <- fn.read_region_areas()}
#' @export
fn.read_region_areas <- function(path     = NULL,
                                 id.col   = "ID",
                                 area.col = "AREA_KM2",
                                 land.ids = -9999){
  if(is.null(path))
    path <- file.path("C:/Users/dchagaris/OneDrive - University of Florida",
                      "WFS Fisheries Ecosystem Modeling/WFS EwE/Ecospace/maps/regions",
                      "combined_regions_5min_attributes.csv")
  if(!file.exists(path)) stop("Region attributes file not found: ", path)
  a <- utils::read.csv(path, stringsAsFactors = FALSE)
  miss <- setdiff(c(id.col, area.col), names(a))
  if(length(miss))
    stop("Region attributes file is missing column(s): ", paste(miss, collapse = ", "))

  a <- a[!(a[[id.col]] %in% land.ids), , drop = FALSE]
  w <- stats::setNames(as.numeric(a[[area.col]]),
                       as.character(as.integer(a[[id.col]])))
  bad <- !is.finite(w) | w <= 0
  if(any(bad)){
    warning("Dropping region(s) with non-positive or non-finite area: ",
            paste(names(w)[bad], collapse = ", "))
    w <- w[!bad]
  }
  if(length(w) == 0) stop("No usable region areas found in: ", path)
  w["0"] <- sum(w)   # region 0 = whole grid in Ecospace output
  w
}

#' @keywords internal
#' @noRd
# Resolve area weights for a set of region IDs being combined into one series.
# Returns equal weights (current behaviour) when there is nothing to weight or
# when region.areas cannot cover the requested IDs, so callers that never set
# region.areas keep working. The warning fires once per session per cause.
.region_weights <- function(region_ids_sel, region.areas = NULL){
  n <- length(region_ids_sel)
  if(n <= 1) return(rep(1, max(n, 1)))
  if(is.null(region.areas)){
    .warn_once("region.areas.missing",
               "Timeseries spans multiple regions but region.areas is not set -- ",
               "combining regions with EQUAL weights. Per-region Ecospace output is ",
               "a density (t/km^2), so this over-weights small regions. Set ",
               "region.areas <- fn.read_region_areas() to fix.")
    return(rep(1, n))
  }
  w <- region.areas[as.character(region_ids_sel)]
  if(anyNA(w) || !all(is.finite(w)) || sum(w, na.rm = TRUE) <= 0){
    .warn_once("region.areas.incomplete",
               "region.areas has no entry for region(s) ",
               paste(region_ids_sel[is.na(w)], collapse = ", "),
               " -- falling back to equal weights for the affected series.")
    return(rep(1, n))
  }
  as.numeric(w)
}

#' @keywords internal
#' @noRd
.warn_once <- function(key, ...){
  store <- getOption("R4EwE.warned", character(0))
  if(key %in% store) return(invisible(NULL))
  options(R4EwE.warned = c(store, key))
  warning(..., call. = FALSE)
  invisible(NULL)
}

#' @keywords internal
#' @noRd
# Per-fleet discard-mortality rate by year from the ts (type 11). Rule: a fleet
# with no type-11 entry defaults to Dmort = 1 (i.e. DT = DD). NA years within an
# existing series are filled with that fleet's mean rate.
.ecospace_dmort_by_year <- function(obs.ts, fleet_poolcodes, years){
  th       <- obs.ts$ts.head
  ts_years <- suppressWarnings(as.numeric(rownames(obs.ts$ts)))
  out <- matrix(1, nrow = length(years), ncol = length(fleet_poolcodes))
  for(j in seq_along(fleet_poolcodes)){
    fc <- fleet_poolcodes[j]
    if(is.na(fc)) next
    rows <- which(th$Type == 11 & th$Poolcode == fc)
    if(length(rows) == 0) next
    series <- obs.ts$ts[, rows[1]]
    series[series < 0] <- NA
    fill <- mean(series, na.rm = TRUE)
    v <- series[match(years, ts_years)]
    v[is.na(v)] <- if(is.finite(fill)) fill else 1
    out[, j] <- v
  }
  out
}

#' @title Read Ecospace per-region prediction arrays.
#' @description Reads the per-region \code{Ecospace_Annual_Average_Region_<n>_<Var>.csv}
#'   (or \code{Ecospace_Average_Region_*} for monthly) outputs from one or more Ecospace run
#'   folders and packs them into a 4-D array
#'   \code{[years, groups (or fleet|group), regions, runs]}. When \code{aggregate_by_group = TRUE}
#'   any \code{"fleet|group"} columns are summed to the group level.
#' @param dir.out Path or vector of paths to Ecospace run folders.
#' @param varname Variable to read: typically \code{"Biomass"}, \code{"Catch"}, or \code{"Landings"}.
#' @param timestep \code{"annual"} or \code{"monthly"}.
#' @param regions Integer vector of region IDs to read.
#' @param styear Calendar start year used to label the year dimension.
#' @param aggregate_by_group Logical. If TRUE, sum \code{"fleet|group"} columns to group level.
#' @return A 4-D array \code{[years, groups, regions, runs]} or \code{NULL} if no matching files found.
#' @export
fn.read_pred_ecospace_wide <- function(dir.out, varname, timestep, regions, styear,
                                       aggregate_by_group = FALSE){
  prefix <- if(timestep == "annual") "Ecospace_Annual_Average_Region_" else "Ecospace_Average_Region_"

  collected <- list()
  for(d in dir.out){
    for(r in regions){
      pat <- paste0("^", prefix, r, "_", varname, "\\.csv$")
      files <- list.files(d, pattern = pat, recursive = TRUE, full.names = TRUE)
      for(f in files){
        run_key <- basename(dirname(f))
        nskip <- .find_year_skip(f, timestep)
        if(length(nskip) == 0 || is.na(nskip)) next
        dt <- utils::read.csv(f, skip = nskip, row.names = 1, check.names = FALSE)
        if(aggregate_by_group && any(grepl("\\|", names(dt)))){
          grp_names <- vapply(strsplit(names(dt), "\\|"), function(x) tail(x, 1), character(1))
          uniq_grps <- unique(grp_names)
          agg <- matrix(0, nrow = nrow(dt), ncol = length(uniq_grps),
                        dimnames = list(rownames(dt), uniq_grps))
          for(g in uniq_grps){
            cg <- which(grp_names == g)
            agg[, g] <- if(length(cg) == 1) dt[, cg] else rowSums(dt[, cg, drop = FALSE], na.rm = TRUE)
          }
          dt <- as.data.frame(agg, check.names = FALSE)
        }
        if(is.null(collected[[run_key]])) collected[[run_key]] <- list()
        collected[[run_key]][[as.character(r)]] <- dt
      }
    }
  }
  if(length(collected) == 0) return(NULL)

  all_groups <- unique(unlist(lapply(collected, function(rl)
    unique(unlist(lapply(rl, names))))))
  template <- collected[[1]][[which(!vapply(collected[[1]], is.null, logical(1)))[1]]]
  n_years <- nrow(template)
  runs <- names(collected)

  arr <- array(NA_real_,
               dim = c(n_years, length(all_groups), length(regions), length(runs)),
               dimnames = list(seq_len(n_years), all_groups, paste0("R", regions), runs))
  for(r_idx in seq_along(runs)){
    for(reg_idx in seq_along(regions)){
      dt <- collected[[runs[r_idx]]][[as.character(regions[reg_idx])]]
      if(is.null(dt)) next
      arr[seq_len(nrow(dt)), names(dt), reg_idx, r_idx] <- as.matrix(dt)
    }
  }
  dimnames(arr)[[1]] <- as.character(styear + seq_len(n_years) - 1L)
  arr
}

#' @title Split Ecospace discards into dead, total, and surviving.
#' @description Derives predicted dead / total / surviving discards from per-region per-fleet-group
#'   Catch and Landings CSVs. Uses \code{DD = Catch - Landings} (dead discards),
#'   \code{DT = DD / Dmort} (total discards, with per-fleet Dmort taken from the type-11
#'   DiscardMortality series in \code{obs.ts}; a fleet with no type-11 entry defaults to
#'   Dmort = 1 so DT = DD), and \code{DS = DT - DD} (surviving discards). Computed at
#'   fleet x group resolution then summed to group level.
#' @param dir.out Path or vector of paths to Ecospace run folders.
#' @param timestep \code{"annual"} or \code{"monthly"}.
#' @param regions Integer vector of region IDs.
#' @param styear Calendar start year used to label the year dimension.
#' @param obs.ts Output of \code{\link{fn.read_ecosim_timeseries}}, providing the per-fleet
#'   type-11 DiscardMortality series.
#' @param fleet.names Character vector indexed by fleet pool code, used to map the
#'   \code{"fleet|group"} column names to fleet pool codes. If \code{NULL}, every fleet
#'   defaults to Dmort = 1 (total discards collapse to dead discards).
#' @return A list with elements \code{dead}, \code{total}, \code{surv} (each
#'   \code{[years, groups, regions, runs]}) and the fleet|group versions
#'   \code{dead_fg}, \code{total_fg}, \code{surv_fg}; or \code{NULL} if Catch/Landings files
#'   are missing.
#' @export
fn.read_pred_ecospace_discards_split <- function(dir.out, timestep, regions, styear,
                                                 obs.ts, fleet.names){
  catch <- fn.read_pred_ecospace_wide(dir.out, "Catch",    timestep, regions, styear, aggregate_by_group = FALSE)
  land  <- fn.read_pred_ecospace_wide(dir.out, "Landings", timestep, regions, styear, aggregate_by_group = FALSE)
  if(is.null(catch) || is.null(land)) return(NULL)

  common <- intersect(dimnames(catch)[[2]], dimnames(land)[[2]])
  common <- common[grepl("\\|", common)]
  if(length(common) == 0) return(NULL)
  catch <- catch[, common, , , drop = FALSE]
  land  <- land [, common, , , drop = FALSE]

  dd <- catch - land
  dd[dd < 0] <- 0

  years    <- as.numeric(dimnames(dd)[[1]])
  parts    <- strsplit(common, "\\|")
  fleet_nm <- vapply(parts, function(x) trimws(x[1]),        character(1))
  group_nm <- vapply(parts, function(x) trimws(tail(x, 1)),  character(1))

  norm  <- function(x) gsub("\\s+", " ", trimws(tolower(x)))
  fl_pc <- if(!is.null(fleet.names)) match(norm(fleet_nm), norm(fleet.names)) else rep(NA_integer_, length(fleet_nm))
  unmatched <- unique(fleet_nm[is.na(fl_pc)])
  if(length(unmatched))
    warning("Ecospace discards: no fleet.names match for fleet(s) ",
            paste(unmatched, collapse = ", "), " - discard mortality defaulted to 1 (DT = DD).")

  uniq_pc   <- unique(fl_pc)
  dmort_mat <- .ecospace_dmort_by_year(obs.ts, uniq_pc, years)

  dt <- dd; ds <- dd
  nreg <- dim(dd)[3]; nrun <- dim(dd)[4]
  for(k in seq_along(common)){
    dmk <- dmort_mat[, match(fl_pc[k], uniq_pc)]
    for(ri in seq_len(nreg)) for(rn in seq_len(nrun)){
      ddv <- dd[, k, ri, rn]
      dtv <- ddv / dmk
      dtv[!is.finite(dtv)] <- 0
      dt[, k, ri, rn] <- dtv
      ds[, k, ri, rn] <- dtv - ddv
    }
  }

  agg_to_group <- function(arr4){
    ug  <- unique(group_nm)
    out <- array(0, dim = c(dim(arr4)[1], length(ug), dim(arr4)[3], dim(arr4)[4]),
                 dimnames = list(dimnames(arr4)[[1]], ug, dimnames(arr4)[[3]], dimnames(arr4)[[4]]))
    for(g in ug){
      cg <- which(group_nm == g)
      out[, g, , ] <- if(length(cg) == 1) arr4[, cg, , ]
                      else apply(arr4[, cg, , , drop = FALSE], c(1, 3, 4), sum, na.rm = TRUE)
    }
    out
  }

  attach_meta <- function(arr4){
    attr(arr4, "fg_fleet_nm") <- setNames(fleet_nm, common)
    attr(arr4, "fg_group_nm") <- setNames(group_nm, common)
    arr4
  }

  list(dead     = agg_to_group(dd),
       total    = agg_to_group(dt),
       surv     = agg_to_group(ds),
       dead_fg  = attach_meta(dd),
       total_fg = attach_meta(dt),
       surv_fg  = attach_meta(ds))
}

#' @title Parse Ecospace \code{fleet|group} column names.
#' @description Splits Ecospace \code{"fleet|group"} column-name strings into a data frame of
#'   fleet and group name parts. Use with the column dimnames of arrays returned by
#'   \code{\link{fn.read_pred_ecospace_wide}} to map columns to fleet and group pool codes.
#' @param fg_names Character vector of \code{"fleet|group"} column names.
#' @return A data frame with columns \code{fg}, \code{fleet_nm}, \code{group_nm}.
#' @export
fn.fg_meta <- function(fg_names){
  parts <- strsplit(fg_names, "\\|")
  data.frame(fg       = fg_names,
             fleet_nm = vapply(parts, function(x) trimws(x[1]),       character(1)),
             group_nm = vapply(parts, function(x) trimws(tail(x, 1)), character(1)),
             stringsAsFactors = FALSE)
}

#' @title Other-mortality loss rate from Ecospace output.
#' @description Computes the annual M0 loss rate (loss / mean biomass) per group from a pair of
#'   Ecospace annual-average output CSVs. Multistanza groups (columns named like
#'   \code{"gag 0"}, \code{"gag 1"}, ..., \code{"gag 5+"}) are detected by suffix; for each
#'   such base group an additional \code{"<base> total"} column is added containing the
#'   biomass-weighted aggregate rate: \code{sum(loss) / sum(biomass)} across stanzas. This
#'   matches the gag/red-grouper Mrt calculations in
#'   \code{red tide sims/get red tide mortality S88.R}, generalized so the function can be
#'   applied to any other-mortality source (red tide, hypoxia, harmful algal bloom, etc.).
#' @param dir.out Character vector of one or more Ecospace run output directories. Each is
#'   expected to contain both \code{file.bio} and \code{file.loss}.
#' @param file.bio Biomass CSV filename. Defaults to the annual-average file.
#' @param file.loss Other-mortality-loss CSV filename. Defaults to the annual-average file.
#' @param out.file Output CSV filename written into each \code{dir.out} when \code{write=TRUE}.
#' @param skip Number of metadata header lines to skip when reading data (default 31, matches
#'   EwE 6.7 Ecospace output). The function also preserves the blank line after the metadata
#'   so the written output mirrors the structure of the source CSV.
#' @param stanza.regex Regular expression matching the stanza-age suffix on a column name. The
#'   default \code{"\\s+[0-9]+\\+?$"} matches the EwE convention of \code{"<group> <age>"} or
#'   \code{"<group> 5+"} for terminal stanzas. The base group name is everything before the
#'   match.
#' @param write Logical; if \code{TRUE} (default), write \code{out.file} into each input dir.
#'   Header info is copied verbatim from the source biomass CSV, with the \code{Data,} line
#'   updated to reflect the new content.
#' @return If \code{length(dir.out) == 1L}, a numeric matrix
#'   \code{[years x (n_groups + n_stanza_totals)]} of rates only: the source year/time column is
#'   dropped from the matrix and its values are carried on \code{rownames()}. Otherwise a 3-D
#'   array \code{[years x groups x runs]} whose \code{dimnames[[1]]} are the years/time index and
#'   whose 3rd dimension is named by \code{basename(dir.out)}. (The on-disk CSV written when
#'   \code{write=TRUE} still re-attaches the year column to mirror the source file structure.)
#' @details Division by zero (biomass = 0) yields \code{NaN}/\code{Inf} in the returned rate,
#'   matching the source mortality scripts; downstream summaries should use \code{na.rm=TRUE}
#'   in percentile / mean calls.
#' @export
fn.M0_loss_rate <- function(dir.out,
                            file.bio  = "Ecospace_Annual_Average_Biomass.csv",
                            file.loss = "Ecospace_Annual_Average_OtherMortalityLoss.csv",
                            out.file  = "M0_loss_rate.csv",
                            skip      = 31,
                            stanza.regex = "\\s+[0-9]+\\+?$",
                            write     = TRUE){
  dir.out  <- as.character(dir.out)
  n.runs   <- length(dir.out)
  if(n.runs < 1L) stop("dir.out must contain at least one directory.")
  out.list <- vector("list", n.runs)

  for(r in seq_len(n.runs)){
    bio.file  <- file.path(dir.out[r], file.bio)
    loss.file <- file.path(dir.out[r], file.loss)
    if(!file.exists(bio.file))  stop("Missing biomass file: ",  bio.file)
    if(!file.exists(loss.file)) stop("Missing loss file: ",     loss.file)

    # Preserve the metadata block plus the blank line that separates it from the
    # column-header row, so the output CSV mirrors the source file's structure.
    header.lines <- readLines(bio.file, n = skip + 1L)

    bio  <- utils::read.csv(bio.file,  as.is = TRUE, skip = skip, check.names = FALSE)
    loss <- utils::read.csv(loss.file, as.is = TRUE, skip = skip, check.names = FALSE)
    if(!identical(dim(bio), dim(loss)))
      stop("Biomass and loss CSVs have different shapes in ", dir.out[r])
    if(!identical(names(bio), names(loss)))
      stop("Biomass and loss CSVs have different column names in ", dir.out[r])

    # First column is the time index (Year / Time step). Peel it off so the
    # returned matrix is rates only, with years carried on dimnames[[1]].
    year.col <- names(bio)[1]
    years    <- as.character(bio[[year.col]])
    grp.cols <- names(bio)[-1]

    rate <- as.matrix(loss[, grp.cols, drop = FALSE]) / as.matrix(bio[, grp.cols, drop = FALSE])
    rownames(rate) <- years

    # Multistanza aggregation: group columns by the part of the name BEFORE the
    # trailing " <age>" suffix and append a biomass-weighted total per base group.
    is.stanza <- grepl(stanza.regex, grp.cols)
    if(any(is.stanza)){
      base.names   <- sub(stanza.regex, "", grp.cols[is.stanza])
      stanza.split <- split(grp.cols[is.stanza], base.names)
      for(bn in names(stanza.split)){
        cols <- stanza.split[[bn]]
        if(length(cols) < 2L) next  # singleton -> not actually multistanza
        bio.sum  <- rowSums(bio[,  cols, drop = FALSE])
        loss.sum <- rowSums(loss[, cols, drop = FALSE])
        rate <- cbind(rate, setNames(data.frame(loss.sum / bio.sum), paste0(bn, " total")))
        rate <- as.matrix(rate); rownames(rate) <- years
      }
    }

    if(write){
      # Re-attach the Year column to the on-disk CSV so the output mirrors the
      # source EwE file structure (header + Year column + group columns); the
      # in-memory return value stays Year-free with years on rownames.
      out.path <- file.path(dir.out[r], out.file)
      data.idx <- grep("^Data,", header.lines)
      if(length(data.idx) == 1L)
        header.lines[data.idx] <- "Data,\"Other mortality loss rate (1/yr)\""
      writeLines(header.lines, out.path)
      csv.df <- data.frame(check.names = FALSE,
                           setNames(list(bio[[year.col]]), year.col))
      for(cn in colnames(rate)) csv.df[[cn]] <- rate[, cn]
      utils::write.table(csv.df, out.path, sep = ",", row.names = FALSE,
                         col.names = TRUE, append = TRUE, quote = FALSE)
    }

    out.list[[r]] <- rate
  }

  if(n.runs == 1L) return(out.list[[1]])

  # Multi-run: stack into a 3-D array. All runs must share shape AND column set;
  # this is the typical case when feeding a directory of GA-final-pop runs.
  dims <- dim(out.list[[1]])
  cols <- colnames(out.list[[1]])
  if(!all(vapply(out.list, function(x) identical(dim(x), dims),       logical(1L))))
    stop("Cannot stack: run output matrices have inconsistent shapes.")
  if(!all(vapply(out.list, function(x) identical(colnames(x), cols),  logical(1L))))
    stop("Cannot stack: run output matrices have inconsistent column sets.")
  arr <- array(NA_real_, dim = c(dims, n.runs),
               dimnames = list(rownames(out.list[[1]]), cols, basename(dir.out)))
  for(r in seq_len(n.runs)) arr[,,r] <- out.list[[r]]
  arr
}

#' @title Ensemble summary of M0 loss rates (long-form for downstream use).
#' @description Reduces an Mrt ensemble to per-year per-group summary statistics
#'   (mean + arbitrary percentiles) across an ensemble of model runs. The output
#'   is a long-form data frame with one row per (year, group) suitable for
#'   plotting (ggplot2) or feeding into stock-assessment projection workflows.
#'   Either pass a pre-computed Mrt array (typically from
#'   \code{\link{fn.M0_loss_rate}}) via \code{mrt}, or a vector of model output
#'   directories via \code{dir.out} (the function will call \code{fn.M0_loss_rate}
#'   first; \code{...} is forwarded).
#' @param dir.out Character vector of Ecospace run output directories. Ignored if
#'   \code{mrt} is supplied. Required otherwise.
#' @param mrt A 3-D array \code{[years x groups x runs]} (multi-run) or a 2-D
#'   matrix \code{[years x groups]} (single run) as returned by
#'   \code{\link{fn.M0_loss_rate}}. If supplied, \code{dir.out} is ignored.
#' @param probs Numeric vector of percentiles to compute across runs. Default
#'   \code{c(0.05, 0.5, 0.95)} (90\% ensemble interval + median). \code{0.5} is
#'   renamed to \code{"median"}; other values become \code{"q05"}, \code{"q95"}, etc.
#' @param styear Optional integer start year. If supplied, the \code{year} column
#'   becomes \code{styear, styear+1, ...}; otherwise the year labels are taken from
#'   \code{dimnames(mrt)[[1]]} (the row names set by \code{\link{fn.M0_loss_rate}}),
#'   falling back to the row index \code{1..N} if those are absent/non-numeric.
#' @param out.file Optional path. If supplied, the long-form summary is written
#'   as CSV (one row per year-group combination).
#' @param ... Forwarded to \code{\link{fn.M0_loss_rate}} when \code{mrt} is NULL
#'   and \code{dir.out} is supplied.
#' @return A data frame with columns \code{year}, \code{group}, \code{mean}, and
#'   one column per requested percentile. For a single-run input, the
#'   percentile columns equal the mean (no ensemble to summarize across).
#' @export
fn.M0_loss_rate_summary <- function(dir.out  = NULL,
                                    mrt      = NULL,
                                    probs    = c(0.05, 0.5, 0.95),
                                    styear   = NULL,
                                    out.file = NULL,
                                    ...){
  if(is.null(mrt)){
    if(is.null(dir.out)) stop("Provide either `dir.out` or a pre-computed `mrt` array.")
    mrt <- fn.M0_loss_rate(dir.out = dir.out, write = FALSE, ...)
  }

  # Promote single-run matrix to a degenerate 3-D array so the summary path is
  # uniform. The percentile columns will collapse to the mean (no ensemble).
  if(length(dim(mrt)) == 2L){
    mrt <- array(mrt,
                 dim = c(dim(mrt), 1L),
                 dimnames = c(dimnames(mrt), list("run1")))
  }
  if(length(dim(mrt)) != 3L)
    stop("`mrt` must be a 2-D matrix (single run) or 3-D array (years x groups x runs).")

  # mrt no longer carries a year column: every column is a group, and the
  # years/time index live on dimnames(mrt)[[1]] (see fn.M0_loss_rate).
  group_cols <- dimnames(mrt)[[2]]

  years <- if(!is.null(styear)){
    styear + seq_len(dim(mrt)[1]) - 1L
  } else {
    rn <- suppressWarnings(as.numeric(dimnames(mrt)[[1]]))
    if(length(rn) == dim(mrt)[1] && !anyNA(rn)) rn else seq_len(dim(mrt)[1])
  }

  sub <- mrt[, group_cols, , drop = FALSE]

  agg <- list(mean = apply(sub, c(1, 2), mean, na.rm = TRUE))
  for(p in probs){
    nm <- if(isTRUE(all.equal(p, 0.5))) "median" else sprintf("q%02d", round(p * 100))
    agg[[nm]] <- apply(sub, c(1, 2), stats::quantile, probs = p, na.rm = TRUE, names = FALSE)
  }

  long <- do.call(rbind, lapply(group_cols, function(g){
    df <- data.frame(year = years, group = g, stringsAsFactors = FALSE)
    for(nm in names(agg)) df[[nm]] <- agg[[nm]][, g]
    df
  }))

  if(!is.null(out.file)){
    utils::write.csv(long, out.file, row.names = FALSE)
  }
  long
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.akaike_weights--------------------------------------------------------------------------------
#' @title Akaike-style weights from a vector of LL values
#' @description
#' Convenience helper to convert a vector of per-run LL (negative log-likelihood)
#' values into Akaike-style weights for use with \code{\link{fn.weighted_mrt_summary}}.
#' The weight for run \code{i} is
#' \code{w_i = exp(-(LL_i - min(LL)) / scale)} normalized to sum to 1.
#'
#' @param LL Numeric vector of LL values, one per run.
#' @param scale Scaling factor in the exponent. Default 2, mimicking the standard
#'   AIC-weight convention that treats LL as \code{-2 log L}. Increase \code{scale}
#'   to flatten the weights (broader ensemble); decrease it to concentrate on the
#'   very best fits. For a wide LL spread you may want \code{scale = 2 * df_eff}
#'   where \code{df_eff} is some effective degrees of freedom.
#' @return Numeric vector of weights summing to 1 (same length as \code{LL}).
#' @examples
#' \dontrun{
#' # From a GA fitness vector (final-gen min_LL = best individual)
#' w <- fn.akaike_weights(fitness_g20, scale = 2)
#' summary(w); sum(w)   # should sum to 1
#' }
#' @export
fn.akaike_weights <- function(LL, scale = 2){
  if(!is.numeric(LL)) stop("`LL` must be numeric.")
  if(scale <= 0)      stop("`scale` must be > 0.")
  delta <- LL - min(LL, na.rm = TRUE)
  w     <- exp(-delta / scale)
  w / sum(w, na.rm = TRUE)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.weighted_mrt_summary--------------------------------------------------------------------------
#' @title Weighted Mrt ensemble summary
#' @description
#' Long-form ensemble summary of per-run Mrt using arbitrary per-run weights
#' (e.g. Akaike-style weights derived from each individual's LL). Mirrors
#' \code{\link{fn.M0_loss_rate_summary}} but replaces the plain ensemble mean
#' and quantiles with weighted versions, so individuals with better fit
#' contribute more to the band.
#'
#' Weighted mean per (year, group):
#'   \code{sum(x_i * w_i) / sum(w_i)} (NA-safe; runs with NA at that cell drop).
#' Weighted quantiles per (year, group): linear interpolation on the empirical
#'   weighted CDF (\code{cumsum(w_sorted) / sum(w_sorted)}).
#'
#' @param mrt 3-D array \code{[years x groups x runs]} of per-run Mrt values
#'   (typically built by stacking outputs of \code{\link{fn.M0_loss_rate}}).
#'   A 2-D \code{[years x groups]} matrix is also accepted (degenerate
#'   single-run case).
#' @param weights Numeric vector of length \code{dim(mrt)[3]} (one weight per
#'   run). Non-positive and NA weights are dropped before summarizing. The
#'   function normalizes the remaining weights internally; you do not need to
#'   normalize beforehand.
#' @param probs Numeric vector of percentiles in [0,1]. Default
#'   \code{c(0.05, 0.5, 0.95)}.
#' @param styear Optional integer start year. If supplied, the \code{year}
#'   column becomes \code{styear, styear+1, ...}; otherwise the year labels
#'   come from \code{dimnames(mrt)[[1]]} (see \code{\link{fn.M0_loss_rate}}).
#' @param out.file Optional path; if supplied, the long-form summary is written
#'   to CSV.
#' @return Data frame with columns \code{year, group, weighted_mean}, then one
#'   column per requested percentile (\code{q05}, \code{median}, \code{q95}, ...).
#' @examples
#' \dontrun{
#' # Stack per-run Mrt into a 3-D array
#' mrt <- abind::abind(lapply(final_pop_dirs, fn.M0_loss_rate, write = FALSE),
#'                     along = 3)
#'
#' # Akaike-style weights from GA fitness vector
#' w <- fn.akaike_weights(fitness_g20, scale = 2)
#'
#' fn.weighted_mrt_summary(mrt, weights = w,
#'                         probs = c(0.05, 0.5, 0.95),
#'                         out.file = "mrt_weighted_summary.csv")
#' }
#' @export
fn.weighted_mrt_summary <- function(mrt,
                                    weights,
                                    probs    = c(0.05, 0.5, 0.95),
                                    styear   = NULL,
                                    out.file = NULL){
  # Promote single-run matrix to a degenerate 3-D array.
  if(length(dim(mrt)) == 2L){
    mrt <- array(mrt,
                 dim = c(dim(mrt), 1L),
                 dimnames = c(dimnames(mrt), list("run1")))
  }
  if(length(dim(mrt)) != 3L)
    stop("`mrt` must be a 2-D matrix or 3-D array (years x groups x runs).")
  if(length(weights) != dim(mrt)[3])
    stop(sprintf("length(weights) (%d) must equal dim(mrt)[3] (%d, number of runs).",
                 length(weights), dim(mrt)[3]))

  ok <- !is.na(weights) & weights > 0
  if(!any(ok)) stop("No positive non-NA weights -- nothing to summarize.")
  mrt_use <- mrt[, , ok, drop = FALSE]
  w       <- weights[ok]

  group_cols <- dimnames(mrt_use)[[2]]
  years <- if(!is.null(styear)){
    styear + seq_len(dim(mrt_use)[1]) - 1L
  } else {
    rn <- suppressWarnings(as.numeric(dimnames(mrt_use)[[1]]))
    if(length(rn) == dim(mrt_use)[1] && !anyNA(rn)) rn else seq_len(dim(mrt_use)[1])
  }

  # Weighted mean per (year, group). NA-safe: skip runs with NA at that cell.
  wtd_mean <- apply(mrt_use, c(1, 2), function(x){
    keep <- !is.na(x)
    if(!any(keep)) return(NA_real_)
    sum(x[keep] * w[keep]) / sum(w[keep])
  })

  # Weighted quantiles per (year, group), all percentiles in one pass per cell.
  q_cube <- apply(mrt_use, c(1, 2), function(x){
    keep <- !is.na(x)
    if(!any(keep)) return(rep(NA_real_, length(probs)))
    xs <- x[keep]; ws <- w[keep]
    ord <- order(xs)
    xs <- xs[ord]; ws <- ws[ord]
    cw <- cumsum(ws) / sum(ws)
    stats::approx(cw, xs, xout = probs, ties = "ordered", rule = 2)$y
  })
  # Force 3-D shape (n_probs, n_years, n_groups) even when length(probs) == 1.
  if(is.matrix(q_cube)) dim(q_cube) <- c(1L, dim(q_cube))

  agg <- list(weighted_mean = wtd_mean)
  for(i in seq_along(probs)){
    p  <- probs[i]
    nm <- if(isTRUE(all.equal(p, 0.5))) "median" else sprintf("q%02d", round(p * 100))
    agg[[nm]] <- q_cube[i, , , drop = TRUE]
  }

  long <- do.call(rbind, lapply(group_cols, function(g){
    df <- data.frame(year = years, group = g, stringsAsFactors = FALSE)
    for(nm in names(agg)) df[[nm]] <- agg[[nm]][, g]
    df
  }))

  if(!is.null(out.file)){
    utils::write.csv(long, out.file, row.names = FALSE)
  }
  long
}


#' Slim-copy a generation of Ecospace GA run folders
#'
#' Copies a directory of Ecospace GA run folders (e.g. `gen30/`) into a sibling
#' output directory (default `<dir.in>_small`) that keeps only annual timeseries
#' plus spatial biomass/catch maps for a chosen set of groups, and optionally
#' effort maps for a chosen set of fleets. Folder and file names are preserved
#' verbatim so downstream R4EwE readers work unchanged. Files are copied byte-
#' for-byte (no re-parse).
#'
#' What gets copied per run folder:
#' \itemize{
#'   \item `Ecospace_Annual_Average_*.csv` — all kept.
#'   \item `Ecospace_Average_*.csv` (monthly averages) — dropped.
#'   \item `EcospaceMap(Biomass|Catch)-<group> <stanza>.csv` — kept for groups in `groups`.
#'   \item `EcospaceMapEffort-<fleet>.csv` — kept for fleets in `fleets` (exact
#'     match on the trailing token, or `fleets = "all"` for every effort map;
#'     NULL drops them all).
#'   \item `cmd.txt` — kept if `keep.cmd`.
#'   \item Anything else (e.g. `M0_loss_rate.csv`) — kept if `keep.misc`.
#' }
#'
#' @param dir.in    Character. Directory of run folders (e.g. `.../gen30`).
#' @param dir.out   Character. Output directory. Default `paste0(dir.in, "_small")`.
#' @param groups    Character vector of group name stems (matched against the
#'   `EcospaceMap...-<stem> <stanza>.csv` filename). Default gag + red grouper.
#' @param fleets    Character vector of exact fleet names for effort maps,
#'   or the string `"all"`, or `NULL` (default) to drop every effort map.
#' @param map.kinds Character vector of spatial map kinds to consider for the
#'   group filter. Default `c("Biomass", "Catch")`. Effort is controlled by
#'   `fleets`, not this argument.
#' @param keep.cmd  Logical. Copy `cmd.txt` per run. Default TRUE.
#' @param keep.misc Logical. Copy any leftover files (e.g. `M0_loss_rate.csv`)
#'   that are neither annual TS nor spatial maps nor cmd.txt. Default TRUE.
#' @param overwrite Logical. If `dir.out` already exists non-empty, delete it
#'   before copying. Default FALSE (function errors out).
#' @param ncores    Integer. Number of worker processes to run the per-subdir
#'   copy in parallel. Default 1 (sequential). On Windows the outer loop is
#'   the bottleneck; setting `ncores = 8`-`16` typically gives near-linear
#'   speedup until disk saturates.
#' @param long.paths Logical. On Windows, prefix `dir.in`/`dir.out` internally
#'   with the `\\?\` extended-length path escape so files under long parent
#'   paths (e.g. deep OneDrive folders) copy past the 259-char MAX_PATH limit.
#'   The user-facing `dir.out` in the return value is left unprefixed. Ignored
#'   on non-Windows. Default TRUE.
#' @param verbose   Logical. Print progress + summary. Default TRUE.
#'
#' @return Invisibly, a list with `dir.out`, `n_dirs`, `n_files_in`,
#'   `n_files_out`, `bytes_in`, `bytes_out`, and any `warnings` collected
#'   (e.g. fleets requested that no run dir contained).
#' @export
fn.gen_slim_copy <- function(dir.in,
                             dir.out    = paste0(dir.in, "_small"),
                             groups     = c("gag", "red grouper"),
                             fleets     = NULL,
                             map.kinds  = c("Biomass", "Catch"),
                             keep.cmd   = TRUE,
                             keep.misc  = TRUE,
                             overwrite  = FALSE,
                             ncores     = 1L,
                             long.paths = TRUE,
                             verbose    = TRUE){

  if(!dir.exists(dir.in))
    stop("dir.in does not exist: ", dir.in)
  if(!length(groups))
    stop("groups must be non-empty")

  subs <- list.dirs(dir.in, recursive = FALSE, full.names = FALSE)
  if(!length(subs))
    stop("dir.in contains no subdirectories: ", dir.in)

  if(dir.exists(dir.out)){
    existing <- list.files(dir.out, all.files = TRUE, no.. = TRUE)
    if(length(existing)){
      if(!isTRUE(overwrite))
        stop("dir.out exists and is non-empty; pass overwrite = TRUE: ", dir.out)
      unlink(dir.out, recursive = TRUE, force = TRUE)
    }
  }
  dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)

  # Extended-length path escape hatch (Windows-only). Turns
  #   C:/deep/long/path -> \\?\C:\deep\long\path  (or \\?\UNC\srv\share\...)
  # so file.copy / list.files / file.info bypass the 259-char MAX_PATH limit.
  use_long <- isTRUE(long.paths) && .Platform$OS.type == "windows"
  .ext_path <- function(p){
    if(!use_long) return(p)
    if(startsWith(p, "\\\\?\\")) return(p)
    ap <- suppressWarnings(normalizePath(p, winslash = "\\", mustWork = FALSE))
    if(startsWith(ap, "\\\\")) paste0("\\\\?\\UNC\\", substring(ap, 3))
    else                        paste0("\\\\?\\", ap)
  }
  # Only dir.out actually exceeds MAX_PATH in practice (a deep OneDrive dest
  # under a short filename-prefix source dir). Prefix ONLY the dest side and
  # keep list.files / file.info on the source with plain paths — R's Windows
  # list.files uses the CRT _wfindfirst, which doesn't accept the \\?\ escape.
  dir_out_ext <- .ext_path(dir.out)
  # \\?\ paths require backslash separators; file.path on Windows uses "/".
  # Under long.paths, join dest parts with paste(sep="\\") to keep them pure.
  join_dst <- if(use_long) function(a, b) paste(a, b, sep = "\\")
              else                          function(a, b) file.path(a, b)

  re_annual <- "^Ecospace_Annual_Average_.*\\.csv$"
  re_spatial_grp <- paste0(
    "^EcospaceMap(", paste(map.kinds, collapse = "|"), ")-(",
    paste(vapply(groups, function(g) gsub("([][{}()+*^$.|?\\\\])", "\\\\\\1", g),
                 character(1)), collapse = "|"),
    ") [0-9]+\\+?\\.csv$"
  )
  effort_all <- identical(fleets, "all")
  effort_set <- if(effort_all || is.null(fleets)) character(0)
                else paste0("EcospaceMapEffort-", fleets, ".csv")

  # Vectorized filter: one grepl per rule across all filenames.
  filter_mask <- function(files){
    m_ann  <- grepl(re_annual,       files)
    m_grp  <- grepl(re_spatial_grp,  files)
    m_eff  <- grepl("^EcospaceMapEffort-.*\\.csv$", files)
    m_cmd  <- isTRUE(keep.cmd)  & files == "cmd.txt"
    m_bcx  <- grepl("^EcospaceMap(Biomass|Catch)-", files)
    m_avg  <- grepl("^Ecospace_Average_", files)
    if(effort_all)             m_eff_keep <- m_eff
    else if(length(effort_set)) m_eff_keep <- m_eff & (files %in% effort_set)
    else                        m_eff_keep <- rep(FALSE, length(files))
    m_misc <- isTRUE(keep.misc) & !(m_ann | m_grp | m_eff | m_cmd | m_bcx | m_avg)
    m_ann | m_grp | m_eff_keep | m_cmd | m_misc
  }

  if(!use_long){
    worst_src  <- max(nchar(list.files(file.path(dir.in, subs[1L]))), na.rm = TRUE)
    worst_path <- nchar(file.path(dir.out, subs[which.max(nchar(subs))], "x")) - 1L + worst_src
    if(worst_path > 259L)
      warning("Worst-case output path is ", worst_path,
              " chars (> Windows MAX_PATH 259). Set long.paths=TRUE or use a shorter dir.out.")
  }

  # Per-subdir worker: returns list of stats (no side channels via <<-,
  # so this is safe under parLapply). Bytes taken from source, not dest,
  # since file.copy is byte-perfect.
  copy_one <- function(sub){
    sub_src <- file.path(dir.in,      sub)   # plain: list.files/file.info
    sub_dst <- join_dst(dir_out_ext,  sub)   # extended-length-safe dest
    dir.create(sub_dst, recursive = TRUE, showWarnings = FALSE)

    files <- list.files(sub_src, full.names = FALSE)
    if(!length(files))
      return(list(n_in = 0L, n_out = 0L, b_in = 0, b_out = 0, missing = character(0)))

    inf   <- file.info(file.path(sub_src, files))
    keep_m <- filter_mask(files)
    kept   <- files[keep_m]

    missing_here <- if(!is.null(fleets) && !effort_all)
                      setdiff(effort_set, files) else character(0)

    if(!length(kept))
      return(list(n_in = length(files), n_out = 0L,
                  b_in = sum(inf$size, na.rm = TRUE), b_out = 0,
                  missing = missing_here))

    ok <- file.copy(file.path(sub_src, kept), join_dst(sub_dst, kept),
                    overwrite = TRUE, copy.date = TRUE)

    list(n_in  = length(files),
         n_out = sum(ok),
         b_in  = sum(inf$size, na.rm = TRUE),
         b_out = sum(inf$size[keep_m][ok], na.rm = TRUE),
         missing = missing_here,
         failed  = sum(!ok))
  }

  # Dispatch: sequential with progress bar, or parallel via parLapply.
  ncores <- max(1L, as.integer(ncores))
  if(ncores == 1L){
    pb <- if(isTRUE(verbose)) utils::txtProgressBar(min = 0, max = length(subs), style = 3) else NULL
    res <- vector("list", length(subs))
    for(i in seq_along(subs)){
      res[[i]] <- copy_one(subs[i])
      if(!is.null(pb)) utils::setTxtProgressBar(pb, i)
    }
    if(!is.null(pb)){ close(pb); cat("\n") }
  } else {
    if(isTRUE(verbose))
      cat(sprintf("[slim-copy] parallel: %d workers x %d dirs\n", ncores, length(subs)))
    cl <- parallel::makeCluster(min(ncores, length(subs)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl,
      c("dir.in", "dir_out_ext", "re_annual", "re_spatial_grp",
        "effort_all", "effort_set", "keep.cmd", "keep.misc",
        "fleets", "filter_mask", "join_dst"),
      envir = environment())
    res <- parallel::parLapply(cl, subs, copy_one)
  }

  # Fold results.
  n_files_in  <- sum(vapply(res, `[[`, integer(1), "n_in"))
  n_files_out <- sum(vapply(res, `[[`, integer(1), "n_out"))
  bytes_in    <- sum(vapply(res, `[[`, numeric(1), "b_in"))
  bytes_out   <- sum(vapply(res, `[[`, numeric(1), "b_out"))
  n_failed    <- sum(vapply(res, function(x) if(is.null(x$failed)) 0L else as.integer(x$failed), integer(1)))
  warn_missing_fleets <- unique(unlist(lapply(res, `[[`, "missing")))

  if(n_failed > 0L)
    warning("Failed to copy ", n_failed, " file(s) across all subdirs.")

  if(length(warn_missing_fleets))
    warning("Requested fleet effort map(s) not present in any run dir: ",
            paste(gsub("^EcospaceMapEffort-|\\.csv$", "", warn_missing_fleets),
                  collapse = ", "))

  if(isTRUE(verbose)){
    mb <- function(x) sprintf("%.1f MB", x / 1024^2)
    cat(sprintf("[slim-copy] %d dirs  %d -> %d files  %s -> %s (%.0f%% saved)\n",
                length(subs), n_files_in, n_files_out,
                mb(bytes_in), mb(bytes_out),
                100 * (1 - bytes_out / max(bytes_in, 1))))
    cat("[slim-copy] out: ", dir.out, "\n", sep = "")
  }

  invisible(list(dir.out = dir.out,
                 n_dirs = length(subs),
                 n_files_in = n_files_in,
                 n_files_out = n_files_out,
                 bytes_in = bytes_in,
                 bytes_out = bytes_out,
                 warnings = warn_missing_fleets))
}

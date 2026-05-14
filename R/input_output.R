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

  # coerce non-Title columns to numeric
  for(j in setdiff(names(hdr), "Title")){
    hdr[[j]] <- suppressWarnings(as.numeric(hdr[[j]]))
  }
  if(!"Poolcode2" %in% names(hdr)) hdr$Poolcode2 <- NA_real_
  if(!"Region"    %in% names(hdr)) hdr$Region    <- NA_integer_

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

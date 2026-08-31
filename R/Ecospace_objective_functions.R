#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#OBJECTIVE FUNCTIONS---------------------------------------------------------------------------------  
#Loop to run objfxn for sensitivity and fill in runlist ////////////////////////
#print(paste("Sensitivity run obj. fxn.", start, "to", stop))

#...............................................................................
##Objective fxn 1-----
#fits to timeseries data only, using the ecosim timeseries csv file
#' @title Ecospace objective function, timeseries only
#' @description Calculates the composite likelihood.
#' @param dir.pred Full directory path to model output folder containing .csv files.
#' @param obs.ts A list object containing reference timeseries information, as that created by 
#' fn.read_ecosim_timeseries().
#' @param run.idx If dir.pred has multiple output folders, use this to select which to calculate.
#' @param fit.abs.catch If TRUE (default), do not rescale predicions.  If FALSE, predictions were 
#' divided by q, where q=mean(pred)/mean(obs)
#' @return A vector containing the total and functional group specific likelihood.  
#' @export
fn.objfxn1 <- function(dir.pred, obs.ts=obs.ts, run.idx=1, fit.abs.catch=TRUE){  
  #get annual timeseries predictions
  #dir.pred = "C:/NWACS MICE/GA output 2026-01-06/GA_Run_20260106_114454/run_00abeb6e03a8897940ee4aa45b593526"
  predB = fn.ecospace_predB_ts2array(dir.out = dir.pred, timestep='annual', n.reg=0)[,,run.idx]
  predC = fn.ecospace_predC_ts2array(dir.out = dir.pred, timestep='annual', n.reg=0)[,,run.idx] #unable to read region catch so far
  dim(predB); dim(predC)
  ###biomass timeseries annual----
  obs.ts.head = obs.ts$obsB.head
  obs.ts.biomass = obs.ts$obsB
  names(obs.ts.head) = c('title','weight','poolcode','type')
  lk.ts.biomass <- cbind(obs.ts.head, dattype='biomass timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=5
    #grp = group.names[j]
    grpnum = obs.ts.head$poolcode[j]
    grp = group.names[grpnum]
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    lk.dat = data.frame(pred=predB[,grpnum], obs=obs.ts.biomass[,j])
    #lk.dat$pred <- ifelse(is.na(lk.dat$pred),1e-6,lk.dat$pred)
    lk.dat$obs = ifelse(lk.dat$obs==0,NA,lk.dat$obs)
    rownames(lk.dat) <- rownames(predB)
    if(obs.ts.head$type[j]==0){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$pred = lk.dat$pred/q
    lk.dat$ll = log(lk.dat$pred/lk.dat$obs)^2/(2*obs.sd^2)
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.biomass$obs.cv[j] = obs.cv
    lk.ts.biomass$obs.sd[j] = obs.sd
    lk.ts.biomass$loglik[j] = lk.sum
  }
  
  ###catch timeseries annual----
  obs.ts.head = obs.ts$obsC.head
  obs.ts.catch = obs.ts$obsC
  names(obs.ts.head) = c('title','weight','poolcode','type')
  lk.ts.catch <- cbind(obs.ts.head, dattype='catch timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=2
    #grp = group.names[j]
    grpnum = obs.ts.head$poolcode[j]
    grp = gsub("/+",".",gsub("-",".",gsub("_", ".", group.names[grpnum])))
    idx_group = which(grepl(grp,colnames(predC)))
    #idx_group = match(grp,colnames(predC))
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    if(length(idx_group)==1) lk.dat = data.frame(pred=predC[,idx_group], obs=obs.ts.catch[,j])
    if(length(idx_group)>1) lk.dat = data.frame(pred=rowSums(predC[,idx_group]), obs=obs.ts.catch[,j])
    #lk.dat$pred <- ifelse(is.na(lk.dat$pred),1e-6,lk.dat$pred)
    lk.dat$obs[lk.dat$obs <=0] <- NA
    rownames(lk.dat) <- rownames(predC)
    
    if(obs.ts.head$type[j]==61 | !fit.abs.catch){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$pred = lk.dat$pred/q
    lk.dat$ll = log(lk.dat$pred/lk.dat$obs)^2/(2*obs.sd^2)
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.catch$obs.cv[j] = obs.cv
    lk.ts.catch$obs.sd[j] = obs.sd
    lk.ts.catch$loglik[j] = lk.sum
  } 
  
  lk.ts <<- rbind(lk.ts.biomass,lk.ts.catch)
  
  ## Make output vector................................
  ##
  lk.vec <- rep(0,length(group.names))
  lk.agg <- aggregate(loglik~poolcode, data=lk.ts, sum,na.rm=T)
  lk.vec[lk.agg$poolcode] <- lk.agg$loglik
  
  outvec = round(c(sum(lk.vec), lk.vec),2)
  names(outvec) = c('LL.total',paste0('LL.',group.names))
  
  #write.csv(outvec,paste0(runlist$dir.out,"/objvals.csv"))
  return(outvec)
}


#fits to timeseries data and some spatial data.  UNDER DEVELOPMENT!!!
#' @title Ecospace objective function, timeseries plut spatial reference data
#' @description  Calculates the composite likelihood, including spatial reference datasets.
#' @param dir.pred Full directory path to model output folder containing .csv files.
#' @param obs.ts A list object containing reference timeseries information, as that created by 
#' fn.read_ecosim_timeseries().
#' @param obs.maps A raster stack???
#' @param obs.maps.meta A dataframe??? Likely nrow=nlyrs(obs.maps) If TRUE (default), do not rescale
#' predicions.  If FALSE, predictions were divided by q, where q=mean(pred)/mean(obs)
#' @param autoweight.LL Logical value.  If TRUE (default), the timeseries and maps likelihoods are 
#' each scaled to their own means LL.w=LL*1/mean(LL)  This puts them on equal footing, but needs more thought.
#' @return A vector containing the total and functional group specific likelihood.  
#' @keywords internal
#' @export
fn.objfxn2 <- function(dir.pred, obs.ts=obs.ts, obs.maps=obs.maps, obs.maps.meta=obs.maps.meta, autoweight.LL=T){
  #get annual timeseries predictions
  predB = fn.ecospace_predB_ts2array(dir.out = dir.pred, timestep='annual', n.reg=0)[,,1]
  predC = fn.ecospace_predC_ts2array(dir.out = dir.pred, timestep='annual', n.reg=0)[,,1] #unable to read region catch so far
  
  ###biomass timeseries annual----
  obs.ts.head = obs.ts$obsB.head
  obs.ts.biomass = obs.ts$obsB
  names(obs.ts.head) = c('title','weight','poolcode','type')
  lk.ts.biomass <- cbind(obs.ts.head, dattype='biomass timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=1
    #grp = group.names[j]
    grpnum = obs.ts.head$poolcode[j]
    grp = group.names[grpnum]
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    lk.dat = data.frame(pred=predB[,grpnum], obs=obs.ts.biomass[,j])
    rownames(lk.dat) <- rownames(predB)
    if(obs.ts.head$type[j]==0){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$pred = lk.dat$pred/q
    lk.dat$ll = log(lk.dat$pred/lk.dat$obs)^2/(2*obs.sd^2)
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.biomass$obs.cv[j] = obs.cv
    lk.ts.biomass$obs.sd[j] = obs.sd
    lk.ts.biomass$loglik[j] = lk.sum
  }
  
  ###catch timeseries annual----
  obs.ts.head = obs.ts$obsC.head
  obs.ts.catch = obs.ts$obsC
  names(obs.ts.head) = c('title','weight','poolcode','type')
  lk.ts.catch <- cbind(obs.ts.head, dattype='catch timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=2
    #grp = group.names[j]
    grpnum = obs.ts.head$poolcode[j]
    grp = gsub("/+",".",gsub("-",".",gsub("_", ".", group.names[grpnum])))
    idx_group = which(grepl(grp,colnames(predC)))
    #idx_group = match(grp,colnames(predC))
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    lk.dat = data.frame(pred=predC[,idx_group], obs=obs.ts.catch[,j])
    lk.dat$obs[lk.dat$obs <=0] <- NA
    rownames(lk.dat) <- rownames(predC)
    
    if(obs.ts.head$type[j]==61){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$pred = lk.dat$pred/q
    lk.dat$ll = log(lk.dat$pred/lk.dat$obs)^2/(2*obs.sd^2)
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.catch$obs.cv[j] = obs.cv
    lk.ts.catch$obs.sd[j] = obs.sd
    lk.ts.catch$loglik[j] = lk.sum
  } 
  
  lk.ts = rbind(lk.ts.biomass,lk.ts.catch)
  
  #fit to raster grids----
  #probability distribution maps (maps should sum to 1)
  lk.maps = cbind(obs.maps.meta, dattype='prob dist map',obs.cv=NA, obs.sd=NA, loglik=NA)
  
  for(i in 1:nlayers(obs.maps)){
    #i=1
    obs.i = rast(obs.maps[[i]])
    spp.i = strsplit(names(obs.i),"_")[[1]][1]
    grpnum.i = as.numeric(strsplit(names(obs.i),"_")[[1]][3])
    grpname.i = gsub("_"," ",df.names$group.names[grpnum.i])
    source.i = obs.maps.meta$source[i]
    files.pred.asc = list.files(path=file.path(dir.pred,'asc'), pattern=paste0("EcospaceMapBiomass-",grpname.i), full.names=T)
    
    if(i==1){
      tsteps <- as.numeric(substr(basename(files.pred.asc),nchar(basename(files.pred.asc))-8,nchar(basename(files.pred.asc))-4)  )
      tsteps.y <- rep(styear:(styear-1+length(tsteps)/12),each=12)
      tsteps.m <- rep(1:12,length(tsteps)/12)
    }
    
    if(source.i=='NEFSC sdm rasters'){
      get.files.i = which(tsteps.y>=obs.maps.meta$yr.start[i] & tsteps.y<=obs.maps.meta$yr.end[i] &
                            tsteps.m>=obs.maps.meta$mo.start[i] & tsteps.m<=obs.maps.meta$mo.end[i] )
      stack.i = terra::rast(files.pred.asc[get.files.i])
      mean.i = terra::mean(stack.i,na.rm=T) #terra::app(stack.i, fun=mean, na.rm=T)
      prop.i <- mean.i / global(mean.i, sum, na.rm=TRUE)[1,1]
    }
    obs_vals <- values(obs.i)
    pred_vals <- values(prop.i)
    valid_idx <- !is.na(obs_vals) & !is.na(pred_vals)
    obs_probs <- obs_vals[valid_idx] / sum(obs_vals[valid_idx])
    pred_probs <- pred_vals[valid_idx] / sum(pred_vals[valid_idx])
    
    # Compute log-likelihood (cross-entropy)
    lk.maps$loglik[i] <- -sum(obs_probs * log(pred_probs))
  }
  
  
  ## combine likelihoods
  if(autoweight.LL){
    lk.ts$loglik.w = lk.ts$loglik*(1/mean(lk.ts$loglik))
    lk.maps$loglik.w = lk.maps$loglik*(1/mean(lk.maps$loglik))
  } else{
    lk.ts$loglik.w = lk.ts$loglik
    lk.maps$loglik.w = lk.maps$loglik
  }
  lk.comb = data.frame(obs.name=c(lk.ts$title, lk.maps$layername),
                       poolcode = as.numeric(c(lk.ts$poolcode, lk.maps$grp.num)),
                       dattype = c(lk.ts$dattype, lk.maps$dattype),
                       loglik = c(lk.ts$loglik, lk.maps$loglik),
                       loglik.w = c(lk.ts$loglik.w, lk.maps$loglik.w))
  
  lk.agg <- aggregate(loglik.w~poolcode, data=lk.comb, sum,na.rm=T)

  ## Make output vector................................
  lk.vec <- rep(0,length(group.names))
  lk.vec[lk.agg$poolcode] <- lk.agg$loglik.w

  outvec = round(c(sum(lk.vec), lk.vec),2)
  names(outvec) = c('LL.total',paste0('LL.',group.names))

  #write.csv(outvec,paste0(runlist$dir.out,"/objvals.csv"))
  return(outvec)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#...............................................................................
##Objective fxn ecospace-----
#' @title Ecospace objective function: timeseries with fleet-specific landings and discards.
#' @description Composite lognormal cost function that fits Ecospace predictions to observed
#'   timeseries in \code{obs.ts} (as produced by \code{\link{fn.read_ecosim_timeseries}}).
#'   Supports biomass (type 0, 1), aggregated catch (6, 61), fleet-specific landings (12),
#'   fleet-specific discards (13, treated as TOTAL discards), and group-level total discards
#'   (19, 20). Region-aware: each obs series is matched to the prediction for its declared
#'   region from the \code{Region} header row in \code{obs.ts}. Reuses the shared prediction
#'   readers \code{\link{fn.read_pred_ecospace_wide}} and
#'   \code{\link{fn.read_pred_ecospace_discards_split}} so cost and plots share matching logic.
#' @param dir.pred Full path (or vector of paths) to Ecospace run folder(s) containing per-region
#'   \code{Ecospace_Annual_Average_Region_<n>_<Var>.csv} outputs.
#' @param obs.ts A list as produced by \code{\link{fn.read_ecosim_timeseries}} (must include
#'   \code{ts.head} with \code{Type}, \code{Poolcode}, \code{Poolcode2}, \code{Region},
#'   \code{Weight}, \code{Absolute}).
#' @param group.names Character vector indexed by group pool code.
#' @param fleet.names Character vector indexed by fleet pool code. Required to fit type-12/13
#'   fleet-specific series; without it those series are skipped.
#' @param run.idx Integer; which prediction run to evaluate when multiple runs are present.
#'   Default 1.
#' @param regions Integer vector of region IDs to consider. \code{NULL} (default) auto-detects
#'   all \code{Region_<n>_Biomass.csv} files under \code{dir.pred}.
#' @param fit.abs.catch Logical. If TRUE (default) absolute series are fit on their literal
#'   scale (\code{q = 1}); only relative series are q-rescaled. If FALSE, every series is
#'   q-rescaled by \code{q = mean(pred[overlap]) / mean(obs[overlap])}.
#' @param types Optional integer vector of obs type codes to fit. Default fits every supported
#'   type present in \code{obs.ts$ts.head$Type}.
#' @param eps Small floor applied to predictions before \code{log()} to avoid \code{log(0)}.
#'   Default \code{.Machine$double.eps^0.5}.
#' @param spatial.obs Optional list of observed spatial maps (per group/pool) added as a
#'   spatial cross-entropy term via \code{\link{fn.spatial_LL}}. Defaults to \code{spatial.obs}
#'   in the global environment, or \code{NULL} to skip the spatial term.
#' @param spatial.weight Scalar weight (prior precision) on the spatial likelihood term.
#'   Defaults to \code{spatial.weight} in the global environment, or 1.
#' @param model_styear First model year, used to align predicted maps with the year window of
#'   each observed map. Defaults to \code{model_styear} in the global environment, or 1985.
#' @param region.areas Named numeric vector of region areas (km^2), names = region ID as
#'   character; see \code{\link{fn.read_region_areas}}. Used to combine a series that spans
#'   several regions as an area-weighted mean of the per-region densities. Defaults to
#'   \code{region.areas} in the global environment, or \code{NULL} -- in which case regions
#'   are combined with equal weights (with a warning), which over-weights small regions.
#' @param return One of \code{"both"} (default), \code{"vector"}, or \code{"detail"}.
#'   \code{"both"} returns \code{list(LL = vector, detail = data.frame)};
#'   \code{"vector"} returns just the per-group cost vector (drop-in for legacy
#'   \code{fn.objfxn1} callers / tight GA worker loops);
#'   \code{"detail"} returns just the per-series data frame. The detail data frame is
#'   \code{obs.ts$ts.head} (Title, Weight, Poolcode, Poolcode2, Type, Region, Type.name,
#'   Absolute, Reference, PoolRole) for the subset of series that were considered, with
#'   the computed columns \code{q}, \code{n}, \code{obs.cv}, \code{obs.sd}, \code{loglik}
#'   appended.
#' @return The shape selected by \code{return}. Series with \code{Weight == 0} are recorded
#'   in the detail table but contribute zero to the cost.
#' @export
fn.ecospace_objfxn <- function(dir.pred,
                               obs.ts,
                               group.names,
                               fleet.names   = NULL,
                               run.idx       = 1,
                               regions       = NULL,
                               fit.abs.catch = TRUE,
                               types         = NULL,
                               eps           = NULL,
                               spatial.obs   = if(exists("spatial.obs",    envir = .GlobalEnv))
                                                 get("spatial.obs",       envir = .GlobalEnv) else NULL,
                               spatial.weight= if(exists("spatial.weight", envir = .GlobalEnv))
                                                 get("spatial.weight",    envir = .GlobalEnv) else 1,
                               model_styear  = if(exists("model_styear",   envir = .GlobalEnv))
                                                 get("model_styear",      envir = .GlobalEnv) else 1985,
                               region.areas  = if(exists("region.areas",   envir = .GlobalEnv))
                                                 get("region.areas",      envir = .GlobalEnv) else NULL,
                               return        = c("both", "vector", "detail")){
  return_mode <- match.arg(return)
  timestep    <- "annual"
  if(is.null(eps)) eps <- .Machine$double.eps^0.5

  # ---- resolve regions ----
  if(is.null(regions)) regions <- .detect_ecospace_regions(dir.pred, timestep)
  if(length(regions) == 0)
    stop("No Ecospace per-region output files found under: ", paste(dir.pred, collapse = ", "))

  # ---- infer styear from obs.ts ----
  obs_yrs <- suppressWarnings(as.numeric(rownames(obs.ts$ts)))
  styear  <- min(obs_yrs[is.finite(obs_yrs)])
  if(!is.finite(styear))
    stop("Could not infer styear from rownames(obs.ts$ts).")

  # ---- pick obs types to fit ----
  supported    <- c(0, 1, 6, 61, 12, 13, 19, 20)
  present      <- intersect(supported, unique(obs.ts$ts.head$Type))
  types_to_fit <- if(is.null(types)) present else intersect(types, present)
  if(length(types_to_fit) == 0)
    stop("No supported obs types present in obs.ts (looking for any of: ",
         paste(supported, collapse = ", "), ").")

  # ---- read predictions once, only what's needed ----
  pred <- list()
  if(any(types_to_fit %in% c(0, 1)))
    pred$biomass <- fn.read_pred_ecospace_wide(dir.pred, "Biomass", timestep, regions, styear,
                                               aggregate_by_group = FALSE)
  if(any(types_to_fit %in% c(6, 61)))
    pred$catch   <- fn.read_pred_ecospace_wide(dir.pred, "Catch",   timestep, regions, styear,
                                               aggregate_by_group = TRUE)
  if(any(types_to_fit %in% c(12)))
    pred$landings_fg <- fn.read_pred_ecospace_wide(dir.pred, "Landings", timestep, regions, styear,
                                                   aggregate_by_group = FALSE)
  if(any(types_to_fit %in% c(13, 19, 20))){
    ds <- fn.read_pred_ecospace_discards_split(dir.pred, timestep, regions, styear,
                                               obs.ts, fleet.names)
    if(!is.null(ds)){
      if(any(types_to_fit %in% c(13)))     pred$disc_total_fg <- ds$total_fg
      if(any(types_to_fit %in% c(19, 20))) pred$disc_total    <- ds$total
    }
  }

  # ---- name normalizers (match plotting.R conventions) ----
  norm_name  <- function(x) gsub("\\s+", "_", trimws(gsub('"', "", as.character(x))))
  norm_fleet <- function(x) gsub("\\s+", " ", trimws(tolower(x)))
  norm_gn    <- norm_name(group.names)
  norm_fn    <- if(!is.null(fleet.names)) norm_fleet(fleet.names) else character(0)

  any_arr       <- pred[[which(!vapply(pred, is.null, logical(1)))[1]]]
  region_labels <- if(!is.null(any_arr)) dimnames(any_arr)[[3]] else paste0("R", regions)
  region_ids    <- as.integer(sub("R", "", region_labels))

  # ---- (fleet_pc, group_pc) -> column index lookup for fleet|group arrays ----
  fg_lookup <- function(arr){
    if(is.null(arr)) return(NULL)
    cn <- dimnames(arr)[[2]]
    fg <- cn[grepl("\\|", cn)]
    if(length(fg) == 0) return(NULL)
    meta <- fn.fg_meta(fg)
    meta$fleet_pc <- if(length(norm_fn)) match(norm_fleet(meta$fleet_nm), norm_fn) else NA_integer_
    meta$group_pc <- match(norm_name(meta$group_nm), norm_gn)
    meta$col_idx  <- match(meta$fg, cn)
    meta
  }
  meta_landings_fg <- fg_lookup(pred$landings_fg)
  meta_disc_fg     <- fg_lookup(pred$disc_total_fg)

  group_col <- function(arr, gpc){
    if(is.null(arr) || is.na(gpc) || gpc < 1 || gpc > length(group.names)) return(NA_integer_)
    match(norm_name(group.names[gpc]), norm_name(dimnames(arr)[[2]]))
  }

  pick_pred <- function(tc){
    switch(as.character(tc),
           "0"  = list(arr = pred$biomass,       kind = "group"),
           "1"  = list(arr = pred$biomass,       kind = "group"),
           "6"  = list(arr = pred$catch,         kind = "group"),
           "61" = list(arr = pred$catch,         kind = "group"),
           "12" = list(arr = pred$landings_fg,   kind = "fg", meta = meta_landings_fg),
           "13" = list(arr = pred$disc_total_fg, kind = "fg", meta = meta_disc_fg),
           "19" = list(arr = pred$disc_total,    kind = "group"),
           "20" = list(arr = pred$disc_total,    kind = "group"),
           NULL)
  }

  # ---- walk obs rows by type and accumulate cost ----
  th      <- obs.ts$ts.head
  ts_data <- obs.ts$ts


  # parallel accumulators (one slot per considered obs series)
  acc_idx    <- integer(0)
  acc_grp_pc <- integer(0)
  acc_q      <- numeric(0)
  acc_n      <- integer(0)
  acc_cv     <- numeric(0)
  acc_sd     <- numeric(0)
  acc_ll     <- numeric(0)

  for(tc in types_to_fit){
    sel <- pick_pred(tc)
    if(is.null(sel) || is.null(sel$arr)) next
    arr_yrs    <- as.numeric(dimnames(sel$arr)[[1]])
    common_idx <- match(obs_yrs, arr_yrs)
    obs_rows   <- which(th$Type == tc)

    for(j in obs_rows){
      wt       <- th$Weight[j]
      abs_flag <- isTRUE(th$Absolute[j])
      pc1      <- th$Poolcode[j]
      pc2      <- th$Poolcode2[j]

      # Multi-pool / multi-region aggregation is enabled for group-kind biomass
      # (types 0, 1). Fleet|group types (12, 13) and single-group catch types
      # (6, 61, 19, 20) keep scalar semantics.
      is_group_agg <- sel$kind == "group" && tc %in% c(0, 1)

      if(is_group_agg){
        # list-column Poolcodes / Regions parsed by fn.read_ecosim_timeseries;
        # back-compat: empty list -> scalar Poolcode / Region.
        poolcodes_j <- if("Poolcodes" %in% names(th)) th$Poolcodes[[j]] else integer(0)
        if(length(poolcodes_j) == 0)
          poolcodes_j <- if(is.na(pc1)) integer(0) else as.integer(pc1)
        regions_j <- if("Regions" %in% names(th)) th$Regions[[j]] else integer(0)
        if(length(regions_j) == 0){
          reg_obs   <- th$Region[j]
          regions_j <- if(is.na(reg_obs)) integer(0) else as.integer(reg_obs)
        }

        # region resolution:
        #   empty      -> whole-grid R0 if present, else first region (legacy)
        #   0          -> whole-grid R0
        #   otherwise  -> match to region_ids
        r0_idx <- if(0L %in% region_ids) which(region_ids == 0L) else 1L
        if(length(regions_j) == 0 || all(is.na(regions_j))){
          reg_idxs <- r0_idx
        } else if(any(regions_j == 0L)){
          reg_idxs <- r0_idx
        } else {
          reg_idxs <- match(regions_j, region_ids)
          if(anyNA(reg_idxs)) next
        }

        # pool code resolution -> group column indices
        cols <- vapply(poolcodes_j,
                       function(p) group_col(sel$arr, p),
                       integer(1))
        if(length(cols) == 0 || anyNA(cols)) next

        # Sum predicted biomass across pool codes (the ages of one stock), then
        # take the AREA-WEIGHTED MEAN across regions. Per-region Ecospace output
        # is a density (t/km^2), so summing or plain-averaging across regions
        # lets a 1-cell region count as much as a 734-cell one. Weights are
        # region areas in km^2; a single region reduces to the old behaviour.
        reg_w  <- .region_weights(region_ids[reg_idxs], region.areas)
        pred_v <- rep(0, length(common_idx))
        for(k in seq_along(reg_idxs)){
          s <- rep(0, length(common_idx))
          for(ci in cols)
            s <- s + sel$arr[, ci, reg_idxs[k], run.idx][common_idx]
          pred_v <- pred_v + reg_w[k] * s
        }
        pred_v <- pred_v / sum(reg_w)

        gpc <- poolcodes_j[1]      # first pool code labels the LL row
        col <- cols[1]
        reg_idx <- reg_idxs[1]
      } else {
        # scalar path for fleet|group and single-group catch types
        reg_obs <- th$Region[j]
        if(is.na(reg_obs)){
          reg_idx <- 1L
        } else {
          reg_idx <- match(reg_obs, region_ids)
          if(is.na(reg_idx)) next
        }
        if(sel$kind == "group"){
          gpc <- pc1
          col <- group_col(sel$arr, gpc)
        } else { # fg: Poolcode = fleet, Poolcode2 = group
          if(is.null(sel$meta)) next
          m <- which(sel$meta$fleet_pc == pc1 & sel$meta$group_pc == pc2)
          col <- if(length(m) == 0) NA_integer_ else sel$meta$col_idx[m[1]]
          gpc <- pc2
        }
        if(is.na(col)) next
        pred_v <- sel$arr[, col, reg_idx, run.idx][common_idx]
      }

      obs_v <- ts_data[, j]
      obs_v[obs_v <= 0] <- NA

      # q rescale decision
      q <- 1
      if(!abs_flag || !fit.abs.catch){
        ok <- !is.na(obs_v) & !is.na(pred_v) & pred_v > 0
        if(any(ok)){
          mp <- mean(pred_v[ok]); mo <- mean(obs_v[ok])
          if(is.finite(mp) && mp > 0 && is.finite(mo) && mo > 0) q <- mp / mo
        }
      }

      n_pts  <- sum(!is.na(obs_v) & !is.na(pred_v))

      if(is.na(wt) || wt == 0){
        # weight-0 series: recorded but does not contribute to cost
        acc_idx    <- c(acc_idx, j)
        acc_grp_pc <- c(acc_grp_pc, gpc)
        acc_q      <- c(acc_q, q)
        acc_n      <- c(acc_n, n_pts)
        acc_cv     <- c(acc_cv, NA_real_)
        acc_sd     <- c(acc_sd, NA_real_)
        acc_ll     <- c(acc_ll, NA_real_)
        next
      }

      obs.cv <- 1 / wt
      obs.sd <- sqrt(log(1 + obs.cv^2))

      pred_scaled <- pmax(pred_v / q, eps)
      ll          <- (log(pred_scaled / obs_v))^2 / (2 * obs.sd^2)
      ll_sum      <- sum(ll, na.rm = TRUE)

      acc_idx    <- c(acc_idx, j)
      acc_grp_pc <- c(acc_grp_pc, gpc)
      acc_q      <- c(acc_q, q)
      acc_n      <- c(acc_n, n_pts)
      acc_cv     <- c(acc_cv, obs.cv)
      acc_sd     <- c(acc_sd, obs.sd)
      acc_ll     <- c(acc_ll, ll_sum)
    }
  }

  # ---- build detail dataframe: obs.ts$ts.head columns + computed columns ----
  if(length(acc_idx) > 0){
    lk.detail <- obs.ts$ts.head[acc_idx, , drop = FALSE]
    lk.detail$q      <- acc_q
    lk.detail$n      <- acc_n
    lk.detail$obs.cv <- acc_cv
    lk.detail$obs.sd <- acc_sd
    lk.detail$loglik <- acc_ll
    rownames(lk.detail) <- NULL
  } else {
    lk.detail <- cbind(obs.ts$ts.head[integer(0), , drop = FALSE],
                       q      = numeric(0),
                       n      = integer(0),
                       obs.cv = numeric(0),
                       obs.sd = numeric(0),
                       loglik = numeric(0))
  }

  # ---- aggregate to per-group cost vector ----
  lk.vec <- rep(0, length(group.names))
  if(length(acc_idx) > 0){
    agg_df <- data.frame(grp_pc = acc_grp_pc, loglik = acc_ll)
    agg    <- aggregate(loglik ~ grp_pc, data = agg_df, sum, na.rm = TRUE)
    ok     <- !is.na(agg$grp_pc) & agg$grp_pc >= 1 & agg$grp_pc <= length(group.names)
    lk.vec[agg$grp_pc[ok]] <- agg$loglik[ok]
  }

  outvec <- round(c(sum(lk.vec, na.rm = TRUE), lk.vec), 2)
  names(outvec) <- c("LL.total", paste0("LL.", group.names))

  # ---- optional spatial LL component (phase 2 task 7c) ----
  spatial_vec <- numeric(0)
  spatial_total <- 0
  if(!is.null(spatial.obs) && is.data.frame(spatial.obs) && nrow(spatial.obs) > 0){
    sp <- tryCatch(
      fn.spatial_LL(dir.pred      = dir.pred[1],
                    spatial.obs   = spatial.obs,
                    spatial.weight = spatial.weight,
                    group.names   = group.names,
                    model_styear  = model_styear),
      error = function(e){
        warning("fn.spatial_LL failed: ", conditionMessage(e))
        NULL
      })
    if(!is.null(sp)){
      spatial_total <- sp$total
      spatial_vec   <- round(sp$per_species, 4)
      names(spatial_vec) <- paste0("LL.spatial.", names(sp$per_species))
    }
  }

  # rebuild outvec with spatial contribution injected
  ll_total  <- round(outvec[["LL.total"]] + spatial_total, 2)
  outvec[["LL.total"]] <- ll_total
  outvec <- c(outvec, spatial_vec)

  switch(return_mode,
         both   = list(LL = outvec, detail = lk.detail),
         vector = outvec,
         detail = lk.detail)
}
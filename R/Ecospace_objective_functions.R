#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#OBJECTIVE FUNCTIONS---------------------------------------------------------------------------------  
#Loop to run objfxn for sensitivity and fill in runlist ////////////////////////
#print(paste("Sensitivity run obj. fxn.", start, "to", stop))

#...............................................................................
##Objective fxn 1-----
#fits to timeseries data only, using the ecosim timeseries csv file
#' @title Ecospace objective function, timeseries only (deprecated)
#' @description Calculates the composite likelihood.
#'
#'   \strong{Deprecated -- use \code{\link{fn.ecospace_objfxn}}.} This function reads the
#'   legacy \code{obs.ts$obsB.head} / \code{obsC.head} sub-objects, is limited to biomass and
#'   aggregated catch on region 0, and assigns \code{lk.ts} into the global environment with
#'   \code{<<-}. It cannot fit regions, fishing mortality (type 104), fleet-specific landings
#'   (12), discards (13, 19, 20), or the spatial cross-entropy term, and it produces no
#'   per-component likelihood report. Retained only so fitness values from pre-2026 GA runs
#'   remain reproducible.
#' @param dir.pred Full directory path to model output folder containing .csv files.
#' @param obs.ts A list object containing reference timeseries information, as that created by
#' fn.read_ecosim_timeseries().
#' @param run.idx If dir.pred has multiple output folders, use this to select which to calculate.
#' @param fit.abs.catch If TRUE (default), do not rescale predicions.  If FALSE, predictions were
#' divided by q, where q=mean(pred)/mean(obs)
#' @return A vector containing the total and functional group specific likelihood.
#' @seealso \code{\link{fn.ecospace_objfxn}}
#' @export
fn.objfxn1 <- function(dir.pred, obs.ts=obs.ts, run.idx=1, fit.abs.catch=TRUE){
  .Deprecated("fn.ecospace_objfxn",
              msg = paste("fn.objfxn1() is deprecated: biomass and aggregated catch on",
                          "region 0 only, and it assigns lk.ts into the global environment.",
                          "Use fn.ecospace_objfxn() (GA obj.fxn = 3)."))
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
#' @param fmort.b.floor Guard for type 104 (fishing mortality), which is derived as
#'   predicted catch / predicted biomass. Within each group and region, years whose
#'   predicted biomass falls below \code{fmort.b.floor} times that group's maximum
#'   predicted biomass are set to \code{NA} rather than scored. Without this, a group
#'   the parameter search has driven toward extinction yields a catch/biomass ratio
#'   built from two numerically meaningless numbers, and because the cost is a squared
#'   log residual those years would dominate the objective and actively pull the search
#'   toward collapsed states. A fraction rather than an absolute floor because group
#'   biomass densities in a typical model span several orders of magnitude. Default
#'   \code{1e-3}; set \code{0} to disable.
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
#' @param ll.file Optional switch to write a per-component likelihood report to disk.
#'   \code{NULL} (default) writes nothing. \code{TRUE} writes \code{LL_components.csv}
#'   into \code{dir.pred[1]}; a character path writes to that path. The report holds one
#'   row per REFERENCE timeseries -- every \code{obs.ts$ts.head} row with
#'   \code{Reference == TRUE}, scored or not -- and one row per \code{spatial.obs} row,
#'   each with its likelihood contribution in \code{loglik}. Spatial values are already
#'   multiplied by \code{spatial.weight}; timeseries values are not. Unscored reference
#'   series carry \code{loglik = NA} with \code{scored = FALSE}; zero-weight series are
#'   \code{scored = TRUE} but also \code{loglik = NA}, since they contribute nothing.
#'   Every row records the scoring settings used -- \code{fit.abs.catch},
#'   \code{region.weights} (\code{"area"} if \code{region.areas} was supplied, else
#'   \code{"equal"}) and \code{spatial.weight} -- and timeseries rows add
#'   \code{abs.scale}, TRUE when the series was fit on its literal scale rather than
#'   q-rescaled. Compare re-scores of historical runs only when these match.
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
                               fmort.b.floor = 1e-3,
                               spatial.obs   = if(exists("spatial.obs",    envir = .GlobalEnv))
                                                 get("spatial.obs",       envir = .GlobalEnv) else NULL,
                               spatial.weight= if(exists("spatial.weight", envir = .GlobalEnv))
                                                 get("spatial.weight",    envir = .GlobalEnv) else 1,
                               model_styear  = if(exists("model_styear",   envir = .GlobalEnv))
                                                 get("model_styear",      envir = .GlobalEnv) else 1985,
                               region.areas  = if(exists("region.areas",   envir = .GlobalEnv))
                                                 get("region.areas",      envir = .GlobalEnv) else NULL,
                               return        = c("both", "vector", "detail"),
                               ll.file       = NULL){
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
  supported    <- c(0, 1, 6, 61, 12, 13, 19, 20, 104)
  present      <- intersect(supported, unique(obs.ts$ts.head$Type))
  types_to_fit <- if(is.null(types)) present else intersect(types, present)
  if(length(types_to_fit) == 0)
    stop("No supported obs types present in obs.ts (looking for any of: ",
         paste(supported, collapse = ", "), ").")

  # ---- read predictions once, only what's needed ----
  # Type 104 (FishingMortalityRef) is not an Ecospace output in its own right; it
  # is derived below as catch / biomass, so it needs both reads.
  pred <- list()
  if(any(types_to_fit %in% c(0, 1, 104)))
    pred$biomass <- fn.read_pred_ecospace_wide(dir.pred, "Biomass", timestep, regions, styear,
                                               aggregate_by_group = FALSE)
  if(any(types_to_fit %in% c(6, 61, 104)))
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

  # ---- derive F = catch / biomass for type 104 ----
  # Ecospace writes no fishing-mortality output, so F is reconstructed from the
  # two arrays already read. Columns are matched on normalized group name because
  # Catch arrives via aggregate_by_group = TRUE (group names rebuilt from
  # "fleet|group" headers) while Biomass arrives with group names verbatim, so the
  # two column sets are not guaranteed to share an order -- or a membership.
  if(104 %in% types_to_fit){
    if(is.null(pred$catch) || is.null(pred$biomass)){
      warning("Type 104 requested but Catch and/or Biomass output is missing; skipping F.")
      types_to_fit <- setdiff(types_to_fit, 104)
    } else {
      pred$fmort <- .derive_fmort_arr(pred$catch, pred$biomass,
                                      b.floor.frac = fmort.b.floor,
                                      norm = norm_name)
      if(is.null(pred$fmort)){
        warning("No group columns shared between Catch and Biomass output; skipping F.")
        types_to_fit <- setdiff(types_to_fit, 104)
      }
    }
  }

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
           # 104 deliberately stays out of the is_group_agg set below: F is a rate,
           # so it is neither additive across pool codes nor area-averageable across
           # regions. It is scored on the single region named by Region (or the first
           # detected region), exactly as types 6/61/19/20 are.
           "104"= list(arr = pred$fmort,         kind = "group"),
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
    lk.detail$ts.row <- acc_idx
    rownames(lk.detail) <- NULL
  } else {
    lk.detail <- cbind(obs.ts$ts.head[integer(0), , drop = FALSE],
                       q      = numeric(0),
                       n      = integer(0),
                       obs.cv = numeric(0),
                       obs.sd = numeric(0),
                       loglik = numeric(0),
                       ts.row = integer(0))
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
  spatial_per <- NULL
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
      spatial_per   <- sp$per_species
      names(spatial_vec) <- paste0("LL.spatial.", names(sp$per_species))
    }
  }

  # rebuild outvec with spatial contribution injected
  ll_total  <- round(outvec[["LL.total"]] + spatial_total, 2)
  outvec[["LL.total"]] <- ll_total
  outvec <- c(outvec, spatial_vec)


  # ---- optional per-component LL report ------------------------------------
  # One row per REFERENCE timeseries and one per spatial obs, each with its
  # likelihood contribution, so a run's total LL can be attributed back to the
  # individual data sources that produced it.
  if(!is.null(ll.file) && !identical(ll.file, FALSE)){
    ll.path <- if(isTRUE(ll.file)) file.path(dir.pred[1], "LL_components.csv")
               else as.character(ll.file)[1]
    rep_df <- .ll_component_report(obs.ts         = obs.ts,
                                   lk.detail      = lk.detail,
                                   spatial.obs    = spatial.obs,
                                   spatial_per    = spatial_per,
                                   spatial.weight = spatial.weight,
                                   fit.abs.catch  = fit.abs.catch,
                                   region.areas   = region.areas)
    ok <- tryCatch({
      d <- dirname(ll.path)
      if(!dir.exists(d)) dir.create(d, recursive = TRUE)
      utils::write.csv(rep_df, ll.path, row.names = FALSE)
      TRUE
    }, error = function(e){
      warning("Could not write LL component report to '", ll.path, "': ",
              conditionMessage(e))
      FALSE
    })
    if(ok) attr(lk.detail, "ll.file") <- ll.path
  }

  switch(return_mode,
         both   = list(LL = outvec, detail = lk.detail),
         vector = outvec,
         detail = lk.detail)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title Per-component likelihood report.
#' @description Builds the data frame written by \code{fn.ecospace_objfxn(ll.file=)}:
#'   one row for every REFERENCE timeseries in \code{obs.ts$ts.head} (whether or not it
#'   was scored) and one row for every \code{spatial.obs} row, each carrying its
#'   likelihood contribution. Timeseries rows keep the full ts.head columns plus the
#'   computed \code{q}, \code{n}, \code{obs.cv}, \code{obs.sd} and \code{abs.scale};
#'   spatial rows carry \code{pc}, \code{species}, \code{styear}, \code{endyear}.
#'   Columns absent from one block are NA-filled in the other.
#'
#'   Every row also records the run-level scoring settings -- \code{fit.abs.catch},
#'   \code{region.weights} and \code{spatial.weight} -- so a report is self-describing.
#'   Re-scoring a historical run under different settings changes the likelihood
#'   silently; these columns are what let you tell a real difference from a setting
#'   mismatch.
#' @param obs.ts Observed timeseries object (needs \code{$ts.head}).
#' @param lk.detail Per-series detail frame built by \code{fn.ecospace_objfxn}; must
#'   carry \code{ts.row} (row index into \code{ts.head}) and \code{loglik}.
#' @param spatial.obs Spatial observation table, or NULL.
#' @title Derive a predicted fishing-mortality array from catch and biomass.
#' @description Builds F = catch / biomass for \code{fn.ecospace_objfxn} obs type 104.
#'   Ecospace has no fishing-mortality writer, so F has to be reconstructed. Group
#'   columns are matched on normalized name, since the catch array is built by
#'   aggregating \code{"fleet|group"} headers while the biomass array carries group
#'   names directly; only groups present in both are returned. Years whose biomass is
#'   below \code{b.floor.frac} of that group's per-region maximum are returned as
#'   \code{NA}, so a collapsed group contributes nothing instead of contributing noise.
#' @param catch_arr,biomass_arr Arrays from \code{\link{fn.read_pred_ecospace_wide}},
#'   dimensioned year x group x region x run.
#' @param b.floor.frac Relative biomass floor; see \code{fn.ecospace_objfxn}.
#' @param norm Name normalizer, passed in so it matches the caller's convention.
#' @return An array shaped like \code{catch_arr} but restricted to the shared group
#'   columns, or \code{NULL} if no groups are shared.
#' @keywords internal
#' @noRd
.derive_fmort_arr <- function(catch_arr, biomass_arr, b.floor.frac = 1e-3, norm = identity){
  cg <- norm(dimnames(catch_arr)[[2]])
  bg <- norm(dimnames(biomass_arr)[[2]])
  shared <- intersect(cg, bg)
  if(length(shared) == 0) return(NULL)
  ci <- match(shared, cg)
  bi <- match(shared, bg)

  # Region and run dimensions must line up; both arrays come from the same reader
  # call pattern, so this is a guard against a future divergence rather than a
  # situation expected in normal use.
  if(!identical(dim(catch_arr)[c(1, 3, 4)], dim(biomass_arr)[c(1, 3, 4)]))
    stop(".derive_fmort_arr(): Catch and Biomass arrays disagree on the year, region ",
         "or run dimension.")

  cc <- catch_arr[, ci, , , drop = FALSE]
  bb <- biomass_arr[, bi, , , drop = FALSE]

  # Mask first, then divide, so the floor is applied on biomass and not on the ratio.
  if(is.finite(b.floor.frac) && b.floor.frac > 0){
    for(g in seq_along(shared)) for(r in seq_len(dim(bb)[3])) for(u in seq_len(dim(bb)[4])){
      v <- bb[, g, r, u]
      mx <- suppressWarnings(max(v, na.rm = TRUE))
      if(!is.finite(mx) || mx <= 0){ bb[, g, r, u] <- NA_real_; next }
      v[v < b.floor.frac * mx] <- NA_real_
      bb[, g, r, u] <- v
    }
  }
  bb[!is.na(bb) & bb <= 0] <- NA_real_

  out <- cc / bb
  # Carry the catch array's own group labels, not the normalized keys, so this array
  # is interchangeable with the other pred arrays when inspected or indexed by name.
  dimnames(out) <- list(dimnames(catch_arr)[[1]], dimnames(catch_arr)[[2]][ci],
                        dimnames(catch_arr)[[3]], dimnames(catch_arr)[[4]])
  out
}

#' @param spatial_per Per-row spatial contributions from \code{fn.spatial_LL}
#'   (\code{$per_species}, already multiplied by \code{spatial.weight}), or NULL.
#' @param spatial.weight Scalar spatial weight used for the run.
#' @param fit.abs.catch The \code{fit.abs.catch} value used for the run.
#' @param region.areas The \code{region.areas} object used for the run; only whether it
#'   was supplied is recorded (\code{"area"}) or not (\code{"equal"}).
#' @return A data frame ordered timeseries-block then spatial-block.
#' @keywords internal
#' @noRd
.ll_component_report <- function(obs.ts, lk.detail, spatial.obs, spatial_per,
                                 spatial.weight = NA_real_,
                                 fit.abs.catch  = NA,
                                 region.areas   = NULL){

  # list-columns (Poolcodes, Regions) cannot go to CSV -- collapse to strings
  flatten <- function(df){
    for(nm in names(df))
      if(is.list(df[[nm]]))
        df[[nm]] <- vapply(df[[nm]],
                           function(x) paste(unlist(x), collapse = ","),
                           character(1))
    df
  }

  th <- obs.ts$ts.head
  ref <- if("Reference" %in% names(th)) which(!is.na(th$Reference) & th$Reference)
         else seq_len(nrow(th))
  ts_df <- flatten(th[ref, , drop = FALSE])

  # exact join on the ts.head row index recorded by the caller
  m <- if(!is.null(lk.detail) && nrow(lk.detail) > 0 && "ts.row" %in% names(lk.detail))
         match(ref, lk.detail$ts.row) else rep(NA_integer_, length(ref))
  pick <- function(col) if(all(is.na(m))) rep(NA_real_, length(ref)) else lk.detail[[col]][m]

  ts_df$q      <- pick("q")
  ts_df$n      <- pick("n")
  ts_df$obs.cv <- pick("obs.cv")
  ts_df$obs.sd <- pick("obs.sd")
  ts_df$loglik <- pick("loglik")
  ts_df$scored <- !is.na(m)
  # Mirrors fn.ecospace_objfxn: a series is q-rescaled when
  # `!Absolute || !fit.abs.catch`, so it is fit on its literal scale only when both
  # hold. Recorded per row because the ts.head Absolute flag alone does not say
  # which scale was actually used.
  absv <- if("Absolute" %in% names(ts_df)) !is.na(ts_df$Absolute) & ts_df$Absolute
          else rep(FALSE, nrow(ts_df))
  ts_df$abs.scale <- isTRUE(fit.abs.catch) & absv
  ts_df$abs.scale[!ts_df$scored] <- NA    # never fit, so no scale applies
  # Fishing mortality gets its own component label. It is derived as catch/biomass,
  # so in log space its residual is the catch residual minus the biomass residual --
  # it reweights information already in the objective rather than adding any. Keeping
  # it separable in the report is what lets you see that trade being made.
  ts_df$component <- "timeseries"
  if("Type" %in% names(ts_df))
    ts_df$component[!is.na(ts_df$Type) & ts_df$Type == 104] <- "fishing_mortality"
  rownames(ts_df) <- NULL

  out <- ts_df
  if(!is.null(spatial.obs) && is.data.frame(spatial.obs) && nrow(spatial.obs) > 0){
    keep  <- intersect(c("pc", "species", "styear", "endyear"), names(spatial.obs))
    sp_df <- flatten(spatial.obs[, keep, drop = FALSE])
    sp_df$loglik <- if(!is.null(spatial_per) && length(spatial_per) == nrow(sp_df))
                      as.numeric(spatial_per) else NA_real_
    sp_df$scored <- !is.na(sp_df$loglik)
    sp_df$component <- "spatial"
    # give the spatial block a Title so both blocks are identifiable in one column
    if("species" %in% names(sp_df)) sp_df$Title <- sp_df$species
    rownames(sp_df) <- NULL

    cols <- union(names(ts_df), names(sp_df))
    fill <- function(df){
      for(nm in setdiff(cols, names(df))) df[[nm]] <- NA
      df[, cols, drop = FALSE]
    }
    out <- rbind(fill(ts_df), fill(sp_df))
  }

  # run-level scoring settings, on every row
  out$fit.abs.catch  <- fit.abs.catch
  out$region.weights <- if(is.null(region.areas)) "equal" else "area"
  out$spatial.weight <- spatial.weight

  # lead with the columns a reader wants first
  front <- intersect(c("component", "Title", "loglik", "scored"), names(out))
  out[, c(front, setdiff(names(out), front)), drop = FALSE]
}

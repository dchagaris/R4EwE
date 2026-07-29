
#prepare parameter vector----
#' @title Construct parameter vectors and bounds for model estimation
#' @description Builds the parameter vector and related metadata (labels, groups, bounds, CVs) for 
#' vulnerabilities, environmental responses, dispersal, and mediation parameters. Several objects 
#' are assigned in the parent/global environment.
#' @param do.vuls Logical; include vulnerability parameters.
#' @param vul_pars Data frame of vulnerability parameters.
#' @param vul.min,vul.max Scalar override bounds for the vulnerability search range.
#'   Only used when \code{vul.bounds.mode = "override"}; ignored under the default
#'   \code{"formula"} mode.
#' @param vul.cv Coefficient of variation for vulnerabilities; controls the
#'   log-scale spread of the formula-derived bounds (effective log-sd =
#'   \code{exp(2 * vul.cv)}).
#' @param vul.bounds.mode Character; either \code{"formula"} (default — bounds
#'   are log-symmetric around each \code{base.val} using \code{vul.cv}) or
#'   \code{"override"} (use \code{[vul.min, vul.max]} as the bounds for every
#'   vul parameter, ignoring \code{base.val} and \code{vul.cv}).
#' @param do.env Logical; include environmental parameters.
#' @param env_pars Data frame of environmental parameters.
#' @param env.min Minimum environmental value.
#' @param env.max Maximum environmental value.
#' @param env.cv Coefficient of variation for environmental parameters.
#' @param do.disp Logical; include dispersal parameters.
#' @param disp_pars Data frame of dispersal parameters.
#' @param disp.min Minimum dispersal (unused placeholder).
#' @param disp.max Maximum dispersal (unused placeholder).
#' @param disp.cv Coefficient of variation for dispersal parameters.
#' @param do.med Logical; include mediation parameters.
#' @param med_pars Data frame of mediation shape parameters.
#' @param med.xbase.cv CV applied to mediation x-base parameters.
#' @param do.fltdyn Logical; include fleet dynamics parameters.
#' @param fltdyn_pars Data frame of fleet parameters with at least \code{base.val}.
#' @param fltdyn.min,fltdyn.max Safety floor/ceiling for fleet parameters. Under
#'   the default \code{"adaptive"} bounds mode they are NOT hard limits — the
#'   per-parameter bounds widen to \code{base.val ± 4 * fltdyn.cv * base.val} so
#'   the full normal bell fits inside, while \code{fltdyn.min}/\code{fltdyn.max}
#'   are used as the lower/upper bound only when the bell is tighter than them.
#'   Under \code{"scalar"} mode every fleet parameter gets these as hard scalar
#'   bounds (legacy behavior; will clip the bell if \code{base.val} exceeds them).
#' @param fltdyn.cv CV for fleet parameter draws; also drives the half-width of
#'   the adaptive bounds (\code{4 * fltdyn.cv * base.val}).
#' @param fltdyn.bounds.mode Either \code{"adaptive"} (default) or \code{"scalar"}.
#' @param do.redtide Logical; include red tide environmental response parameters.
#' @param redtide_pars Data frame of red tide response parameters (shape 11 logistic4params).
#' @param redtide.min,redtide.max Scalar multiplier bounds applied to both inflection.adj
#'   and slope.adj for red tide responses. Used as the fallback default when the
#'   per-parameter overrides below are \code{NULL}.
#' @param redtide.inf.min,redtide.inf.max Optional per-parameter bounds for the red tide
#'   \code{inflection.adj} multiplier (\code{NULL} → fall back to redtide.min/max).
#' @param redtide.slope.min,redtide.slope.max Optional per-parameter bounds for the red tide
#'   \code{slope.adj} multiplier (\code{NULL} → fall back to redtide.min/max).
#' @param redtide.cv CV used for red tide parameter draws when sampling with a normal
#'   distribution.
#' @return Invisibly returns \code{NULL}. Several objects are set globally.
#' @export
fn.makeparvec <- function(
  do.vuls = TRUE,
  vul_pars = predprey_pars,
  vul.min = 1.01,
  vul.max = 1e6,
  vul.cv = 0.4,
  vul.bounds.mode = "formula",
  
  do.env = FALSE, 
  env_pars = NULL,
  env.min=0.25,
  env.max=2,
  env.cv=0.2,
  
  do.disp = TRUE,
  disp_pars = disp_pars,
  disp.cv=0.2,
  disp.min=3,
  disp.max=3000,
  
  do.med = FALSE,
  med_pars = NULL,
  med.xbase.cv = .1,
  
  do.fltdyn = FALSE,
  fltdyn_pars=NULL,
  fltdyn.min=0,
  fltdyn.max=5,
  fltdyn.cv=0.2,
  fltdyn.bounds.mode = "adaptive",

  do.redtide = FALSE,
  redtide_pars = NULL,
  redtide.min = 0.5,
  redtide.max = 1.5,
  redtide.inf.min   = NULL,
  redtide.inf.max   = NULL,
  redtide.slope.min = NULL,
  redtide.slope.max = NULL,
  redtide.cv = 0.2){
  
  # do.vuls=fit.vuls; vul_pars=predprey_pars; vul.min=1.01; vul.max=1e6; vul.cv=0.4; 
  # do.env=fit.env; env_pars=env_pars; env.min=0.1; env.max=2; env.cv=0.3;
  # do.disp=fit.disp; disp_pars=disp_pars; disp.cv=0.2;
  # do.med = fit.med; med_pars=NULL; med.xbase.cv = .1;
  # do.fltdyn = fit.fltdyn; fltdyn_pars=fltdyn_pars; fltdyn.min=0.1; fltdyn.max=5; fltdyn.cv=0.3
  
  #validate setup----
  if(do.vuls){
    pdcol = unique(vul_pars$pred[is.na(vul_pars$prey)])
    pdpy = c(0,unique(vul_pars$pred[!is.na(vul_pars$prey)]))
    if(length(which(pdcol %in% pdpy))>0){
      stop("Trying to estimate ki and kij for at least on predator. Check vul_pars input.")
    }
  }
  
  # Normalize column conventions: accept either canonical (shape.idx + parN) or
  # sensitivity-results (shape + Param.N) schemas for env_pars and redtide_pars.
  if(do.env)     env_pars     <- .normalize_env_pars(env_pars)
  if(do.redtide) redtide_pars <- .normalize_env_pars(redtide_pars)

  #build parameter vector----
  n_vuls <- n_env <- n_disp <- n_med <- n_fleet <- n_redtide <- 0
  if(do.vuls)    n_vuls    <- nrow(vul_pars)
  if(do.env)     n_env     <- nrow(env_pars)
  if(do.disp)    n_disp    <- nrow(disp_pars)
  if(do.med)     n_med     <- nrow(med_pars)
  if(do.fltdyn)  n_fleet   <- nrow(fltdyn_pars)
  if(do.redtide) n_redtide <- nrow(redtide_pars)

  n_pars <- n_vuls + n_env + n_redtide + n_disp + n_med + n_fleet
  if (n_pars == 0) stop("n_pars==0: No predator-prey pairs or environmental responses parameters to estimate.  Check inputs.")

  #no longer taking log of parameters
  vuln_vec <- env_vec <- redtide_vec <- disp_vec <- med_vec <- fleet_vec <- numeric()
  par.idx <- par.labels <- par.groups <- character()
  #log.par.idx <- logical()
  
  if(do.vuls) {
    #log_vuln_vec <- log(vul_pars$base.val)
    vuln_vec <- vul_pars$base.val
    par.idx <- c(par.idx,rep('vul',n_vuls))
    par.labels <- c(par.labels,paste('vul',vul_pars$pred,vul_pars$prey,sep="_"))
    par.groups <- c(par.groups,paste('vul',vul_pars$pred,vul_pars$prey,sep="_"))
    #log.par.idx <- c(log.par.idx,rep(TRUE,length(log_vuln_vec)))
  }
  
  if(do.env) {
    env_idx_vec <- .env_resp_idx(env_pars)   # accepts fxn_num or resp.idx
    for(i in 1:n_env){
      shape_type.i <- env_pars$shape.idx[i]
      resp_id.i    <- env_idx_vec[i]
      # All supported shapes use 2 multiplicative GA pars centered on 1.
      adj_names.i <- switch(as.character(shape_type.i),
        "6"  = c('mean.adj',       'width.adj'),    # normal: Mean (par4); SDLeft (par1) & SDRight (par3) jointly
        "9"  = c('mid.adj',        'width.adj'),    # trapezoid: mid=(LT+RT)/2; pref=RT-LT + shoulders
        "10" = c('xmid.adj',       'slope.adj'),    # sigmoid: XMid (par3); Slope (par5)
        "11" = c('inflection.adj', 'slope.adj'),    # logistic4params: Inflection (par3); Slope (par4)
        stop("Unsupported env response shape.idx: ", shape_type.i,
             " (supported: 6 normal, 9 trapezoid, 10 sigmoid, 11 logistic4params).")
      )
      par_vec.i    <- rep(1, length(adj_names.i))
      par.labels.i <- paste(paste0('env',  resp_id.i),
                            paste0('type', shape_type.i),
                            adj_names.i, sep = "_")
      par.groups.i <- rep(paste0('env', resp_id.i), length(par_vec.i))

      env_vec    <- c(env_vec, par_vec.i)
      par.labels <- c(par.labels, par.labels.i)
      par.groups <- c(par.groups, par.groups.i)
    } #end loop
    par.idx <- c(par.idx, rep('env', length(env_vec)))
  } #end do.env

  if(do.redtide){
    redtide_idx_vec <- .env_resp_idx(redtide_pars)
    for(i in 1:n_redtide){
      shape_type.i <- redtide_pars$shape.idx[i]
      resp_id.i    <- redtide_idx_vec[i]
      # All red tide responses are shape 11 (logistic4params): GA pars inflection.adj, slope.adj.
      if(shape_type.i != 11)
        stop("Red tide response (fxn_num=", resp_id.i, ") has shape.idx=", shape_type.i,
             "; expected 11 (logistic4params).")
      adj_names.i  <- c('inflection.adj', 'slope.adj')
      par_vec.i    <- rep(1, length(adj_names.i))
      par.labels.i <- paste(paste0('rt',   resp_id.i),
                            paste0('type', shape_type.i),
                            adj_names.i, sep = "_")
      par.groups.i <- rep(paste0('rt', resp_id.i), length(par_vec.i))

      redtide_vec <- c(redtide_vec, par_vec.i)
      par.labels  <- c(par.labels, par.labels.i)
      par.groups  <- c(par.groups, par.groups.i)
    }
    par.idx <- c(par.idx, rep('redtide', length(redtide_vec)))
  } #end do.redtide
  
  if(do.disp){
    #log_disp_vec = log(disp_pars$base.val)
    disp_vec = disp_pars$base.val
    par.idx = c(par.idx,rep('disp',n_disp))
    par.labels = c(par.labels,paste('disp',disp_pars$pred,sep="_"))
    par.groups = c(par.groups,paste('disp',disp_pars$pred,sep="_"))
    #log.par.idx <- c(log.par.idx,rep(TRUE,length(log_disp_vec)))
  }
  
  med.cv = numeric()
  if(do.med){
    for(i in 1:n_med){
      #i=2
      #HYPERBOLIC SHAPE
      if(med_pars$shape_type[i]==3){
        #par_vec.i <- c(log(med_pars$x_base[i]), med_pars$par1[i], med_pars$par2[i], med_pars$par3[i])
        par_vec.i <- c(med_pars$x_base[i], med_pars$par1[i], med_pars$par2[i], med_pars$par3[i])
        par.labels.i <- paste(paste0('med',rep(med_pars$med_index[i],length(par_vec.i))), paste0('type',rep(med_pars$shape_type[i],length(par_vec.i))), c('xbase','Yzero','Yend','Ybase'),sep="_")
        par.groups.i <- paste0('med',rep(med_pars$med_index[i],length(par_vec.i)))
      }
      #LINEAR SHAPE
      if(med_pars$shape_type[i]==1){
        slope.i <- med_pars$par2[i] - med_pars$par1[i]
        #par_vec.i <- c(log(med_pars$x_base[i]), slope.i)
        par_vec.i <- c(med_pars$x_base[i], slope.i)
        par.labels.i <- paste(paste0('med',rep(med_pars$med_index[i],length(par_vec.i))), paste0('type',rep(med_pars$shape_type[i],length(par_vec.i))), c('xbase','slope'),sep="_")
        par.groups.i <- paste0('med',rep(med_pars$med_index[i],length(par_vec.i)))
      }
      med_vec <- c(med_vec,par_vec.i)
      par.labels <- c(par.labels,par.labels.i)
      par.groups <- c(par.groups,par.groups.i)
      #log.par.idx <- c(log.par.idx,c(TRUE,rep(FALSE,length(par_vec.i)-1)))
    } #end loop
    par.idx <- c(par.idx,rep('med',length(med_vec)))
  }
  
  if(do.fltdyn){
    #log_disp_vec = log(disp_pars$base.val)
    fleet_vec = fltdyn_pars$base.val
    par.idx = c(par.idx,rep('fleetdyn',length(fleet_vec)))
    tag.short = ifelse(fltdyn_pars$tag.type=='effective power','pow','mult')
    par.labels = c(par.labels,paste('flt',fltdyn_pars$fleet.idx,tag.short,sep="_"))
    par.groups = c(par.groups,paste('flt',fltdyn_pars$fleet.idx,tag.short,sep="_"))
    #log.par.idx <- c(log.par.idx,rep(TRUE,length(log_disp_vec)))
  }
  
  
  
  n_pars <<- length(vuln_vec) + length(env_vec) + length(redtide_vec) + length(disp_vec) + length(med_vec) + length(fleet_vec)
  if (n_pars == 0) stop("n_pars==0: No parameters to estimate.  Check inputs.")

  est_par_vec <- c(vuln_vec, env_vec, redtide_vec, disp_vec, med_vec, fleet_vec)
  names(est_par_vec) <- par.labels

  # index parameter types
  vul.par.idx     <<- which(par.idx=='vul')
  env.par.idx     <<- which(par.idx=='env')
  redtide.par.idx <<- which(par.idx=='redtide')
  disp.par.idx    <<- which(par.idx=='disp')
  med.par.idx     <<- which(par.idx=='med')
  flt.par.idx     <<- which(par.idx=='fleetdyn')

  # super assignment so these are available
  vul_pars     <<- vul_pars
  env_pars     <<- env_pars
  redtide_pars <<- redtide_pars
  disp_pars    <<- disp_pars
  med_pars     <<- med_pars
  fltdyn_pars  <<- fltdyn_pars
  
  # set parameter bounds--------------------------------------------------------
  ##vulnerabilities - log-scale bounds.
  ## vul.bounds.mode controls how the bounds are derived:
  ##   "formula"  (default, back-compat): log-symmetric around each base.val using
  ##              an effective log-sd of exp(2*vul.cv); the scalar vul.min/vul.max
  ##              are ignored.
  ##   "override": use [vul.min, vul.max] as the bounds for every vul parameter,
  ##              ignoring base.val and vul.cv. Useful when you want a single wide
  ##              search range (e.g., [1.01, 1e6]) rather than per-base ranges.
  lower.vuls <- upper.vuls <- numeric()
  if(do.vuls){
    if(identical(vul.bounds.mode, "override")){
      lower.vuls <- rep(vul.min, length(vul_pars$base.val))
      upper.vuls <- rep(vul.max, length(vul_pars$base.val))
    } else {
      lower.vuls <- exp(log(vul_pars$base.val) + exp(vul.cv)^2 * qnorm(0.05))
      upper.vuls <- exp(log(vul_pars$base.val) + exp(vul.cv)^2 * qnorm(0.95))
    }
    lower.vuls <- pmax(lower.vuls, 1.01)
  }
  ##environmental shapes — all supported shapes use uniform [env.min, env.max] multiplier bounds
  lower.env <- upper.env <- numeric()
  if(do.env){
  if(n_env>0){
    for(i in 1:length(env_vec)){
      label.i <- par.labels[env.par.idx[i]]
      type.i  <- as.numeric(gsub('type','', unlist(strsplit(label.i, "_"))[2]))
      if(type.i %in% c(6, 9, 10, 11)){
        lower.i <- env.min
        upper.i <- env.max
      } else {
        stop("Unsupported env response shape.idx: ", type.i,
             " (supported: 6 normal, 9 trapezoid, 10 sigmoid, 11 logistic4params).")
      }
      lower.env <- c(lower.env, lower.i)
      upper.env <- c(upper.env, upper.i)
    } #end i loop
  } #end if env
  }
  
  ##red tide — per-parameter [min, max] multiplier bounds.
  ## inflection and slope each get their own bounds; redtide.min/max are the
  ## scalar fallback for back-compat when the per-parameter overrides are NULL.
  lower.redtide <- upper.redtide <- numeric()
  if(do.redtide && length(redtide_vec)>0){
    inf_min   <- if(is.null(redtide.inf.min))   redtide.min else redtide.inf.min
    inf_max   <- if(is.null(redtide.inf.max))   redtide.max else redtide.inf.max
    slope_min <- if(is.null(redtide.slope.min)) redtide.min else redtide.slope.min
    slope_max <- if(is.null(redtide.slope.max)) redtide.max else redtide.slope.max

    rt_local_idx <- which(par.idx == 'redtide')
    rt_labels    <- par.labels[rt_local_idx]
    is_inf   <- grepl("inflection.adj", rt_labels, fixed = TRUE)
    is_slope <- grepl("slope.adj",      rt_labels, fixed = TRUE)

    lower.redtide <- ifelse(is_inf, inf_min,
                     ifelse(is_slope, slope_min, redtide.min))
    upper.redtide <- ifelse(is_inf, inf_max,
                     ifelse(is_slope, slope_max, redtide.max))
  }

  ##dispersal
  lower.disp <- upper.disp <- numeric()
  if(do.disp){
  lower.disp <- disp_vec-2*disp.cv*disp_vec
  upper.disp <- disp_vec+2*disp.cv*disp_vec
  }
  
  ##fleet dynamics -- bounds derived per-parameter so a base.val above the scalar
  ##fltdyn.max doesn't get the entire normal bell clipped against the ceiling.
  ## "adaptive" (default): bounds = base.val +/- 4 * fltdyn.cv * base.val,
  ##   widened (never shrunk) to honor [fltdyn.min, fltdyn.max] as safety bounds.
  ## "scalar" (legacy): every fleet param gets the same scalar [fltdyn.min, fltdyn.max].
  lower.flt <- upper.flt <- numeric()
  if(do.fltdyn){
    if(identical(fltdyn.bounds.mode, "scalar")){
      lower.flt <- rep(fltdyn.min, length(fleet_vec))
      upper.flt <- rep(fltdyn.max, length(fleet_vec))
    } else {
      sd_flt    <- abs(fltdyn.cv * fleet_vec)
      # lower: at least fltdyn.min (safety floor), wider only if base-4sd is even higher.
      # upper: at least fltdyn.max (safety ceiling), wider when base+4sd needs more room.
      lower.flt <- pmax(fleet_vec - 4 * sd_flt, fltdyn.min)
      upper.flt <- pmax(fleet_vec + 4 * sd_flt, fltdyn.max)
    }
  }
  ##mediation
  lower.med <- upper.med <- numeric()  
  if(do.med){
  med.cv = numeric()  
  for(i in 1:length(med_vec)){
    #i=1
    label.i <- par.labels[med.par.idx[i]]
    shape.i <- as.numeric(gsub('med','',unlist(strsplit(label.i,"_"))[1]))
    type.i <- as.numeric(gsub('type','',unlist(strsplit(label.i,"_"))[2]))
    par.i <- unlist(strsplit(label.i,"_"))[3]
    val.i <- med_vec[i]
    if(par.i=="xbase"){
      cv.i <- med.xbase.cv
      lower.i <- floor(val.i-2*val.i*med.xbase.cv)
      upper.i <- ceiling(val.i+2*val.i*med.xbase.cv)
      #cv.i = med.xbase.cv
    } else{
      cv.i <- med_pars$cv[med_pars$med_index==shape.i]
        if(val.i>0) lower.i <- val.i-2*cv.i*val.i
        if(val.i<0) lower.i <- val.i+2*cv.i*val.i
        if(val.i>0) upper.i <- val.i+2*cv.i*val.i
        if(val.i<0) upper.i <- val.i-2*cv.i*val.i
      }
    lower.med <- c(lower.med,lower.i)
    upper.med <- c(upper.med,upper.i)
    med.cv <- c(med.cv,cv.i)
  }
  }
  
  #super assignment
  est_par_vec <<- round(est_par_vec,2)
  L.bounds <<- round(c(lower.vuls, lower.env, lower.redtide, lower.disp, lower.med, lower.flt), 2)
  U.bounds <<- round(c(upper.vuls, upper.env, upper.redtide, upper.disp, upper.med, upper.flt), 2)
  n_vuls    <<- n_vuls
  n_env     <<- length(env_vec)
  n_redtide <<- length(redtide_vec)
  n_disp    <<- n_disp
  n_med     <<- length(med_vec)
  n_flt     <<- n_fleet
  par.groups <<- par.groups
  par.labels <<- par.labels

  #create vector of cv for rnorm() draws
  par_cv_vec <- est_par_vec
  par_cv_vec[vul.par.idx]     <- vul.cv
  par_cv_vec[env.par.idx]     <- env.cv
  par_cv_vec[redtide.par.idx] <- redtide.cv
  par_cv_vec[disp.par.idx]    <- disp.cv
  par_cv_vec[med.par.idx]     <- med.cv
  par_cv_vec[flt.par.idx]     <- fltdyn.cv
  par_cv_vec <<- par_cv_vec
  
} #end function

#make GA populations----

#' Generate an initial GA population matrix
#'
#' Creates the initial population for a genetic algorithm using either
#' uniform draws within parameter bounds or (optionally) normal/gamma
#' draws centered on existing parameter estimates. Uses several global
#' objects created elsewhere.
#' @param run_config List containing GA settings; must include \code{popSize}.
#' @param pardist Character string; distribution type, either"uniform" (uniform between bounds) or
#' "normal" (normal/gamma hybrid). Default "uniform"
#' @param vuldist Character string; distribution type for vulnerability parameters.
#' One of \code{"log-uniform"} (default — \code{exp(runif(log(L), log(U)))}, evenly samples
#' each order of magnitude so the sensitive low end gets equal density), \code{"uniform"}
#' (linear uniform between bounds; over-weights the high end when bounds span multiple
#' orders of magnitude), or \code{"gamma"} (shape=2, rate=2/base.val; mode near
#' \code{base/2}, mass concentrated around the base value).
#' @param rtdist Character string; distribution type for red tide parameters. Either
#'   \code{"log-uniform"} (default — symmetric in log space around the baseline of 1; appropriate
#'   when bounds span more than one order of magnitude) or \code{"uniform"} (same linear
#'   uniform as the other params, follows \code{pardist}).
#' @param fltdist Character string; distribution type for fleet dynamics parameters
#'   (effective power and effort multiplier). One of \code{"normal"} (default —
#'   \code{rnorm(mean = base.val, sd = fltdyn.cv * base.val)} centered on each parameter's
#'   base value, then clipped to \code{[L.bounds, U.bounds]}) or \code{"uniform"}
#'   (linear uniform between the bounds, follows whatever \code{pardist} produced).
#' @details
#' When red tide parameters are present, the first 6 rows of the returned matrix are
#' overwritten with explicit anchor combinations to guarantee corner-case coverage in
#' the initial population:
#' \enumerate{
#'   \item baseline (all RT adj = 1)
#'   \item gentle, early-onset (inflection.adj = min, slope.adj = min)
#'   \item steep, late-onset (inflection.adj = max, slope.adj = max)
#'   \item early crash (inflection.adj = min, slope.adj = max)
#'   \item late & gradual (inflection.adj = max, slope.adj = min)
#'   \item RT effectively off (slope.adj = min; inflection.adj as drawn)
#' }
#' Non-RT parameters in those rows keep their random draws so the anchors still vary
#' on the rest of the parameter vector.
#' @return A numeric matrix of dimension run_config$popSize x n_pars.
#' @export
fn.GApop = function(run_config=myconfig, pardist='uniform', vuldist='log-uniform',
                    rtdist='log-uniform', fltdist='normal'){
  # run_config <- myconfig
  # pardist='uniform'
  # vuldist='gamma'
  # rtdist='log-uniform'
  mat <- matrix(NA,nrow=run_config$popSize, ncol=n_pars)
  colnames(mat) <- par.labels

  #inspect the parameter set and their bounds
  data.frame(L.bounds, est_par_vec, U.bounds, par_cv_vec,par.groups)

  #population of vulnerabilities, use gamma to get more values at low end with a long tail to some higher values

  #other parameters
  if(pardist=='uniform'){
    pop.pars <- matrix(runif(run_config$popSize*(n_pars), min=L.bounds, max=U.bounds),nrow=run_config$popSize, byrow=T)
    mat[] <- pop.pars
  }
  if(pardist=='normal') {
    par_sd <- abs(par_cv_vec*est_par_vec)
    pop.pars <- matrix(rnorm(run_config$popSize*n_pars, mean=est_par_vec, sd=par_sd),nrow=run_config$popSize, byrow=T)
  }
  if(vuldist=='gamma'){
    pop.vuls <- round(matrix(rgamma(run_config$popSize*length(vul.par.idx), shape=2, rate=2/est_par_vec[vul.par.idx]),nrow=run_config$popSize, byrow=T),4)
  }
  if(vuldist=='uniform'){
    pop.vuls <- round(matrix(runif(run_config$popSize*length(vul.par.idx), min=L.bounds[vul.par.idx],max=U.bounds[vul.par.idx]),nrow=run_config$popSize, byrow=T),4)
  }
  if(vuldist=='log-uniform'){
    lo <- L.bounds[vul.par.idx]
    hi <- U.bounds[vul.par.idx]
    if(any(lo <= 0)){
      warning("vuldist='log-uniform' requires positive lower bounds; clipping to 1.01.")
      lo <- pmax(lo, 1.01)
    }
    pop.vuls <- round(matrix(exp(runif(run_config$popSize * length(vul.par.idx),
                                       min = log(lo), max = log(hi))),
                             nrow = run_config$popSize, byrow = TRUE), 4)
  }
  pop.vuls[pop.vuls<1.01] <- 1.01

  mat[] <- pop.pars
  mat[,vul.par.idx] <- pop.vuls

  # --- Red tide log-uniform sampling (option A) -------------------------------
  # Multipliers spanning multiple orders of magnitude (e.g., [0.01, 100]) sampled
  # with plain runif are biased toward the high end -- ~99% of draws end up above
  # the baseline of 1. Sampling on the log scale puts the baseline at the geometric
  # midpoint and gives symmetric coverage above and below 1.
  if(rtdist == 'log-uniform' && length(redtide.par.idx) > 0L){
    lo <- L.bounds[redtide.par.idx]
    hi <- U.bounds[redtide.par.idx]
    if(any(lo <= 0)){
      warning("rtdist='log-uniform' requires positive lower bounds; clipping to 1e-6.")
      lo <- pmax(lo, 1e-6)
    }
    rt_draws <- matrix(exp(runif(run_config$popSize * length(redtide.par.idx),
                                 min = log(lo), max = log(hi))),
                       nrow = run_config$popSize, byrow = TRUE)
    mat[, redtide.par.idx] <- rt_draws
  }

  # --- Fleet dynamics: normal around base value -------------------------------
  # Independent of pardist so the fleet shape stays normal-around-base regardless
  # of what the other parameters do. par_cv_vec[flt.par.idx] is the user's
  # fltdyn.cv; sd = fltdyn.cv * base.val. The downstream pmin/pmax clip enforces
  # the declared [fltdyn.min, fltdyn.max] bounds.
  if(fltdist == 'normal' && length(flt.par.idx) > 0L){
    mean_flt <- est_par_vec[flt.par.idx]
    sd_flt   <- abs(par_cv_vec[flt.par.idx] * mean_flt)
    flt_draws <- matrix(rnorm(run_config$popSize * length(flt.par.idx),
                              mean = mean_flt, sd = sd_flt),
                        nrow = run_config$popSize, byrow = TRUE)
    mat[, flt.par.idx] <- flt_draws
  }

  mat <- round(mat,4)

  integer.idx <- grep("xbase",par.labels)
  if(length(integer.idx)>0) mat[,integer.idx] <- round(mat[,integer.idx])
  # Clip fleet and env draws to their declared per-parameter bounds. Naively
  # passing a length-k bounds vector to pmin/pmax(matrix, vec) silently recycles
  # element-wise through the column-major matrix (so cell (r, c) gets bounds at
  # an index that depends on r AND c, not just c). When bounds are heterogeneous
  # across columns -- as they are under fltdyn.bounds.mode = "adaptive" -- this
  # scrambles the per-column bound assignment and produces spurious pile-ups at
  # the wrong bound. Build full bounds matrices to enforce per-column clipping.
  if(length(flt.par.idx) > 0L){
    N    <- nrow(mat)
    Lmat <- matrix(rep(L.bounds[flt.par.idx], each = N), nrow = N)
    Umat <- matrix(rep(U.bounds[flt.par.idx], each = N), nrow = N)
    mat[, flt.par.idx] <- pmin(pmax(mat[, flt.par.idx], Lmat), Umat)
  }
  if(length(env.par.idx) > 0L){
    N    <- nrow(mat)
    Lmat <- matrix(rep(L.bounds[env.par.idx], each = N), nrow = N)
    Umat <- matrix(rep(U.bounds[env.par.idx], each = N), nrow = N)
    mat[, env.par.idx] <- pmin(pmax(mat[, env.par.idx], Lmat), Umat)
  }

  # --- Red tide anchor seeding (option D) -------------------------------------
  # Inject explicit corner-case combinations into the first rows so the GA explores
  # baseline, gentle, steep, early-crash, late-gradual, and RT-off regimes from the
  # start instead of drifting there.
  if(length(redtide.par.idx) > 0L && run_config$popSize >= 6L){
    rt_labels <- par.labels[redtide.par.idx]
    is_inf    <- grepl("inflection.adj", rt_labels, fixed = TRUE)
    is_slope  <- grepl("slope.adj",      rt_labels, fixed = TRUE)
    inf_cols   <- redtide.par.idx[is_inf]
    slope_cols <- redtide.par.idx[is_slope]

    if(length(inf_cols) > 0L && length(slope_cols) > 0L){
      inf_lo <- L.bounds[inf_cols]; inf_hi <- U.bounds[inf_cols]
      slo_lo <- L.bounds[slope_cols]; slo_hi <- U.bounds[slope_cols]

      # 1: baseline (no perturbation)
      mat[1, redtide.par.idx] <- 1
      # 2: gentle, early-onset
      mat[2, inf_cols]   <- inf_lo; mat[2, slope_cols] <- slo_lo
      # 3: steep, late-onset
      mat[3, inf_cols]   <- inf_hi; mat[3, slope_cols] <- slo_hi
      # 4: early crash
      mat[4, inf_cols]   <- inf_lo; mat[4, slope_cols] <- slo_hi
      # 5: late & gradual
      mat[5, inf_cols]   <- inf_hi; mat[5, slope_cols] <- slo_lo
      # 6: RT effectively off (slope at min; inflection left as drawn)
      mat[6, slope_cols] <- slo_lo

      mat[1:6, ] <- round(mat[1:6, ], 4)
    }
  }

  # Note: the base run is included in the initial population by fn.GA (which
  # overwrites gapop[1, ] with est_par_vec after this function returns), so we
  # do not seed a baseline row here. All rows returned from this function are
  # random draws (subject only to the red-tide corner-case anchors in rows 1-6
  # if the RT block above fired).

  return(mat)
} #eof


#' @keywords internal
#' @noRd
safe_make_run_dir_tempfile <- function(base_dir,
                                       par_vec,
                                       gen = gen,
                                       idx = idx,
                                       prefix = "run",
                                       random_len = 8) {

  
  # Create a compact, reproducible tag for the parameter vector
  # Use a hash to avoid very long names; include gen/idx if provided
  par_tag <- tryCatch({
    # stringify with limited precision to keep it small and stable
    s <- paste(round(par_vec, 6), collapse = "_")
    # hash the string; if digest not available, fallback to simple checksum
    if (requireNamespace("digest", quietly = TRUE)) {
      digest::digest(s, algo = "xxhash64")
    } else {
      # simple fallback: integer checksum
      sprintf("%08X", as.integer(sum(charToRaw(s))))
    }
  }, error = function(e) {
    sprintf("%08X", as.integer(runif(1, 0, .Machine$integer.max)))
  })
  
  # Add generation/index if provided
  gid <- if (!is.null(gen)) sprintf("g%03d", as.integer(gen)) else NULL
  iid <- if (!is.null(idx)) sprintf("i%04d", as.integer(idx)) else NULL
  
  # Add process id for extra uniqueness
  pid <- tryCatch(Sys.getpid(), error = function(e) sample.int(1e6, 1))
  
  # Random slug to avoid race collisions (tempfile-like behavior)
  rand <- paste(sample(c(letters, LETTERS, 0:9), random_len, replace = TRUE), collapse = "")
  
  # Assemble a safe folder name: prefix + tags
  parts <- c(prefix, gid, iid, paste0("pid", pid), paste0("h", par_tag), rand)
  folder_name <- paste(parts[!sapply(parts, is.null)], collapse = "_")
  
  # Sanitize (just in case)
  folder_name <- gsub("[^A-Za-z0-9_.-]", "_", folder_name)
  
  run_path <- file.path(base_dir, folder_name)
  
  # Create (retry if rare collision)
  if (!dir.exists(run_path)) {
    ok <- dir.create(run_path, recursive = TRUE, showWarnings = FALSE)
    if (!ok) {
      # one more try with a new random slug
      rand2 <- paste(sample(c(letters, LETTERS, 0:9), random_len, replace = TRUE), collapse = "")
      run_path <- file.path(base_dir, paste0(folder_name, "_", rand2))
      dir.create(run_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(run_path)
}



#write command files----

#' @title Write command file from parameter vector
#' @description Takes as input a parameter vector and creates taglines and command files.
#' @param par_vec Parameter vector.
#' @param g Generation number, used for file tracking in a genetic algorithm loop. Default 0
#' @param idx Index for run identifier.
#' @param out_dir Output directory, where the command file will be saved along with model output.#'
#' @keywords internal
#' @noRd
.gen_dir <- function(base_run_dir, g)
  file.path(base_run_dir, sprintf("gen%02d", as.integer(g)))

#' @keywords internal
#' @noRd
#' Move elite run directories forward into the next generation's folder,
#' preserving their original folder names (which encode gen-of-origin +
#' offspring index for traceability), and rewriting the
#' <ECOSPACE_OUTPUT_DIR> line inside each moved cmd.txt so the cmd file
#' remains re-runnable from its new location.
#'
#' @param elite_cmd_paths Character vector of cmd.txt paths for the elites
#'   in their CURRENT (source) generation folder.
#' @param dst_gen_dir Destination generation folder (created if missing).
#' @return Character vector of new cmd.txt paths, aligned with input.
.move_elite_dirs <- function(elite_cmd_paths, dst_gen_dir){
  if(!dir.exists(dst_gen_dir))
    dir.create(dst_gen_dir, recursive = TRUE, showWarnings = FALSE)

  new_paths <- character(length(elite_cmd_paths))
  for(k in seq_along(elite_cmd_paths)){
    old_dir <- dirname(elite_cmd_paths[k])
    if(!dir.exists(old_dir)){
      warning(sprintf("Elite dir missing on disk: %s (leaving cmd path unchanged)", old_dir))
      new_paths[k] <- elite_cmd_paths[k]
      next
    }
    new_dir <- file.path(dst_gen_dir, basename(old_dir))
    if(dir.exists(new_dir)){
      # Same-name dir already in destination — shouldn't happen given the
      # random-hash suffix on run folders, but handle defensively.
      warning(sprintf("Destination exists (leaving in place): %s", new_dir))
      new_paths[k] <- elite_cmd_paths[k]
      next
    }
    moved <- suppressWarnings(file.rename(old_dir, new_dir))
    if(!moved){
      # Cross-device fallback (rare — gens should share one filesystem).
      dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)
      copied <- file.copy(list.files(old_dir, full.names = TRUE, recursive = FALSE),
                          new_dir, overwrite = TRUE, recursive = TRUE)
      if(all(copied)){
        unlink(old_dir, recursive = TRUE)
      } else {
        stop(sprintf("Failed to move elite dir from %s to %s", old_dir, new_dir))
      }
    }
    new_cmd <- file.path(new_dir, "cmd.txt")
    new_paths[k] <- new_cmd
    if(file.exists(new_cmd)){
      L <- readLines(new_cmd)
      new_out <- normalizePath(new_dir, winslash = "\\", mustWork = FALSE)
      out_line <- startsWith(L, "<ECOSPACE_OUTPUT_DIR>")
      if(any(out_line))
        L[out_line] <- sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", new_out)
      writeLines(L, new_cmd)
    }
  }
  new_paths
}

#' @return Path to the written command file.
#' @export
fn.parvec2cmd <- function(par_vec=est_par_vec, g=0, idx=0,
                          out_dir=.gen_dir(run_dir, g)){
  if(!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # par_vec = est_par_vec
  # g = 999
  # idx = 0
  # out_dir=run_dir
  
  #output directory----
  run_path <- safe_make_run_dir_tempfile(base_dir=out_dir,
                                         par_vec=par_vec,
                                         gen = g,
                                         idx = idx,
                                         prefix = "run",
                                         random_len = 8)
  
  #parameter tags----
  tags.vul <- tags.env <- tags.redtide <- tags.disp <- tags.med <- tags.flt <- character()
  if(length(vul.par.idx)>0){
    vuln_vec <- par_vec[vul.par.idx]
    tags.vul <- paste0("<ECOSIM_VULNERABILITIES_INDEXED>(", vul_pars$pred, " ", ifelse(is.na(vul_pars$prey),"",vul_pars$prey),
                       "), ", sprintf("%.2f", vuln_vec), ", Indexed.Single")
  }

  if(length(disp.par.idx)>0){
    disp_vec <- par_vec[disp.par.idx]
    tags.disp <- paste0("<ECOSPACE_DISPERSAL_RATE_INDEXED>(", disp_pars$pred,"), ",disp_vec, ", Indexed.Single")
  }
  
  if(length(env.par.idx)>0){
    env_groups <- unique(par.groups[env.par.idx])
    env_idx_vec <- .env_resp_idx(env_pars)         # accepts fxn_num or resp.idx
    tags.env <- character()
    for(i in 1:length(env_groups)){
      par_set.i    <- par_vec[which(par.groups == env_groups[i])]
      par_labels.i <- par.labels[which(par.groups == env_groups[i])]
      type.i       <- as.numeric(gsub('type','', strsplit(par_labels.i[1], "_")[[1]][2]))
      resp_id.i    <- as.numeric(gsub('env','',  strsplit(par_labels.i[1], "_")[[1]][1]))
      row.i        <- which(env_idx_vec == resp_id.i)[1]
      if(is.na(row.i))
        stop("fn.parvec2cmd: no env_pars row matches response index ", resp_id.i)

      if(type.i == 9){
        ##TRAPEZOID — par1=LB, par2=LT, par3=RT, par4=RB
        ## GA pars: mid.adj, width.adj
        LB.base   <- env_pars[row.i, 'par1']
        LT.base   <- env_pars[row.i, 'par2']
        RT.base   <- env_pars[row.i, 'par3']
        RB.base   <- env_pars[row.i, 'par4']
        mid.base  <- (RT.base + LT.base) / 2
        pref.base <- RT.base - LT.base
        mid.i     <- mid.base  * par_set.i[1]      # mid.adj
        pref.i    <- pref.base * par_set.i[2]      # width.adj
        LT.i      <- round(pmax(0, mid.i - pref.i / 2), 2)
        RT.i      <- round(           mid.i + pref.i / 2,  2)
        LB.i      <- round(LT.i - (LT.base - LB.base) * par_set.i[2], 2)
        RB.i      <- round(RT.i + (RB.base - RT.base) * par_set.i[2], 2)
        par.string.i <- paste(type.i, LB.i, LT.i, RT.i, RB.i, sep = " ")

      } else if(type.i == 6){
        ##NORMAL — par1=SDLeft, par2=DataWidth(auto), par3=SDRight, par4=Mean, par5=Max
        ## GA pars: mean.adj (multiplies Mean), width.adj (multiplies SDLeft and SDRight jointly)
        SDLeft.base    <- env_pars[row.i, 'par1']
        SDRight.base   <- env_pars[row.i, 'par3']
        Mean.base      <- env_pars[row.i, 'par4']
        Max.base       <- env_pars[row.i, 'par5']
        Mean.i         <- round(Mean.base    * par_set.i[1], 4)    # mean.adj
        SDLeft.i       <- round(SDLeft.base  * par_set.i[2], 4)    # width.adj (symmetric)
        SDRight.i      <- round(SDRight.base * par_set.i[2], 4)
        DataWidth.i    <- round(5 * SDLeft.i + 5 * SDRight.i, 4)   # auto per cNormalShapeFunction.vb
        par.string.i   <- paste(type.i, SDLeft.i, DataWidth.i, SDRight.i, Mean.i, Max.base, sep = " ")

      } else if(type.i == 10){
        ##SIGMOID — par1=XMin, par2=XMax, par3=XMid, par4=XOpt, par5=Slope, par6=Scalar
        ## GA pars: xmid.adj (multiplies XMid), slope.adj (multiplies Slope)
        XMin.base    <- env_pars[row.i, 'par1']
        XMax.base    <- env_pars[row.i, 'par2']
        XMid.base    <- env_pars[row.i, 'par3']
        XOpt.base    <- env_pars[row.i, 'par4']
        Slope.base   <- env_pars[row.i, 'par5']
        Scalar.base  <- env_pars[row.i, 'par6']
        XMid.i       <- round(XMid.base  * par_set.i[1], 4)        # xmid.adj
        Slope.i      <- round(Slope.base * par_set.i[2], 4)        # slope.adj
        par.string.i <- paste(type.i, XMin.base, XMax.base, XMid.i, XOpt.base, Slope.i, Scalar.base, sep = " ")

      } else if(type.i == 11){
        ##LOGISTIC_4PARAMS — par1=XMin, par2=XMax, par3=Inflection, par4=Slope
        ## GA pars: inflection.adj (multiplies Inflection), slope.adj (multiplies Slope)
        XMin.base    <- env_pars[row.i, 'par1']
        XMax.base    <- env_pars[row.i, 'par2']
        Inf.base     <- env_pars[row.i, 'par3']
        Slope.base   <- env_pars[row.i, 'par4']
        Inf.i        <- round(Inf.base   * par_set.i[1], 4)        # inflection.adj
        Slope.i      <- round(Slope.base * par_set.i[2], 4)        # slope.adj
        par.string.i <- paste(type.i, XMin.base, XMax.base, Inf.i, Slope.i, sep = " ")

      } else {
        stop("fn.parvec2cmd: unsupported env response shape.idx: ", type.i,
             " (supported: 6 normal, 9 trapezoid, 10 sigmoid, 11 logistic4params).")
      }

      tag.i <- paste0("<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(", resp_id.i, "), ",
                      par.string.i, ", Indexed.Single[]")
      tags.env <- c(tags.env, tag.i)
    }
  }

  if(length(redtide.par.idx)>0){
    redtide_pars   <- .normalize_env_pars(redtide_pars)
    redtide_groups <- unique(par.groups[redtide.par.idx])
    rt_idx_vec     <- .env_resp_idx(redtide_pars)
    for(i in 1:length(redtide_groups)){
      par_set.i    <- par_vec[which(par.groups == redtide_groups[i])]
      par_labels.i <- par.labels[which(par.groups == redtide_groups[i])]
      type.i       <- as.numeric(gsub('type','', strsplit(par_labels.i[1], "_")[[1]][2]))
      resp_id.i    <- as.numeric(gsub('rt','',   strsplit(par_labels.i[1], "_")[[1]][1]))
      row.i        <- which(rt_idx_vec == resp_id.i)[1]
      if(is.na(row.i))
        stop("fn.parvec2cmd: no redtide_pars row matches response index ", resp_id.i)

      # All red tide responses are shape 11 (logistic4params).
      # par1=XMin, par2=XMax, par3=Inflection, par4=Slope
      # GA pars: inflection.adj, slope.adj
      XMin.base  <- redtide_pars[row.i, 'par1']
      XMax.base  <- redtide_pars[row.i, 'par2']
      Inf.base   <- redtide_pars[row.i, 'par3']
      Slope.base <- redtide_pars[row.i, 'par4']
      Inf.i      <- round(Inf.base   * par_set.i[1], 4)    # inflection.adj
      Slope.i    <- round(Slope.base * par_set.i[2], 4)    # slope.adj
      par.string.i <- paste(type.i, XMin.base, XMax.base, Inf.i, Slope.i, sep = " ")
      tag.i <- paste0("<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(", resp_id.i, "), ",
                      par.string.i, ", Indexed.Single[]")
      tags.redtide <- c(tags.redtide, tag.i)
    }
  }

  if(length(med.par.idx)>0){
    med_groups <- unique(par.groups[med.par.idx])
    tags.med <- character()
    for(i in 1:length(med_groups)){
      #i=1
      par_set.i <- par_vec[which(par.groups==med_groups[i])]
      par_labels.i <- par.labels[which(par.groups==med_groups[i])]
      type.i <- as.numeric(gsub('type','',strsplit(par_labels.i[1],"_")[[1]][2]))
      shape.i <- as.numeric(gsub('med','',strsplit(par_labels.i[1],"_")[[1]][1]))
      #LINEAR SHAPE
      if(type.i==1){
        par.string.i <- paste(par_set.i[1], type.i, 1-0.5*par_set.i[2], 1+0.5*par_set.i[2], sep=" ")
        tag.i <- paste0("<MEDIATION_FUNCTION_INDEXED>(",shape.i,"), ",par.string.i,", Indexed.Single[]")
      }
      #HYPERBOLIC SHAPE
      if(type.i==3){
        par.string.i <- paste(par_set.i[1], type.i, par_set.i[2], par_set.i[3], par_set.i[4], 1, sep=" ")
        tag.i <- paste0("<MEDIATION_FUNCTION_INDEXED>(",shape.i,"), ",par.string.i,", Indexed.Single[]")
      }
      
      tags.med <- c(tags.med,tag.i)
    }
  }
  
  if(length(flt.par.idx)>0){
    flt_vec <- par_vec[flt.par.idx]
    flt.par.type = sapply(strsplit(names(flt_vec),"_"),"[",3)
    flt.par.type = ifelse(flt.par.type=='pow','EFFECTIVEPOWER','TOTAL_EFFORTMULT')
    fleetnum = sapply(strsplit(names(flt_vec),"_"),"[",2)
    tags.flt <- paste0("<ECOSPACE_FLEET_",flt.par.type,"_INDEXED>(", fleetnum,"), ",flt_vec, ", Indexed.Single")
  }
  
  
  tags = c(tags.vul, tags.env, tags.redtide, tags.disp, tags.med, tags.flt)
  
  
  #command files----
  # EwE silently mis-parses mixed-separator paths (e.g. "C:\WFS MICE\foo/bar")
  # on Windows -- the forward slashes get treated as filename characters, so
  # output ends up in a flat sibling directory with dashes instead of slashes.
  # Force native Windows separators so the user can construct dir.out.ga any way
  # they like without this trap.
  output_dir_for_ewe <- normalizePath(run_path, winslash = "\\", mustWork = FALSE)
  cmd_j <- cmd_base
  cmd_j[startsWith(cmd_j, "<ECOSPACE_OUTPUT_DIR>")] <-
    sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", output_dir_for_ewe)
  cmd_j <- c(cmd_j, tags)
  cmd_file <- file.path(run_path, "cmd.txt")
  writeLines(cmd_j, cmd_file)
  return(cmd_file)
  
}

#' @title Plot GA population: response-function curves for env/redtide, histograms otherwise.
#' @description For each unique parameter group, plots a panel:
#'   \itemize{
#'     \item Env response groups (\code{env<N>_type<S>}) and red-tide groups (\code{rt<N>_type<S>}):
#'           overlay one response curve per population member, with a highlighted curve
#'           drawn in blue. By default that curve is the baseline; pass \code{highlight_pars}
#'           to highlight a different parameter set (e.g. the min-LL individual from the
#'           final generation).
#'     \item All other groups (vul, disp, fleetdyn, med): one histogram per parameter with
#'           a vertical dashed line at the highlighted value.
#'   }
#'   Writes the panels to a multipage PDF.
#' @param gapop Numeric matrix \code{[popSize x n_pars]} of parameter values.
#' @param est_par_vec Named numeric vector of baseline parameter values.
#' @param par.groups Character vector mapping each column of \code{gapop} to a group.
#' @param par.labels Character vector of column labels (typically the colnames of \code{gapop}).
#' @param env_pars Data frame of env response baselines (canonical or sensitivity-results
#'   schema, see \code{\link{fn.makeparvec}}). May be \code{NULL}.
#' @param redtide_pars Same as \code{env_pars} for red-tide responses. May be \code{NULL}.
#' @param file Output PDF path.
#' @param dims Panel grid \code{c(rows, cols)}. Default \code{c(3, 3)}.
#' @param highlight_pars Numeric vector of length \code{ncol(gapop)} drawn in BLUE
#'   (dashed vline on histograms; overlaid response curve on env/rt panels). If
#'   \code{NULL} (default — initial-pop plot use case), only the red baseline is drawn.
#'   Pass the min-LL individual (e.g. \code{gapop[which.min(fitness), ]}) for the
#'   final-population plot to overlay best-fit on top of the baseline + ensemble.
#' @section Markers:
#'   RED: \code{est_par_vec} — baseline / starting values; on env/rt panels the red
#'   curve uses adj = 1 (i.e. the natural EwE curve).
#'   BLUE: \code{highlight_pars} (only when supplied) — typically the gen-N best-fit.
#' @return Invisibly returns \code{file}.
#' @export
fn.plot_ga_pop_responses <- function(gapop, est_par_vec, par.groups, par.labels,
                                     env_pars = NULL, redtide_pars = NULL,
                                     file, dims = c(3, 3),
                                     highlight_pars = NULL){
  if(!is.null(env_pars))     env_pars     <- .normalize_env_pars(env_pars)
  if(!is.null(redtide_pars)) redtide_pars <- .normalize_env_pars(redtide_pars)

  # Defensive coercion: callers may pass highlight_pars / est_par_vec as a 1-row
  # data.frame (e.g. gapop.final[which.min(fitness), ] when gapop.final is a
  # data.frame). Indexing such an object with [j] returns a 1-col data.frame
  # instead of a scalar, which silently breaks the env/rt response-curve math
  # (numeric_vec * 1-col-df preserves the data.frame class and yields a length-1
  # column, mismatching x in the lines() call).
  if(!is.null(highlight_pars)) highlight_pars <- unlist(highlight_pars, use.names = TRUE)
  if(!is.null(est_par_vec))    est_par_vec    <- unlist(est_par_vec,    use.names = TRUE)

  # Each panel gets up to two reference markers:
  #   RED  = est_par_vec (the baseline / starting values)
  #   BLUE = highlight_pars (e.g. the min-LL individual from the final gen)
  # Initial-pop plots pass highlight_pars = NULL -- only the red baseline is drawn.
  # Final-pop plots pass highlight_pars = best-fit individual -- both red+blue drawn.
  show_blue <- !is.null(highlight_pars)
  hl <- if(show_blue) highlight_pars else est_par_vec

  alpha_col <- grDevices::rgb(0, 0, 0, alpha = 0.15)

  # Small inline legend for each panel. Line style matches how the marker is
  # drawn: dashed vline on histograms, solid line on response curves. Only
  # includes the blue "best-fit" entry when highlight_pars was supplied.
  .add_legend <- function(kind = c("hist", "curve"),
                          where = "topright"){
    kind <- match.arg(kind)
    lty  <- if(kind == "hist") 2   else 1
    lwd  <- if(kind == "hist") 1.5 else 2
    if(show_blue){
      graphics::legend(where,
                       legend = c("baseline (est_par_vec)", "best-fit"),
                       col = c("red", "blue"), lty = lty, lwd = lwd,
                       bty = "o", box.col = "gray70", bg = "white",
                       cex = 0.65, inset = 0.01)
    } else {
      graphics::legend(where,
                       legend = "baseline (est_par_vec)",
                       col = "red", lty = lty, lwd = lwd,
                       bty = "o", box.col = "gray70", bg = "white",
                       cex = 0.65, inset = 0.01)
    }
  }

  grDevices::pdf(file, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = dims, mar = c(3, 3, 2, 1), mgp = c(1.6, 0.5, 0))

  draw_hists <- function(col_idx){
    for(j in col_idx){
      graphics::hist(gapop[, j], breaks = 50, main = par.labels[j], xlab = 'value')
      graphics::abline(v = est_par_vec[j], col = 'red',  lty = 2, lwd = 1.5)
      if(show_blue)
        graphics::abline(v = hl[j],        col = 'blue', lty = 2, lwd = 1.5)
      .add_legend("hist")
    }
  }

  unique_groups <- unique(par.groups)
  for(grp in unique_groups){
    col_idx <- which(par.groups == grp)
    is_env  <- startsWith(grp, "env")
    is_rt   <- startsWith(grp, "rt")
    pars_df <- if(is_rt) redtide_pars else if(is_env) env_pars else NULL

    if(is.null(pars_df)){ draw_hists(col_idx); next }

    parts      <- strsplit(par.labels[col_idx[1]], "_")[[1]]
    shape_type <- suppressWarnings(as.numeric(gsub("type", "", parts[2])))
    resp_id    <- suppressWarnings(as.numeric(gsub("^(env|rt)", "", parts[1])))
    row.i      <- which(.env_resp_idx(pars_df) == resp_id)[1]
    if(is.na(row.i) || is.na(shape_type)){ draw_hists(col_idx); next }

    adj1 <- gapop[, col_idx[1]]
    adj2 <- if(length(col_idx) >= 2) gapop[, col_idx[2]] else rep(1, nrow(gapop))

    # Baseline curve uses adj=1 (env/rt GA pars are MULTIPLIERS on the natural curve),
    # so red curve == natural EwE curve. Blue curve uses highlight_pars when supplied.
    h_adj1 <- if(show_blue) highlight_pars[col_idx[1]] else 1
    h_adj2 <- if(show_blue && length(col_idx) >= 2) highlight_pars[col_idx[2]] else 1

    if(shape_type == 11){
      # Logistic 4-params: par1=XMin, par2=XMax, par3=Inflection, par4=Slope
      XMin       <- pars_df[row.i, 'par1']
      XMax       <- pars_df[row.i, 'par2']
      Inf.base   <- pars_df[row.i, 'par3']
      Slope.base <- pars_df[row.i, 'par4']
      x <- seq(XMin, XMax, length.out = 200)
      ycurve <- function(Inf.i, Slope.i) 1 - 1 / (1 + (x / Inf.i)^Slope.i)
      graphics::plot(NA, xlim = c(XMin, XMax), ylim = c(0, 1.05),
                     main = grp, xlab = "x", ylab = "response")
      for(k in seq_len(nrow(gapop)))
        graphics::lines(x, ycurve(Inf.base * adj1[k], Slope.base * adj2[k]), col = alpha_col)
      graphics::lines(x, ycurve(Inf.base, Slope.base), col = 'red',  lwd = 2)
      if(show_blue)
        graphics::lines(x, ycurve(Inf.base * h_adj1, Slope.base * h_adj2), col = 'blue', lwd = 2)
      .add_legend("curve", where = "topleft")

    } else if(shape_type == 9){
      # Trapezoid: par1=LB, par2=LT, par3=RT, par4=RB. GA: mid.adj, width.adj.
      LB <- pars_df[row.i, 'par1']; LT <- pars_df[row.i, 'par2']
      RT <- pars_df[row.i, 'par3']; RB <- pars_df[row.i, 'par4']
      mid.base  <- (RT + LT) / 2
      pref.base <- RT - LT
      xlim <- range(c(LB, LT, RT, RB), na.rm = TRUE)
      xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
      graphics::plot(NA, xlim = xlim, ylim = c(0, 1.05),
                     main = grp, xlab = "x", ylab = "response")
      trap_pts <- function(mid.adj, w.adj){
        mid.i  <- mid.base  * mid.adj
        pref.i <- pref.base * w.adj
        LT.i <- pmax(0, mid.i - pref.i/2); RT.i <- mid.i + pref.i/2
        LB.i <- LT.i - (LT - LB) * w.adj;  RB.i <- RT.i + (RB - RT) * w.adj
        list(x = c(LB.i, LT.i, RT.i, RB.i), y = c(0, 1, 1, 0))
      }
      for(k in seq_len(nrow(gapop))){
        p <- trap_pts(adj1[k], adj2[k]); graphics::lines(p$x, p$y, col = alpha_col)
      }
      rp <- trap_pts(1, 1)
      graphics::lines(rp$x, rp$y, col = 'red', lwd = 2)
      if(show_blue){
        hp <- trap_pts(h_adj1, h_adj2)
        graphics::lines(hp$x, hp$y, col = 'blue', lwd = 2)
      }
      .add_legend("curve", where = "bottomright")

    } else if(shape_type == 6){
      # Normal: par1=SDLeft, par3=SDRight, par4=Mean, par5=Max. GA: mean.adj, width.adj.
      SDLeft.base  <- pars_df[row.i, 'par1']
      SDRight.base <- pars_df[row.i, 'par3']
      Mean.base    <- pars_df[row.i, 'par4']
      Max.base     <- pars_df[row.i, 'par5']
      DataWidth.base <- 5 * SDLeft.base + 5 * SDRight.base
      x <- seq(Mean.base - DataWidth.base/2, Mean.base + DataWidth.base/2, length.out = 200)
      norm_y <- function(Mean.i, SDL, SDR, Mx){
        y <- numeric(length(x)); lf <- x < Mean.i
        y[lf]  <- Mx * exp(-0.5 * ((x[lf]  - Mean.i) / SDL)^2)
        y[!lf] <- Mx * exp(-0.5 * ((x[!lf] - Mean.i) / SDR)^2)
        y
      }
      graphics::plot(NA, xlim = range(x), ylim = c(0, Max.base * 1.2),
                     main = grp, xlab = "x", ylab = "response")
      for(k in seq_len(nrow(gapop)))
        graphics::lines(x, norm_y(Mean.base * adj1[k],
                                  SDLeft.base * adj2[k],
                                  SDRight.base * adj2[k],
                                  Max.base), col = alpha_col)
      graphics::lines(x, norm_y(Mean.base, SDLeft.base, SDRight.base, Max.base),
                      col = 'red',  lwd = 2)
      if(show_blue)
        graphics::lines(x, norm_y(Mean.base * h_adj1,
                                  SDLeft.base * h_adj2,
                                  SDRight.base * h_adj2,
                                  Max.base),
                        col = 'blue', lwd = 2)
      .add_legend("curve", where = "topright")

    } else if(shape_type == 10){
      # Sigmoid: par1=XMin, par2=XMax, par3=XMid, par5=Slope, par6=Scalar. GA: xmid.adj, slope.adj.
      XMin        <- pars_df[row.i, 'par1']
      XMax        <- pars_df[row.i, 'par2']
      XMid.base   <- pars_df[row.i, 'par3']
      Slope.base  <- pars_df[row.i, 'par5']
      Scalar.base <- pars_df[row.i, 'par6']
      x <- seq(XMin, XMax, length.out = 200)
      ycurve <- function(XMid.i, Slope.i) Scalar.base / (1 + exp(-Slope.i * (x - XMid.i)))
      graphics::plot(NA, xlim = c(XMin, XMax), ylim = c(0, Scalar.base * 1.1),
                     main = grp, xlab = "x", ylab = "response")
      for(k in seq_len(nrow(gapop)))
        graphics::lines(x, ycurve(XMid.base * adj1[k], Slope.base * adj2[k]), col = alpha_col)
      graphics::lines(x, ycurve(XMid.base, Slope.base), col = 'red', lwd = 2)
      if(show_blue)
        graphics::lines(x, ycurve(XMid.base * h_adj1, Slope.base * h_adj2), col = 'blue', lwd = 2)
      .add_legend("curve", where = "topleft")

    } else {
      draw_hists(col_idx)
    }
  }
  invisible(file)
}


#' @keywords internal
#' @noRd
safe_runEwE <- function(cmdfile, do.obj, fit.abs.catch = TRUE) {
  on.exit(gc(), add = TRUE)   # ensure cleanup
  out <- tryCatch(
    fn.runEwE(cmdfile = cmdfile, do.obj = do.obj, fit.abs.catch = fit.abs.catch),
    error = function(e) NA_real_)
  if (is.numeric(out) && length(out) >= 1 && is.finite(out[1])){
    # Reject artificially low scores from Ecospace runs that aborted mid-write
    # (only Region_0 CSVs on disk). fn.ecospace_objfxn quietly skips obs series
    # whose region file is missing, so the returned fitness is finite but
    # systematically low. When `expected_regions` is defined in .GlobalEnv
    # (e.g. 0:16 for phase-3 WFS MICE), require every listed region's annual
    # biomass CSV to be on disk with real content; otherwise return NA so the
    # gapop retry loop swaps in a fresh random parvec.
    expected_regions <- if(exists("expected_regions", envir = .GlobalEnv))
                          get("expected_regions", envir = .GlobalEnv) else NULL
    if(!.check_ecospace_output(dirname(cmdfile), expected_regions))
      return(NA_real_)
    out[1]
  } else {
    NA_real_
  }
}#eof


#run the population of models----
#' @title Run Ecospace in parallel, for GA calibration.
#' @description
#' Executes Ecospace runs from a vector of command file paths using parallel
#' workers, retries any failed runs up to a limit, optionally deletes outputs,
#' and returns the fitness (e.g., log-likelihood) vector.
#'
#' @param files.cmd Character vector of paths to Ecospace command files.
#' @param obj.fxn Integer or flag selecting which objective function to use.
#' @param delete.output Logical; if true, delete generated output folders after runs.
#' @param max_tries Maximum number of retry rounds for failed runs.
#'
#' @details
#' Uses foreach with a parallel backend to run fn_runEwE on each command file.
#' Failed runs (NA results) are identified and resubmitted up to max_tries times.
#' If delete.output is true and run_dir exists in the environment, subdirectories
#' under run_dir are removed after completion.
#'
#' @return Numeric vector of fitness values (same length as files.cmd).
#'

fn.runEwE.gapop <-  function(
    files.cmd,
    obj.fxn=1,
    fit.abs.catch=TRUE,
    delete.output=T,
    max_tries=50,
    gen=NULL,
    pardist = if(exists("gapop.pardist", envir=.GlobalEnv)) get("gapop.pardist", envir=.GlobalEnv) else "uniform",
    vuldist = if(exists("gapop.vuldist", envir=.GlobalEnv)) get("gapop.vuldist", envir=.GlobalEnv) else "uniform",
    rtdist  = if(exists("gapop.rtdist",  envir=.GlobalEnv)) get("gapop.rtdist",  envir=.GlobalEnv) else "log-uniform",
    fltdist = if(exists("gapop.fltdist", envir=.GlobalEnv)) get("gapop.fltdist", envir=.GlobalEnv) else "normal"
){
  t1 <- Sys.time()

  # Windowed progress bar so the per-generation summary lines remain readable in
  # the console. Title carries the generation tag; label shows the live count.
  # Requires doSNOW backend (see ensure_cluster).
  gen_label <- if(is.null(gen)) "GA pop" else sprintf("Gen %s", as.character(gen))
  n_runs    <- length(files.cmd)
  pbar <- utils::winProgressBar(
    title = sprintf("Ecospace GA - %s  (%d runs)", gen_label, n_runs),
    label = sprintf("0 / %d completed (0%%)", n_runs),
    min   = 0, max = n_runs, initial = 0, width = 420)
  on.exit(try(close(pbar), silent = TRUE), add = TRUE)
  prog <- function(i){
    utils::setWinProgressBar(pbar, i,
      label = sprintf("%d / %d completed (%.1f%%)", i, n_runs, 100 * i / n_runs))
  }
  opts <- list(progress = prog)

  runLL.out <- foreach(i = seq_along(files.cmd),
                       .errorhandling = 'pass',
                       .options.snow  = opts) %dopar% {
    safe_runEwE(cmdfile=files.cmd[i], do.obj = obj.fxn, fit.abs.catch = fit.abs.catch)
  }
  close(pbar)
  runLL <- unlist(runLL.out)
  
  ##missing runs----------------------------------------------------------------
  # Retry generates a FRESH random parvec (via fn.GApop) for each failed run —
  # rationale: failures usually indicate a bad parameter combination, so a new
  # draw is more likely to succeed than a rerun of the same parvec. The new
  # cmd file replaces the failed one in files.cmd, and the substituted parvec
  # is recorded so the caller can update its population matrix accordingly.
  substituted <- list()   # names = character(index); values = numeric parvecs
  erruns <- which(is.na(runLL))
  # Delete any dirs from the initial (first-attempt) failures before retrying.
  if(length(erruns) > 0L) unlink(dirname(files.cmd[erruns]), recursive = TRUE)

  tries  <- 0
  while(length(erruns) >= 1L && tries < max_tries){
    tries <- tries + 1L
    n_err <- length(erruns)

    # Fresh random population; only the rows at `erruns` indices are used.
    gapop.err <- fn.GApop(pardist=pardist, vuldist=vuldist, rtdist=rtdist, fltdist=fltdist)
    # Retry cmd files land in the CURRENT gen's folder so they live alongside
    # their siblings and get cleaned up together. If `gen` isn't numeric (e.g.
    # base run label -1 or a string), fall back to 0.
    retry_g <- if(is.numeric(gen) && length(gen) == 1L && is.finite(gen)) as.integer(gen) else 0L
    new_cmd <- vapply(erruns, function(i)
      fn.parvec2cmd(par_vec = gapop.err[i, ], g = retry_g, idx = i),
      character(1L))

    cat(sprintf("\n[%s retry %d] %d failed runs (fresh random parvec)\n",
                gen_label, tries, n_err))
    pbar.err <- utils::winProgressBar(
      title = sprintf("Ecospace GA - %s retry %d  (%d runs)", gen_label, tries, n_err),
      label = sprintf("0 / %d completed (0%%)", n_err),
      min   = 0, max = n_err, initial = 0, width = 420)
    on.exit(try(close(pbar.err), silent = TRUE), add = TRUE)
    prog.err <- function(i){
      utils::setWinProgressBar(pbar.err, i,
        label = sprintf("%d / %d completed (%.1f%%)", i, n_err, 100 * i / n_err))
    }
    opts.err <- list(progress = prog.err)

    runLL.err.out <- foreach(k = seq_along(erruns),
                             .errorhandling = 'pass',
                             .options.snow  = opts.err) %dopar% {
      safe_runEwE(cmdfile = new_cmd[k], do.obj = obj.fxn, fit.abs.catch = fit.abs.catch)
    }
    close(pbar.err)
    runLL.err <- unlist(runLL.err.out)

    # Update state: swap in the retry cmd path for every slot we tried this
    # round (so a subsequent retry deletes the CURRENT attempt's dir, not the
    # original). For slots that succeeded, record the substituted parvec.
    for(k in seq_along(erruns)){
      idx <- erruns[k]
      runLL[idx]      <- runLL.err[k]
      files.cmd[idx]  <- new_cmd[k]
      if(!is.na(runLL.err[k]))
        substituted[[as.character(idx)]] <- gapop.err[idx, ]
    }

    still_failing <- which(is.na(runLL[erruns]))
    if(length(still_failing) > 0L)
      unlink(dirname(files.cmd[erruns[still_failing]]), recursive = TRUE)
    erruns <- which(is.na(runLL))
  }

  # Optional: warn if some runs still missing after retries
  if (length(erruns) > 0) {
    warning(sprintf("Incomplete: %d runs still missing after %d retries at indices: %s",
                    length(erruns), tries, paste(erruns, collapse = ", ")))
  }
  
  # Delete only the directories this call created (dirname of each cmd file),
  # not everything under run_dir. Prevents blowing away sibling gen<N>/ folders
  # or elite dirs that were moved forward for the next generation.
  if (delete.output && length(files.cmd) > 0L) {
    unlink(dirname(files.cmd), recursive = TRUE)
  }
  
  # Return fitness with two sync attributes so the caller (fn.GA) can keep
  # its parallel arrays consistent after any retry substitutions:
  #   cmd.files     -- updated paths (retry cmd files replace failed originals)
  #   substituted   -- named list; names are character(index), values are the
  #                    fresh random parvecs that replaced failed originals.
  #                    Empty list if no substitutions occurred.
  fitness <- runLL
  attr(fitness, "cmd.files")   <- files.cmd
  attr(fitness, "substituted") <- substituted
  # print(paste('Run time', round(as.numeric(Sys.time() - t1), 2)))
  return(fitness)
} # eof


# === Selection ===
#' @title Select parents by rank-based sampling (minimization)
#' @description
#' Ranks individuals by fitness for a MINIMIZATION objective (lower = better),
#' converts ranks to selection probabilities, and samples a new parent population
#' with replacement. The lowest-fitness individual receives the highest rank and
#' therefore the highest probability of being chosen.
#'
#' @param gapop Matrix or data frame of individuals (rows = individuals, columns = parameters).
#' @param fitness Numeric vector of fitness values (length equals nrow(gapop)); lower is better.
#'
#' @details
#' Ranks are computed via \code{rank(-fitness, ties.method='random')} so the smallest
#' fitness gets rank \code{nrow(gapop)} (i.e., the highest selection weight). Selection
#' pressure is roughly linear: the best individual is only ~2x more likely than the
#' median to be picked. Use \code{\link{tournament_select}} for stronger pressure.
#'
#' @return Matrix of selected parents with the same dimensions as \code{gapop}.
#' @export
select_parents <- function(gapop, fitness) {
  ranks <- rank(-fitness, ties.method='random')
  probs <- ranks / sum(ranks)
  selected <- gapop[sample(1:nrow(gapop), nrow(gapop), replace = TRUE, prob = probs), ]
  return(selected)
}

#' @title Select parents by tournament (minimization)
#' @description
#' Tournament selection for a MINIMIZATION objective. For each parent slot, draw
#' \code{k} individuals at random (without replacement) and keep the one with the
#' lowest fitness. Stronger and more tunable selection pressure than rank-based
#' sampling: larger \code{k} pushes harder toward the current best, \code{k=2} is
#' near-rank-equivalent.
#'
#' @param gapop Matrix or data frame of individuals (rows = individuals, columns = parameters).
#' @param fitness Numeric vector of fitness values (lower is better). NA fitnesses
#'   are treated as worst (lose every tournament unless all contenders are NA, in
#'   which case a random contender wins).
#' @param k Tournament size. Default 3.
#'
#' @return Matrix of selected parents with the same dimensions as \code{gapop}.
#' @export
tournament_select <- function(gapop, fitness, k = 3){
  n <- nrow(gapop)
  if(k < 1L) stop("tournament_select: k must be >= 1")
  if(k > n)  k <- n
  parents <- gapop[integer(0L), , drop = FALSE]
  parents <- gapop[rep(1L, n), , drop = FALSE]   # preserves colnames/types
  for(i in seq_len(n)){
    contenders <- sample.int(n, k, replace = FALSE)
    f <- fitness[contenders]
    if(all(is.na(f))){
      winner <- contenders[sample.int(k, 1L)]
    } else {
      winner <- contenders[which.min(f)]          # ignores NAs
    }
    parents[i, ] <- gapop[winner, ]
  }
  parents
}

# === Crossover ===
#' @title One-point crossover for paired parents
#' @description
#' Applies one-point crossover to pairs of parents to create offspring.
#'
#' @param parents Matrix of parent individuals (rows = individuals, columns = parameters).
#'
#' @details
#' Processes rows in pairs (1–2, 3–4, …). With probability 0.8, a single crossover
#' point is drawn uniformly from positions 1 to n_pars - 1, and the trailing gene
#' segments are swapped between the pair.
#'
#' Requires global variables `pop_size` (number of rows in `parents`) and `n_pars`
#' (number of columns) to be defined in the calling environment.
#'
#' @return Matrix of offspring with the same dimensions as `parents`.
#' @export
crossover <- function(parents) {
  offspring <- parents
  for (i in seq(1, pop_size - 1, by = 2)) {
    if (runif(1) < 0.8) {
      point <- sample(1:(n_pars - 1), 1)
      temp <- offspring[i, (point + 1):n_pars]
      offspring[i, (point + 1):n_pars] <- offspring[i + 1, (point + 1):n_pars]
      offspring[i + 1, (point + 1):n_pars] <- temp
    }
  }
  return(offspring)
}


#' @title Uniform crossover with grouped (linked) parameters
#' @description
#' Performs uniform crossover on pairs of parents, but treats specified columns
#' as linked "groups" that must be swapped together.
#'
#' @param parents Numeric matrix. Rows = individuals, columns = parameters.
#' @param group_id Vector of length ncol(parents) giving the group ID of each column.
#'   Columns with the same ID form a linkage block and will be swapped together.
#'   Use a unique ID per column to leave it ungrouped. NA IDs are treated as unique.
#' @param p_cross Probability of attempting crossover for a parent pair.
#' @param p_group Probability a group is swapped when crossover occurs.
#'   Either a scalar (applied to all groups) or a vector of length = number of groups.
#'
#' @details
#' Pairs rows (1&2, 3&4, ...) and for each pair attempts crossover with probability
#' `p_cross`. If crossover happens, a binary mask is drawn over groups with probability
#' `p_group`, and all columns in the selected groups are swapped.
#'
#' @return Matrix of offspring with the same dimensions as `parents`.
#' @export
crossover_uniform <- function(parents, group_id=1:ncol(parents), p_cross = 0.8, p_group = 0.5) {
  if (!is.matrix(parents)) stop("`parents` must be a matrix.")
  n_pars <- ncol(parents)
  if (length(group_id) != n_pars)
    stop("`group_id` must have length equal to ncol(parents).")
  
  # Treat NA group IDs as unique (i.e., their own single-gene group)
  if (anyNA(group_id)) {
    na_idx <- which(is.na(group_id))
    group_id[na_idx] <- paste0("__NA__", na_idx)
  }
  
  # Build list of column indices per group; drop ensures no empty groups
  groups <- split(seq_len(n_pars), as.factor(group_id), drop = TRUE)
  n_groups <- length(groups)
  
  # Handle scalar or group-wise p_group
  if (length(p_group) == 1L) {
    p_group_vec <- rep(p_group, n_groups)
  } else if (length(p_group) == n_groups) {
    p_group_vec <- as.numeric(p_group)
  } else {
    stop("`p_group` must be length 1 or equal to the number of groups (", n_groups, ").")
  }
  
  offspring <- parents
  pop_size <- nrow(parents)
  
  for (i in seq(1L, pop_size - 1L, by = 2L)) {
    if (runif(1L) < p_cross) {
      swap_groups <- runif(n_groups) < p_group_vec
      if (any(swap_groups)) {
        idx <- unlist(groups[swap_groups], use.names = FALSE)
        tmp <- offspring[i, idx, drop = FALSE]
        offspring[i, idx] <- offspring[i + 1L, idx]
        offspring[i + 1L, idx] <- tmp
      }
    }
  }
  offspring
}


# === Mutation ===
#' @title Mutate population with parameter-wise random resets
#' @description
#' For each gene of each individual, with probability \code{mutation_rate} redraw
#' the value within parameter-specific bounds. The vulnerability bounds adapt to
#' the current population mean (so the search drifts up or down with the population)
#' and vul and red tide draws are log-uniform to match the initial-population
#' sampling philosophy. All other parameters use fixed \code{L.bounds}/\code{U.bounds}
#' with linear-uniform draws.
#'
#' @param population Matrix of individuals (rows = individuals, columns = parameters).
#'
#' @details
#' Adaptive vul bounds: \code{Gamma(shape=2, mean=pop_mean)}; bounds are the 1\% and
#' 99\% quantiles, clipped to \code{[1.01, 1e6]}. So as selection pulls the population
#' mean upward (or downward), so does the mutation search window. Vul draws inside
#' those bounds are log-uniform.
#'
#' Red tide draws are log-uniform over the fixed \code{[L.bounds, U.bounds]} for
#' \code{redtide.par.idx} (symmetric around the baseline of 1 in log space, matching
#' \code{fn.GApop(rtdist='log-uniform')}).
#'
#' Requires global variables: \code{n_pars}, \code{mutation_rate}, \code{vul.par.idx},
#' \code{redtide.par.idx}, \code{par.labels}, \code{L.bounds}, \code{U.bounds}.
#'
#' @return Mutated population matrix (same dimensions as \code{population}).
#' @export
mutate <- function(population) {
  low <- L.bounds
  upp <- U.bounds

  # --- Adaptive vul bounds (drift with population mean) ----------------------
  # Gamma(shape=2, rate=2/avg) has mean=avg and a right-skewed shape with a wide
  # 1%-99% range (~avg/40 to ~avg*3.3). Clipped to [1.01, 1e6] for safety.
  if(length(vul.par.idx) > 0L){
    avg <- colMeans(population)
    low[vul.par.idx] <- qgamma(0.01, shape = 2, rate = 2/avg[vul.par.idx])
    upp[vul.par.idx] <- qgamma(0.99, shape = 2, rate = 2/avg[vul.par.idx])
    low[vul.par.idx] <- pmax(low[vul.par.idx], 1.01)
    upp[vul.par.idx] <- pmin(upp[vul.par.idx], 1e6)
  }

  # --- Per-individual mutation -----------------------------------------------
  # Default: linear-uniform within [low, upp]. For vuls and red tide, override
  # to log-uniform so the draws aren't biased toward the high end of bounds that
  # span multiple orders of magnitude.
  for (i in seq_len(nrow(population))) {
    mask <- runif(n_pars) < mutation_rate
    if(!any(mask)) next

    # baseline linear-uniform draws for all masked genes
    draws <- runif(sum(mask), low[mask], upp[mask])

    # locate which masked genes are vul / RT (positions within the masked subset)
    mask_idx <- which(mask)
    is_vul_in_mask <- mask_idx %in% vul.par.idx
    is_rt_in_mask  <- mask_idx %in% redtide.par.idx

    if(any(is_vul_in_mask)){
      vul_genes <- mask_idx[is_vul_in_mask]
      lo <- pmax(low[vul_genes], 1.01)
      hi <- pmax(upp[vul_genes], lo * 1.0001)
      draws[is_vul_in_mask] <- exp(runif(length(vul_genes), log(lo), log(hi)))
    }
    if(any(is_rt_in_mask)){
      rt_genes <- mask_idx[is_rt_in_mask]
      lo <- pmax(low[rt_genes], 1e-6)
      hi <- pmax(upp[rt_genes], lo * 1.0001)
      draws[is_rt_in_mask] <- exp(runif(length(rt_genes), log(lo), log(hi)))
    }

    population[i, mask] <- round(draws, 4)
  }

  integer.idx <- grep("xbase", par.labels)
  if(length(integer.idx) > 0L) population[, integer.idx] <- round(population[, integer.idx])
  return(population)
} #eof

#' @keywords internal
#' @noRd
cluster_is_ok <- function() {
  ok <- tryCatch({
    res <- foreach(i = 1:2, .combine = c) %dopar% { Sys.getpid() }
    length(res) == 2
  }, error = function(e) FALSE)
  ok
}

#' @title Canonical list of objects to clusterExport for R4EwE GA/sensitivity workers.
#' @description Returns the character vector of object names that every worker process
#'   needs in its global environment to call \code{safe_runEwE} → \code{fn.runEwE} →
#'   any of the supported objective functions (\code{fn.objfxn1}, \code{fn.objfxn2},
#'   \code{fn.ecospace_objfxn}). Use it in your top-level \code{clusterExport()} call
#'   so the three clusterExport sites (script, \code{ensure_cluster}, \code{fn.GA}'s
#'   every-5-gen reset) stay in sync.
#' @return Character vector of object names.
#' @export
.r4ewe_worker_exports <- function(){
  required <- c(
    # console + objective functions
    "file.console",
    "safe_runEwE", "fn.runEwE", ".check_ecospace_output",
    "fn.objfxn1", "fn.objfxn2", "fn.ecospace_objfxn",
    # per-region prediction readers and their internal helpers
    "fn.read_pred_ecospace_wide", "fn.read_pred_ecospace_discards_split", "fn.fg_meta",
    ".find_year_skip", ".detect_ecospace_regions", ".ecospace_dmort_by_year",
    # legacy single-region readers (still used by fn.objfxn1)
    "fn.ecospace_predB_ts2array", "fn.ecospace_predC_ts2array",
    # spatial LL package helpers (phase 2 task 7c) - always needed for the
    # worker to be able to run the spatial component
    "fn.spatial_LL", "fn.ecospace_spatial_pred",
    ".parse_ecospace_map_csv", ".read_asc", ".normalise_map",
    # data the objective functions need
    "obs.ts", "group.names", "fleet.names", "df.names",
    # model years and run context
    "styear", "enyear",
    "cmd_base"
  )
  # Optional session objects: only export if they exist in .GlobalEnv (or
  # search path). Lets a phase-1 session (no spatial fit configured) or a
  # sensitivity-only session (no GA myconfig / run_dir) run without needing
  # to define these.
  optional <- c("spatial.obs", "spatial.weight", "model_styear",
                "myconfig", "run_dir", "expected_regions")
  optional <- optional[vapply(optional,
                              function(nm) exists(nm, envir = .GlobalEnv),
                              logical(1))]
  c(required, optional)
}

#' @keywords internal
#' @noRd
ensure_cluster <- function(workers = NULL) {
  if (!cluster_is_ok()) {
    if(exists("cl", envir = .GlobalEnv)){
      try(stopCluster(get("cl", envir = .GlobalEnv)), silent = TRUE)
      rm(list = "cl", envir = .GlobalEnv)
    }
    gc()
    if(is.null(workers)) workers <- max(1L, floor(detectCores() / 4))
    cl <<- parallel::makePSOCKcluster(workers, outfile = "cluster_workers.log")
    doSNOW::registerDoSNOW(cl)
    clusterExport(cl, .r4ewe_worker_exports())
  }
}#eof

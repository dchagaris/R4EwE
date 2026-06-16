
#prepare parameter vector----
#' @title Construct parameter vectors and bounds for model estimation
#' @description Builds the parameter vector and related metadata (labels, groups, bounds, CVs) for 
#' vulnerabilities, environmental responses, dispersal, and mediation parameters. Several objects 
#' are assigned in the parent/global environment.
#' @param do.vuls Logical; include vulnerability parameters.
#' @param vul_pars Data frame of vulnerability parameters.
#' @param vul.min Minimum vulnerability (unused placeholder).
#' @param vul.max Maximum vulnerability (unused placeholder).
#' @param vul.cv Coefficient of variation for vulnerabilities.
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
#' @return Invisibly returns \code{NULL}. Several objects are set globally.
#' @export
fn.makeparvec <- function(
  do.vuls = TRUE, 
  vul_pars = predprey_pars,
  vul.min = 1.01,
  vul.max = 1e6,
  vul.cv = 0.4,
  
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

  do.redtide = FALSE,
  redtide_pars = NULL,
  redtide.min = 0.5,
  redtide.max = 1.5,
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
  ##vulnerabilities - do these on the log scale
  lower.vuls <- upper.vuls <- numeric()
  if(do.vuls){
  lower.vuls <- exp(log(vul_pars$base.val) + exp(vul.cv)^2 * qnorm(0.05))
  lower.vuls <- ifelse(lower.vuls<=1.0,1.01,lower.vuls)
  upper.vuls <- exp(log(vul_pars$base.val) + exp(vul.cv)^2 * qnorm(0.95))
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
  
  ##red tide — uniform [redtide.min, redtide.max] multiplier bounds
  lower.redtide <- upper.redtide <- numeric()
  if(do.redtide && length(redtide_vec)>0){
    lower.redtide <- rep(redtide.min, length(redtide_vec))
    upper.redtide <- rep(redtide.max, length(redtide_vec))
  }

  ##dispersal
  lower.disp <- upper.disp <- numeric()
  if(do.disp){
  lower.disp <- disp_vec-2*disp.cv*disp_vec
  upper.disp <- disp_vec+2*disp.cv*disp_vec
  }
  
  ##fleet dynamics
  lower.flt <- upper.flt <- numeric()
  if(do.fltdyn){
  lower.flt <- rep(fltdyn.min,length(fleet_vec))
  upper.flt <- rep(fltdyn.max,length(fleet_vec))
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
#' @param vuldist Character string; distribution type for vulnerability parameters, either"uniform" (uniform between bounds) or 
#' "gamma" withc shape=2 and rate=shape/par. Default "uniform"
#' @return A numeric matrix of dimension run_config$popSize x n_pars.
#' @export
fn.GApop = function(run_config=myconfig, pardist='uniform', vuldist='uniform'){
  # run_config <- myconfig
  # pardist='uniform'
  # vuldist='gamma'
  #two things: fix vul so it goes below 2; fix mediation base_x values
  mat <- matrix(NA,nrow=run_config$popSize, ncol=n_pars)
  colnames(mat) <- par.labels
  
  #inspect the parameter set and their bounds
  data.frame(L.bounds, est_par_vec, U.bounds, par_cv_vec,par.groups)
  #get standard deviation from cv
  #par_sd = sqrt(par_cv_vec+1)

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
  pop.vuls[pop.vuls<1.01] <- 1.01

  mat[] <- pop.pars  
  mat[,vul.par.idx] <- pop.vuls
  mat <- round(mat,4)
  
  integer.idx <- grep("xbase",par.labels)
  if(length(integer.idx)>0) mat[,integer.idx] <- round(mat[,integer.idx])
  mat[,flt.par.idx] <- ifelse(mat[,flt.par.idx]<0,0,mat[,flt.par.idx])
  mat[,env.par.idx] <- ifelse(mat[,env.par.idx]<0,0,mat[,env.par.idx])
    
  return(mat)
  #matrix(rep(log_par_vec, run_config$popSize), nrow=run_config$popSize, byrow=T)
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
#' @return Path to the written command file.
#' @export
fn.parvec2cmd <- function(par_vec=est_par_vec, g=0, idx=0, out_dir=run_dir){
  
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
#' @param highlight_pars Numeric vector of length \code{ncol(gapop)} whose values get the
#'   blue highlight (vline on histograms, overlaid curve on env/rt panels). If
#'   \code{NULL} (default), \code{est_par_vec} is used — i.e. the baseline is highlighted,
#'   matching the original behavior used for the initial-population plot. Pass the
#'   min-LL individual (e.g. \code{gapop[which.min(fitness), ]}) for the final-population
#'   posterior plot.
#' @return Invisibly returns \code{file}.
#' @export
fn.plot_ga_pop_responses <- function(gapop, est_par_vec, par.groups, par.labels,
                                     env_pars = NULL, redtide_pars = NULL,
                                     file, dims = c(3, 3),
                                     highlight_pars = NULL){
  if(!is.null(env_pars))     env_pars     <- .normalize_env_pars(env_pars)
  if(!is.null(redtide_pars)) redtide_pars <- .normalize_env_pars(redtide_pars)

  # When no highlight set is provided, fall back to the baseline est_par_vec so existing
  # callers (initial-pop plot) keep their current behavior. For env/rt curves a NULL
  # highlight means the blue curve is the BASE curve (adj1 = adj2 = 1), which matches
  # the original "baseline curve highlighted in blue" behavior.
  hl <- if(is.null(highlight_pars)) est_par_vec else highlight_pars

  alpha_col <- grDevices::rgb(0, 0, 0, alpha = 0.15)

  grDevices::pdf(file, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = dims, mar = c(3, 3, 2, 1), mgp = c(1.6, 0.5, 0))

  draw_hists <- function(col_idx){
    for(j in col_idx){
      graphics::hist(gapop[, j], breaks = 50, main = par.labels[j], xlab = 'value')
      graphics::abline(v = hl[j], col = 'blue', lty = 2)
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

    # Highlight curve adj values: 1,1 when no override (= base curve), else from highlight_pars.
    h_adj1 <- if(is.null(highlight_pars)) 1 else highlight_pars[col_idx[1]]
    h_adj2 <- if(is.null(highlight_pars) || length(col_idx) < 2) 1 else highlight_pars[col_idx[2]]

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
      graphics::lines(x, ycurve(Inf.base * h_adj1, Slope.base * h_adj2), col = 'blue', lwd = 2)

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
      hp <- trap_pts(h_adj1, h_adj2)
      graphics::lines(hp$x, hp$y, col = 'blue', lwd = 2)

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
      graphics::lines(x, norm_y(Mean.base * h_adj1,
                                SDLeft.base * h_adj2,
                                SDRight.base * h_adj2,
                                Max.base),
                      col = 'blue', lwd = 2)

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
      graphics::lines(x, ycurve(XMid.base * h_adj1, Slope.base * h_adj2), col = 'blue', lwd = 2)

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
    vuldist = if(exists("gapop.vuldist", envir=.GlobalEnv)) get("gapop.vuldist", envir=.GlobalEnv) else "uniform"
){
  t1 <- Sys.time()

  # console progress bar; requires doSNOW backend (see ensure_cluster)
  gen_label <- if(is.null(gen)) "GA pop" else sprintf("Gen %d", gen)
  cat(sprintf("\n[%s] %d simulations\n", gen_label, length(files.cmd)))
  pbar <- utils::txtProgressBar(min = 0, max = length(files.cmd), style = 3)
  prog <- function(n) utils::setTxtProgressBar(pbar, n)
  opts <- list(progress = prog)

  runLL.out <- foreach(i = seq_along(files.cmd),
                       .errorhandling = 'pass',
                       .options.snow  = opts) %dopar% {
    safe_runEwE(cmdfile=files.cmd[i], do.obj = obj.fxn, fit.abs.catch = fit.abs.catch)
  }
  close(pbar)
  runLL <- unlist(runLL.out)
  
  ##missing runs----------------------------------------------------------------
  #tmp.errs = sample(1:340,10)
  #runLL[11:15] <- NA_real_ 
  erruns <- which(is.na(runLL))
  #filecheck <- sapply(dirname(files.cmd),FUN=function(x)length(list.files(x)))
  unlink(dirname(files.cmd[erruns]), recursive=T)
  #erruns <- which(filecheck<=1)  
  tries <- 0
  while(length(erruns)>=1 && tries<max_tries){
    tries <- tries+1

    gapop.err <- fn.GApop(pardist=pardist, vuldist=vuldist)
    files.cmd.err <- lapply(erruns,function(i) fn.parvec2cmd(par_vec=gapop.err[i,], g=0, idx=i))
    files.cmd.err <- unlist(files.cmd.err, use.names=F)

    cat(sprintf("\n[%s retry %d] %d failed runs\n", gen_label, tries, length(erruns)))
    pbar.err <- utils::txtProgressBar(min = 0, max = length(erruns), style = 3)
    prog.err <- function(n) utils::setTxtProgressBar(pbar.err, n)
    opts.err <- list(progress = prog.err)

    runLL.err.out <- foreach(i = seq_along(erruns),
                             .errorhandling = 'pass',
                             .options.snow  = opts.err) %dopar% {
      safe_runEwE(cmdfile=files.cmd.err[i], do.obj = obj.fxn, fit.abs.catch = fit.abs.catch)
    }
    close(pbar.err)
    runLL.err <- unlist(runLL.err.out)
    
    #for(k in 1:length(erruns)) runs[[erruns[k]]] <- unlist(runs.erruns[k])
    for (k in seq_along(erruns)) runLL[erruns[k]] <- runLL.err[k]

    # Preserve the index map: point files.cmd at the new retry cmd files for the
    # indices we just re-ran, so the next round's unlink targets the correct run dirs.
    # Do NOT re-list run_dir here — list.files() returns alphabetical order, which
    # breaks the runLL <-> files.cmd index alignment.
    for (k in seq_along(erruns)) files.cmd[erruns[k]] <- files.cmd.err[k]
    erruns <- which(is.na(runLL))
    if(length(erruns)>0) unlink(dirname(files.cmd[erruns]), recursive=T)
  }
  
  
  # Optional: warn if some runs still missing after retries
  if (length(erruns) > 0) {
    warning(sprintf("Incomplete: %d runs still missing after %d retries at indices: %s",
                    length(erruns), tries, paste(erruns, collapse = ", ")))
  }
  
  # Clean up any generated output directories, if applicable
  if (delete.output && exists("run_dir")) {
    unlink(list.dirs(run_dir, full.names = TRUE, recursive = FALSE), recursive = TRUE)
  }
  
  # Return fitness vector (total log-likelihoods to minimize)
  fitness <- runLL
  # print(paste('Run time', round(as.numeric(Sys.time() - t1), 2)))
  return(fitness)
} # eof


# === Selection ===
#' @title Select parents by rank-based sampling
#' @description
#' Ranks individuals by fitness (higher is better), converts ranks to selection
#' probabilities, and samples a new parent population with replacement.
#'
#' @param gapop Matrix or data frame of individuals (rows = individuals, columns = parameters).
#' @param fitness Numeric vector of fitness values (length equals nrow(gapop)); higher is better.
#'
#' @details
#' Ranks are computed with ties broken at random, then normalized to probabilities.
#' Parents are sampled with replacement to keep population size constant.
#'
#' @return Matrix of selected parents with the same dimensions as `gapop`.
#' @export
select_parents <- function(gapop, fitness) {
  ranks <- rank(-fitness, ties.method='random')
  probs <- ranks / sum(ranks)
  selected <- gapop[sample(1:nrow(gapop), nrow(gapop), replace = TRUE, prob = probs), ]
  return(selected)
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
#' Mutates individuals by redrawing selected genes from a range informed by the
#' current population, with a configurable margin.
#'
#' @param population Matrix of individuals on the log scale (rows = individuals, columns = 
#' parameters).
#' @param margin Fractional expansion applied to the min–max range per parameter (on orig scale).
#'
#' @details
#' The algorithm:
#' 1) Convert each parameter to original scale (exp), compute min and max across the population.  
#' 2) Expand bounds by `margin` (lower reduced, upper increased); enforce vulnerability
#'    parameters `vul.par.idx` to be strictly > 1.  
#' 3) Transform bounds back to log scale and, for each individual, redraw genes indicated
#'    by a Bernoulli mask (`mutation_rate`) uniformly within the parameter-specific bounds.
#'
#' Requires global variables: `n_pars`, `mutation_rate`, and `vul.par.idx`.
#'
#' @return Mutated population matrix (same dimensions as `population`).
#' @export
mutate <- function(population, margin=0.2) {
  low = L.bounds
  upp = U.bounds
  #for vuls allow to creep up/down with gamma draws
  avg = apply(population,2,mean)
  low[vul.par.idx] <- qgamma(0.01, shape = 2, rate = 2/avg[vul.par.idx])
  upp[vul.par.idx] <- qgamma(0.99, shape = 2, rate = 2/avg[vul.par.idx])
  
  low[vul.par.idx] <- ifelse(low[vul.par.idx]<=1.0,1.01,low[vul.par.idx])
  upp[vul.par.idx] <- ifelse(upp[vul.par.idx]>1e6,1e6,upp[vul.par.idx])
  #back to log scale
  #low = log(low)
  #upp = log(upp)
  for (i in 1:nrow(population)) {
        mask <- runif(n_pars) < mutation_rate #mask are the parameters to mutate
        population[i, mask] <- round(runif(sum(mask), low[mask], upp[mask]),4)
  }
  integer.idx <- grep("xbase",par.labels)
  population[,integer.idx] <- round(population[,integer.idx])
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
  c(
    # console + objective functions
    "file.console",
    "safe_runEwE", "fn.runEwE",
    "fn.objfxn1", "fn.objfxn2", "fn.ecospace_objfxn",
    # per-region prediction readers and their internal helpers
    "fn.read_pred_ecospace_wide", "fn.read_pred_ecospace_discards_split", "fn.fg_meta",
    ".find_year_skip", ".detect_ecospace_regions", ".ecospace_dmort_by_year",
    # legacy single-region readers (still used by fn.objfxn1)
    "fn.ecospace_predB_ts2array", "fn.ecospace_predC_ts2array",
    # data the objective functions need
    "obs.ts", "group.names", "fleet.names", "df.names",
    # model years and run context
    "styear", "enyear",
    "cmd_base", "myconfig", "run_dir"
  )
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

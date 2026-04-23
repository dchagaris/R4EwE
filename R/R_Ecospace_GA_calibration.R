
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
  vul.cv = 0.3,
  
  do.env = FALSE, 
  env_pars = NULL,
  env.min=0.25,
  env.max=2,
  env.cv=0.2,
  
  do.disp = TRUE,
  disp_pars = disp_pars,
  disp.cv=0.1,
  disp.min=3,
  disp.max=3000,
  
  do.med = FALSE,
  med_pars = NULL,
  med.xbase.cv = .1,
  
  do.fltdyn = FALSE,
  fltdyn_pars=NULL,
  fltdyn.min=0,
  fltdyn.max=5,
  fltdyn.cv=0.2){
  
  # do.vuls=TRUE; vul_pars=predprey_pars; vul.min=1.01; vul.max=1e6; vul.cv=0.4;
  # do.env=TRUE; env_pars=env_pars; env.min=0.1; env.max=2; env.cv=0.2;
  # do.disp=TRUE; disp_pars=disp_pars; disp.cv=0.2;
  # do.med = FALSE; med_pars = NULL; med.xbase.cv = .1;
  # do.fltdyn = TRUE; fltdyn_pars=fltdyn_pars; fltdyn.min=0; fltdyn.max=5; fltdyn.cv=0.2
  
  #validate setup----
  if(do.vuls){
    pdcol = unique(vul_pars$pred[is.na(vul_pars$prey)])
    pdpy = c(0,unique(vul_pars$pred[!is.na(vul_pars$prey)]))
    if(length(which(pdcol %in% pdpy))>0){
      stop("Trying to estimate ki and kij for at least on predator. Check vul_pars input.")
    }
  }
  
  #build parameter vector----
  n_vuls <- n_env <- n_disp <- n_med <- n_fleet <- 0
  if(do.vuls) n_vuls <- nrow(vul_pars)
  if(do.env) n_env <- nrow(env_pars)  #need to melt the env resp fxn dataframe
  if(do.disp) n_disp <- nrow(disp_pars)
  if(do.med) n_med <- nrow(med_pars)
  if(do.fltdyn) n_fleet <- nrow(fltdyn_pars)
  #medations functions
  
  n_pars <- n_vuls + n_env + n_disp + n_med + n_fleet
  if (n_pars == 0) stop("n_pars==0: No predator-prey pairs or environmental responses parameters to estimate.  Check inputs.")
  
  #no longer taking log of parameters
  vuln_vec <- env_vec <- disp_vec <- med_vec <- fleet_vec <- numeric()
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
    for(i in 1:n_env){
      #i=1
      #TRAPEZOID SHAPE
      if(env_pars$shape.idx[i]==9){
        #par_vec.i <- c(env_pars$par1[i],env_pars$par2[i],env_pars$par3[i],env_pars$par4[i])  
        par_vec.i <- rep(1,2)
        par.labels.i <- paste(paste0('env',rep(env_pars$resp.idx[i],length(par_vec.i))),paste0('type',rep(env_pars$shape.idx[i],length(par_vec.i))),c('mid.adj','width.adj'),sep="_")
        par.groups.i <- paste0('env',rep(env_pars$resp.idx[i],length(par_vec.i)))
      } #end trapezoid
      env_vec <- c(env_vec,par_vec.i)
      par.labels <- c(par.labels,par.labels.i)
      par.groups <- c(par.groups,par.groups.i)
    } #end loop
    par.idx <- c(par.idx,rep('env',length(env_vec)))
  } #end do.env
  
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
    par.groups = c(par.groups,paste('flt',fltdyn_pars$fleet.idx,sep="_"))
    #log.par.idx <- c(log.par.idx,rep(TRUE,length(log_disp_vec)))
  }
  
  
  
  #n_pars <<- length(log_vuln_vec) + length(env_vec) + length(log_disp_vec) + length(med_vec)
  n_pars <<- length(vuln_vec) + length(env_vec) + length(disp_vec) + length(med_vec) + length(fleet_vec)
  if (n_pars == 0) stop("n_pars==0: No parameters to estimate.  Check inputs.")
  
  #est_par_vec <- c(log_vuln_vec, env_vec, log_disp_vec, med_vec)
  est_par_vec <- c(vuln_vec, env_vec, disp_vec, med_vec, fleet_vec)
  names(est_par_vec) <- par.labels
  
  # index parameter types
  vul.par.idx <<- which(par.idx=='vul')
  env.par.idx <<- which(par.idx=='env')
  disp.par.idx <<- which(par.idx=='disp')
  med.par.idx <<- which(par.idx=='med')
  flt.par.idx <<- which(par.idx=='fleetdyn')
  #log.par.idx <<- log.par.idx
  
  # super assignment so these are available
  vul_pars <<- vul_pars
  env_pars <<- env_pars
  disp_pars <<- disp_pars
  med_pars <<- med_pars
  fltdyn_pars <<- fltdyn_pars
  
  # set parameter bounds--------------------------------------------------------
  ##vulnerabilities - do these on the log scale
  lower.vuls <- upper.vuls <- numeric()
  lower.vuls <- exp(log(vul_pars$base.val) + exp(vul.cv)^2 * qnorm(0.05))
  lower.vuls <- ifelse(lower.vuls<=1.0,1.01,lower.vuls)
  upper.vuls <- exp(log(vul_pars$base.val) + exp(vul.cv)^2 * qnorm(0.95))
  
  ##environmental shapes - TBD
  lower.env <- upper.env <- numeric()
  if(n_env>0){
    for(i in 1:length(env_vec)){
      #i=1
      label.i <- par.labels[env.par.idx[i]]
      shape.i <- as.numeric(gsub('env','',unlist(strsplit(label.i,"_"))[1]))
      type.i <- as.numeric(gsub('type','',unlist(strsplit(label.i,"_"))[2]))
      par.i <- unlist(strsplit(label.i,"_"))[3]
      val.i <- env_vec[i]
      if(type.i==9){
        lower.i <- env.min
        upper.i <- env.max
        
      } #end trapezoid
      lower.env <- c(lower.env,lower.i)
      upper.env <- c(upper.env,upper.i)
    } #end i loop
  } #end if env

  
  ##dispersal
  lower.disp <- upper.disp <- numeric()
  lower.disp <- disp_vec-2*disp.cv*disp_vec
  upper.disp <- disp_vec+2*disp.cv*disp_vec
  
  ##fleet dynamics
  lower.flt <- upper.flt <- numeric()
  lower.flt <- rep(fltdyn.min,length(fleet_vec))
  upper.flt <- rep(fltdyn.max,length(fleet_vec))
  
  ##mediation
  lower.med <- upper.med <- numeric()  
  if(n_med>0){
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
  L.bounds <<- round(c(lower.vuls,lower.env,lower.disp,lower.med,lower.flt),2)
  U.bounds <<- round(c(upper.vuls, upper.env, upper.disp,upper.med,upper.flt),2)
  n_vuls <<- n_vuls
  n_env <<- length(env_vec)
  n_disp <<- n_disp
  n_med <<- length(med_vec)
  n_flt <<- n_fleet
  par.groups <<- par.groups
  par.labels <<- par.labels
  
  #create vector of cv for rnorm() draws
  par_cv_vec <- est_par_vec
  par_cv_vec[vul.par.idx] <- vul.cv
  par_cv_vec[env.par.idx] <- env.cv
  par_cv_vec[disp.par.idx] <- disp.cv
  par_cv_vec[med.par.idx] <- med.cv
  par_cv_vec[flt.par.idx] <- fltdyn.cv
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
  tags.vul <- tags.env <- tags.disp <- tags.med <- tags.flt <- character()
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
    tags.env <- character()
    for(i in 1:length(env_groups)){
      #i=1
      par_set.i <- par_vec[which(par.groups==env_groups[i])]
      par_labels.i <- par.labels[which(par.groups==env_groups[i])]
      type.i <- as.numeric(gsub('type','',strsplit(par_labels.i[1],"_")[[1]][2]))
      shape.i <- as.numeric(gsub('env','',strsplit(par_labels.i[1],"_")[[1]][1]))
      #TRAPEZOID
      if(type.i==9){
        LB.base <- env_pars[env_pars$resp.idx==shape.i,'par1']
        LT.base <- env_pars[env_pars$resp.idx==shape.i,'par2']
        RT.base <- env_pars[env_pars$resp.idx==shape.i,'par3']
        RB.base <- env_pars[env_pars$resp.idx==shape.i,'par4']
        mid.base <- (RT.base+LT.base)/2
        pref.base <- RT.base-LT.base
        mid.i <- mid.base*par_set.i[1]
        pref.i <- pref.base*par_set.i[2]
        LT.i <- round(pmax(0,mid.i-pref.i/2),2)
        RT.i <- round(mid.i+pref.i/2,2)
        LB.i <- round(pmax(0,LT.i-(LT.base-LB.base)*par_set.i[2]),2)
        RB.i <- round(RT.i+(RB.base-RT.base)*par_set.i[2],2)
        par.string.i <- paste(type.i, LB.i, LT.i, RT.i, RB.i, sep=" ")
        tag.i <- paste0("<ECOSPACE_ENVIRONMENTAL_RESPONSE_INDEXED>(",shape.i,"), ",par.string.i,", Indexed.Single[]")
      }
      
      tags.env <- c(tags.env,tag.i)
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
  
  
  tags = c(tags.vul,tags.env,tags.disp, tags.med, tags.flt)
  
  
  #command files----
  cmd_j <- cmd_base
  cmd_j[startsWith(cmd_j, "<ECOSPACE_OUTPUT_DIR>")] <-
    sprintf("<ECOSPACE_OUTPUT_DIR>, %s, System.String, Updated", run_path)
  cmd_j <- c(cmd_j, tags)
  cmd_file <- file.path(run_path, "cmd.txt")
  writeLines(cmd_j, cmd_file)
  return(cmd_file)
  
}

#' @keywords internal
#' @noRd

safe_runEwE <- function(cmdfile, do.obj) {
  on.exit(gc(), add = TRUE)
  #out <- tryCatch({
    # example: quiet system call
    # system2("EwE.exe", args = c(cmdfile), stdout = FALSE, stderr = FALSE, wait = TRUE)
  out <- fn.runEwE(cmdfile = cmdfile, do.obj = do.obj)
  #}, error = function(e) NA_real_)
  if (is.numeric(out) && length(out) >= 1 && is.finite(out[1])) out[1] else NA_real_

# safe_runEwE <- function(cmdfile, do.obj, bug=F) {
#   if (!bug){
#     on.exit({
#       gc()
#       # Aggressively kill any lingering EwE console processes spawned by this worker
#       # This prevents zombies if the timeout is triggered
#       try(system("taskkill /F /IM EwE.exe /T", show.output.on.console = FALSE), silent = TRUE)
#       try(system("taskkill /F /IM Ecospace.exe /T", show.output.on.console = FALSE), silent = TRUE)
#     }, add = TRUE)
#     
#     out <- tryCatch({
#       # Wrap the original function in a strict timeout (e.g., 20 minutes / 1200 seconds)
#       if (requireNamespace("R.utils", quietly = TRUE)) {
#         R.utils::withTimeout({
#           fn.runEwE(cmdfile = cmdfile, do.obj = do.obj)
#         }, timeout = 1200, onTimeout = "error") # throw error on timeout to trigger NA_real_
#       } else {
#         # Fallback if R.utils isn't installed
#         fn.runEwE(cmdfile = cmdfile, do.obj = do.obj)
#       }
#     }, error = function(e) {
#       # If it times out or fails, return NA to let the GA know this run failed
#       return(NA_real_)
#     })
#     
#     if (is.numeric(out) && length(out) >= 1 && is.finite(out[1])) {
#       return(out[1]) 
#     } else {
#       return(NA_real_)
#     }
#   } else {
#     on.exit(gc(), add = TRUE)
#     run_dir <- dirname(cmdfile)
#     
#     # Ensure clean slate for output detection
#     target_file <- file.path(run_dir, "Ecospace_Annual_Average_Biomass.csv")
#     if(file.exists(target_file)) file.remove(target_file)
#     
#     # Normalize paths and handle spaces in "King Mackerel" via PowerShell triple-quotes
#     console_path <- normalizePath(file.console, winslash = "\\", mustWork = TRUE)
#     cmd_path     <- normalizePath(cmdfile, winslash = "\\", mustWork = TRUE)
#     
#     # PowerShell command construction
#     ps_cmd <- paste0("powershell -Command \"Start-Process -FilePath '", console_path, 
#                      "' -ArgumentList '\"\"", cmd_path, "\"\"' -PassThru | Select-Object -ExpandProperty Id\"")
#     
#     # Jitter to prevent concurrent Access DB collisions
#     Sys.sleep(runif(1, 0.1, 5.0)) 
#     
#     # Launch and capture PID, trimming whitespace/newlines to prevent taskkill failure
#     raw_pid <- tryCatch(suppressWarnings(system(ps_cmd, intern = TRUE)), error = function(e) NA)
#     pid <- if(is.character(raw_pid) && length(raw_pid) > 0) trimws(raw_pid[1]) else NA
#     
#     start_time <- Sys.time()
#     model_finished <- FALSE
#     
#     # Monitoring loop with 20-minute timeout
#     while(difftime(Sys.time(), start_time, units="mins") < 20) {
#       Sys.sleep(5) 
#       if(file.exists(target_file)) {
#         if(file.info(target_file)$size > 100) {
#           Sys.sleep(2) # Buffer for disk writing
#           model_finished <- TRUE
#           break
#         }
#       }
#     }
#     
#     # Force kill the specific PID and its child tree (/T) to release file locks
#     if(!is.na(pid) && pid != "") {
#       try(system(paste0("taskkill /T /F /PID ", pid), show.output.on.console = FALSE), silent = TRUE)
#     }
#     
#     if(model_finished) {
#       out <- tryCatch({
#         # Pass obs.ts explicitly from the worker's environment 
#         res <- fn.objfxn1(dir.pred = run_dir, obs.ts = obs.ts)
#         # Save results to local text file in the run directory
#         write.table(t(res), file = file.path(run_dir, "obj_results.txt"), 
#                     sep = ",", row.names = FALSE, col.names = TRUE)
#         
#         res # Return the full vector for the logic below
#       }, error = function(e) {
#         # Log the error to a file if objective function fails
#         writeLines(as.character(e), file.path(run_dir, "obj_error.txt"))
#         return(NA_real_)
#       })
#       
#       # Return the total likelihood (first element) to the GA 
#       if(is.numeric(out) && length(out) >= 1) return(out[1])
#     }
#     return(NA_real_)
#   }
# }#eof


# safe_runEwE <- function(cmdfile, do.obj, timeout_sec = 1) {
#   if (requireNamespace("R.utils", quietly = TRUE)) {
#     R.utils::withTimeout({
#       fn.runEwE(cmdfile = cmdfile, do.obj = do.obj)
#     }, timeout = timeout_sec, onTimeout = "silent")
#   } else {
#     # Fallback without timeout (consider installing R.utils)
#     fn.runEwE(cmdfile = cmdfile, do.obj = do.obj)
#   }
# }#eof


# safe_runEwE <- function(cmdfile, do.obj) {
#   on.exit(gc(), add = TRUE)   # ensure cleanup
#   out <- tryCatch(
#     fn.runEwE(cmdfile = cmdfile, do.obj = do.obj),
#     error = function(e) NA_real_)
#     if (is.numeric(out) && length(out) >= 1 && is.finite(out[1])){
#       out[1]
#     } else {
#       NA_real_
#     }
# }#eof


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
    delete.output=T,
    max_tries=5
){
  #source(file.setup)

  #clusterExport(cl,append(cl.export,c("file.console", "fn.runEwE", "safe_runEwE","fn.objfxn1", "fn.objfxn2","styear","enyear","group.names","df.names","fn.ecospace_predB_ts2array","fn.ecospace_predC_ts2array")))

  #pbar <- winProgressBar("Running Ecospace GApop",label=paste0("Simulation 0 of ",length(files.cmd)),max=100)
  #prog <- function(n) setWinProgressBar(pbar,(n/length(files.cmd)*100),label=paste("Simulation Run", n,"of",length(files.cmd),"Completed"))
  #opts <- list(progress=prog)
  
  #print(paste('Running',length(files.cmd),'Ecospace simulations'))
  t1 <- Sys.time()
  runLL.out <- foreach(i = seq_along(files.cmd), .errorhandling = 'pass') %dopar% { #, .options.snow = opts
    #fn.runEwE(cmdfile=files.cmd[i], do.obj=obj.fxn)
    #i=1
    safe_runEwE(cmdfile=files.cmd[i], do.obj = obj.fxn)
    #fn.runEwE(cmdfile=files.cmd[i], do.obj=obj.fxn)
  }
  runLL <- unlist(runLL.out)
  #close(pbar)
  #print(paste('Run time',round(as.numeric(Sys.time()-t1),2)))
  
  ##missing runs----------------------------------------------------------------
  #tmp.errs = sample(1:340,10)
  #runLL[11:15] <- NA_real_ 
  erruns <- which(is.na(runLL))
  #filecheck <- sapply(dirname(files.cmd),FUN=function(x)length(list.files(x)))
  #erruns <- which(filecheck<=1)  
  tries <- 0
  while(length(erruns)>=1 && tries<max_tries){
    tries <- tries+1
    #print(paste0('Redo missing runs: n=',length(erruns)))
    
    #pbar <- winProgressBar("Running Ecospace GApop: Missing Runs",label=paste0("Simulation 0 of ",length(erruns)),max=100)
    #prog <- function(n) setWinProgressBar(pbar,(n/length(erruns)*100),label=paste("Simulation Run", n,"of", length(erruns),"Completed"))
    #opts <- list(progress=prog)
    
    runLL.err.out <- foreach(i=1:length(erruns),.errorhandling='pass') %dopar% {#,.options.snow=opts
      safe_runEwE(cmdfile=files.cmd[i], do.obj = obj.fxn)
      #fn.runEwE(cmdfile=files.cmd[i], do.obj=obj.fxn)
    }
    runLL.err <- unlist(runLL.err.out)
    #close(pbar)
    
    #for(k in 1:length(erruns)) runs[[erruns[k]]] <- unlist(runs.erruns[k])
    for (k in seq_along(erruns)) runLL[erruns[k]] <- runLL.err[k]
    
    #filecheck <- sapply(dirname(files.cmd),FUN=function(x)length(list.files(x)))
    #erruns <- which(filecheck<=1)  #if there are many missing runs, then need to do this in parallel
    #runLL[16] <- NA_real_
    erruns <- which(is.na(runLL))
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

#' @keywords internal
#' @noRd
ensure_cluster <- function() {
  if (!cluster_is_ok()) {
    # rebuild
    try(silent = TRUE, stopCluster(cl))
    rm(list = "cl", envir = .GlobalEnv)
    gc()
    closeAllConnections()
    workers <- max(1L, floor(detectCores() / 4))
    cl <<- parallel::makePSOCKcluster(workers, outfile = "cluster_workers.log")
    doParallel::registerDoParallel(cl)
    clusterExport(cl, c(
      # everything used from worker side:
      "safe_runEwE", "fn.runEwE", "fn.objfxn1", "fn.objfxn2",
      "styear", "enyear", "group.names", "df.names", "obs.ts",
      "cmd_base", "myconfig", "run_dir"
    ))
  }
}
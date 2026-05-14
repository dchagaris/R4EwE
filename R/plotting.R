#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title Plot Ecosim predicted timeseries.
#' @description Makes a multipanel plot of 1 or more ecopace runs with observed data where available.
#' @param predB An array of biomass predictions, as produced by fn.ecospace_ts2array()
#' @param predC An array of catch predictions, as produced by fn.ecospace_ts2array()
#' @param obs.ts A list object containing reference timeseries information, as that created by 
#' fn.read_ecosim_timeseries().
#' @param plt.obs Logical. Should observed data in obs.ts be plotted?
#' @param timestep Read 'annual' or 'monthly' data.
#' @param scale2run If multiple models, which to scale the observed values to.
#' @param pltB.dims Numeric vector of length 2 specifying the number of rows and columns for multipanel plotting.
#' @param pltC.dims Numeric vector of length 2 specifying the number of rows and columns for multipanel plotting.
#' @param scaleCatch Logical. TRUE will scale the predictions so the mean is equal to that of the observed.  Default is FALSE.
#' @param plt.cols Vector of line colors.
#' @param plot2pdf Logical.  TRUE (default) will save files to the specified directory.
#' @param dir.plts  Location to save pdf plots. Defaults to dir.pred.
#' @return Generate a figure in the active plotting device or saves as pdf file.
#' @examples
#' # example code:
#' \dontrun{result <- fn.runEwE.parallel(runlist=myrunlist, obj.fxn=1, cl.export=list('myrunlist','obs.ts'))}
#' @export
fn.ecosim_plot_ts <- function(predB=predB, predC=predC, timestep='annual',obs.ts=obs.ts, scale2run=1, pltB.dims=c(3,3), pltC.dims=c(3,3), 
                                scaleCatch=FALSE,plt.cols=1:dim(predB)[3], plt.obs=TRUE, sim.labels=NULL,
                                dir.plts = dir.pred, plot2pdf=FALSE, run.label=Sys.Date()){
  
  # timestep='annual'
  # scale2run=1
  # pltB.dims=c(4,3)
  # pltC.dims=c(4,3)
  # scaleCatch=TRUE
  # plt.cols=1:dim(predB)[3]
  # dir.plts = dir.pred
  # plot2pdf=FALSE
  
  if(is.null(sim.labels)) sim.labels = dimnames(predB)[[3]]
  #biomass----
  xtime = as.numeric(dimnames(predB)[[1]])
  if(plot2pdf) pdf(file.path(dir.plts,paste0("biomass timeseries fits_",run.label,".pdf")), onefile=T)
  par(mfrow=pltB.dims, mar=c(2,4,2,1), oma=c(4,0,0,1),xpd=F)
  
  for(s in 1:dim(predB)[2]){
    #s=2
    par(xpd=F)
    has.obs = ifelse(s%in%obs.ts$obsB.head$Poolcode,TRUE,FALSE)
    pred.s = predB[,s,]
    pred.base.s = predB[,s,scale2run]
    obs.s=NULL
    plt.ylims = c(0,max(pred.s,na.rm=T)*1.2)
    if(has.obs & plt.obs){
      obs.s = as.matrix(obs.ts$obsB[,which(obs.ts$obsB.head$Poolcode==s)])
      colnames(obs.s) = obs.ts$obsB.head$Title[which(obs.ts$obsB.head$Poolcode==s)]
      #rescale obs to EwE units
      obs.q = colMeans(obs.s, na.rm=T)/mean(pred.base.s, na.rm=T)  
      obs.scaled.s = sweep(obs.s,2,obs.q,"/")
      plt.ylims = c(0,max(pred.s,obs.scaled.s,na.rm=T)*1.2)
    }
    matplot(xtime,pred.s, type='l', lty=1, lwd=2, ylim=plt.ylims, main=dimnames(predB)[[2]][s], col=plt.cols, xlab='',ylab='biomass')
    lines(xtime,pred.base.s,lwd=2,lty=1,col='black')
    if(has.obs & plt.obs){
      legend('topleft',legend=colnames(obs.scaled.s), pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=.8)    
      matpoints(xtime,obs.scaled.s, type='p', pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=0.75)
    }
    
    if(s%in%c(prod(pltB.dims)*(1:ceiling(dim(predB)[2]/prod(pltB.dims))),dim(predB)[2])){
      
      # --- Draw the legend in the bottom outer margin ---
      par(xpd = NA)  # allow drawing outside current panel (no clipping)
      
      # Convert normalized device coordinates (0..1) to current user coords
      # Here we want bottom center: x = 0.5, y = a bit above 0
      x_ndc <- 0.5
      y_ndc <- 0.02   # tweak this if needed (0 = very bottom of device)
      
      x_dev <- grconvertX(x_ndc, from = "ndc", to = "user")
      y_dev <- grconvertY(y_ndc, from = "ndc", to = "user")
      
      legend(x = x_dev, y = y_dev, legend = sim.labels, lty=1, col=plt.cols, bty = "n", xpd = NA, lwd=2,
             xjust = 0.5, yjust = 0, ncol=4)
      
    }
  }
  if(plot2pdf) dev.off()
  
  #catch
  if(plot2pdf) pdf(file.path(dir.plts,paste0("catch timeseries fits_",run.label,".pdf")), onefile=T)
  par(mfrow=pltC.dims, mar=c(2,4,2,1))
  for(s in 1:dim(predC)[2]){
    #s=1
    #pred.name = paste(unique(strsplit(dimnames(predC)[[2]][s],split="\\.")[[1]][-c(1,2)],fromLast=T),collapse="_")
    #ss = match(pred.name, gsub("\\+","",gsub("-","_",df.names$group.names)))
    #pred.name = gsub(" ","_",dimnames(predC.agg)[[2]][s])
    ss = s #match(pred.name, df.names$group.names)
    pred.s = predC[,s,]
    has.pred = ifelse(sum(pred.s)>0,TRUE,FALSE)
    has.obs = ifelse(ss%in%obs.ts$obsC.head$Poolcode,TRUE,FALSE)
    pred.base.s = predC[,s,scale2run]
    obs.s=NULL
    plt.ylims = c(0,max(pred.s,na.rm=T)*1.2)
    if(has.pred & has.obs & plt.obs){
      obs.s = as.matrix(obs.ts$obsC[,which(obs.ts$obsC.head$Poolcode==ss)])
      obs.type.s = obs.ts$obsC.head$Type[which(obs.ts$obsC.head$Poolcode==ss)]
      colnames(obs.s) = obs.ts$obsC.head$Title[which(obs.ts$obsC.head$Poolcode==ss)]
      #rescale obs to EwE units
      if(scaleCatch){ obs.q = colMeans(obs.s, na.rm=T)/mean(pred.base.s, na.rm=T) 
      obs.q[obs.type.s==6] = 1
      }
      if(!scaleCatch) obs.q = rep(1,ncol(obs.s))
      obs.scaled.s = sweep(obs.s,2,obs.q,"/")
      plt.ylims = c(0,max(pred.s,obs.scaled.s,na.rm=T)*1.2)
    }
    if(has.pred){
      matplot(xtime,pred.s, type='l', lty=1, lwd=2, ylim=plt.ylims, main=dimnames(predC)[[2]][s], col=plt.cols, xlab='',ylab='catch')
      lines(xtime,pred.base.s,lty=1,lwd=2,col='black')
    }
    if(has.pred & has.obs & plt.obs){
      legend('topleft',legend=colnames(obs.scaled.s), pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=.8)    
      matpoints(xtime,obs.scaled.s, type='p', pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=0.75)
    } 
    
    if(s%in%c(prod(pltC.dims)*(1:ceiling(dim(predC)[2]/prod(pltC.dims))),dim(predC)[2])){
      
      # --- Draw the legend in the bottom outer margin ---
      par(xpd = NA)  # allow drawing outside current panel (no clipping)
      
      # Convert normalized device coordinates (0..1) to current user coords
      # Here we want bottom center: x = 0.5, y = a bit above 0
      x_ndc <- 0.5
      y_ndc <- 0.02   # tweak this if needed (0 = very bottom of device)
      
      x_dev <- grconvertX(x_ndc, from = "ndc", to = "user")
      y_dev <- grconvertY(y_ndc, from = "ndc", to = "user")
      
      legend(x = x_dev, y = y_dev, legend = sim.labels, lty=1, col=plt.cols, bty = "n", xpd = NA, lwd=2,
             xjust = 0.5, yjust = 0, ncol=4)
      
    }
  }
  if(plot2pdf) dev.off()
  
}#eof





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title Plot Ecospace predicted timeseries.
#' @description Makes a multipanel plot of 1 or more ecopace runs with observed data where available.
#' @param predB An array of biomass predictions, as produced by fn.ecospace_ts2array()
#' @param predC An array of catch predictions, as produced by fn.ecospace_ts2array()
#' @param obs.ts A list object containing reference timeseries information, as that created by 
#' fn.read_ecosim_timeseries().
#' @param plt.obs Logical. Should observed data in obs.ts be plotted?
#' @param timestep Read 'annual' or 'monthly' data.
#' @param scale2run If multiple models, which to scale the observed values to.
#' @param pltB.dims Numeric vector of length 2 specifying the number of rows and columns for multipanel plotting.
#' @param pltC.dims Numeric vector of length 2 specifying the number of rows and columns for multipanel plotting.
#' @param scaleCatch Logical. TRUE will scale the predictions so the mean is equal to that of the observed.  Default is FALSE.
#' @param plt.cols Vector of line colors.
#' @param plot2pdf Logical.  TRUE (default) will save files to the specified directory.
#' @param dir.plts  Location to save pdf plots. Defaults to dir.pred.
#' @return Generate a figure in the active plotting device or saves as pdf file.
#' @examples
#' # example code:
#' \dontrun{result <- fn.runEwE.parallel(runlist=myrunlist, obj.fxn=1, cl.export=list('myrunlist','obs.ts'))}
#' @export
fn.ecospace_plot_ts <- function(predB=predB, predC=predC, timestep='annual',obs.ts=obs.ts, scale2run=1, pltB.dims=c(3,3), pltC.dims=c(1,1), 
                                scaleCatch=FALSE,plt.cols=1:dim(predB)[3], plt.obs=TRUE,
                                dir.plts = dir.pred, plot2pdf=TRUE, run.label=Sys.Date()){
  
  # timestep='annual'
  # scale2run=1
  # pltB.dims=c(4,3)
  # pltC.dims=c(4,3)
  # scaleCatch=TRUE
  # plt.cols=1:dim(predB)[3]
  # dir.plts = dir.pred
  # plot2pdf=FALSE
  
  xtime = as.numeric(dimnames(predB)[[1]])
  if(plot2pdf) pdf(file.path(dir.plts,paste0("biomass timeseries fits_",run.label,".pdf")), onefile=T)
  par(mfrow=pltB.dims, mar=c(2,4,2,1), oma=c(4,0,0,1),xpd=F)
  
  for(s in 1:dim(predB)[2]){
    #s=2
    par(xpd=F)
    has.obs = ifelse(s%in%obs.ts$obsB.head$Poolcode,TRUE,FALSE)
    pred.s = predB[,s,]
    pred.base.s = predB[,s,scale2run]
    obs.s=NULL
    plt.ylims = c(0,max(pred.s,na.rm=T)*1.2)
    if(has.obs & plt.obs){
      obs.s = as.matrix(obs.ts$obsB[,which(obs.ts$obsB.head$Poolcode==s)])
      colnames(obs.s) = obs.ts$obsB.head$Title[which(obs.ts$obsB.head$Poolcode==s)]
      #rescale obs to EwE units
      obs.q = colMeans(obs.s, na.rm=T)/mean(pred.base.s, na.rm=T)  
      obs.scaled.s = sweep(obs.s,2,obs.q,"/")
      plt.ylims = c(0,max(pred.s,obs.scaled.s,na.rm=T)*1.2)
    }
    matplot(xtime,pred.s, type='l', lty=1, lwd=2, ylim=plt.ylims, main=dimnames(predB)[[2]][s], col=plt.cols, xlab='',ylab='biomass')
    if(has.obs & plt.obs){
      legend('topleft',legend=colnames(obs.scaled.s), pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=.8)    
      matpoints(xtime,obs.scaled.s, type='p', pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=0.75)
    }
    
    if(s%in%c(prod(pltB.dims)*(1:ceiling(dim(predB)[2]/prod(pltB.dims))),dim(predB)[2])){
      #add legend
    }
  }
  if(plot2pdf) dev.off()
  
  
  if(plot2pdf) pdf(file.path(dir.plts,paste0("catch timeseries fits_",run.label,".pdf")), onefile=T)
  par(mfrow=pltC.dims, mar=c(2,4,2,1))
  predC.agg <- fn.agg_catch_by_group(predC)
  dimnames(predC.agg)
  for(s in 1:dim(predC.agg)[2]){
    #s=26
    #pred.name = paste(unique(strsplit(dimnames(predC)[[2]][s],split="\\.")[[1]][-c(1,2)],fromLast=T),collapse="_")
    #ss = match(pred.name, gsub("\\+","",gsub("-","_",df.names$group.names)))
    pred.name = gsub(" ","_",dimnames(predC.agg)[[2]][s])
    ss = match(pred.name, df.names$group.names)
    has.obs = ifelse(ss%in%obs.ts$obsC.head$Poolcode,TRUE,FALSE)
    pred.s = predC.agg[,s,]
    pred.base.s = predC.agg[,s,scale2run]
    obs.s=NULL
    plt.ylims = c(0,max(pred.s,na.rm=T)*1.2)
    if(has.obs & plt.obs){
      obs.s = as.matrix(obs.ts$obsC[,which(obs.ts$obsC.head$Poolcode==ss)])
      obs.type.s = obs.ts$obsC.head$Type[which(obs.ts$obsC.head$Poolcode==ss)]
      colnames(obs.s) = obs.ts$obsC.head$Title[which(obs.ts$obsC.head$Poolcode==ss)]
      #rescale obs to EwE units
      if(scaleCatch){ obs.q = colMeans(obs.s, na.rm=T)/mean(pred.base.s, na.rm=T) 
      obs.q[obs.type.s==6] = 1
      }
      if(!scaleCatch) obs.q = rep(1,ncol(obs.s))
      obs.scaled.s = sweep(obs.s,2,obs.q,"/")
      plt.ylims = c(0,max(pred.s,obs.scaled.s,na.rm=T)*1.2)
    }
    matplot(xtime,pred.s, type='l', lty=1, lwd=2, ylim=plt.ylims, main=pred.name, col=plt.cols, xlab='',ylab='catch')
    if(has.obs & plt.obs){
      legend('topleft',legend=colnames(obs.scaled.s), pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=.8)    
      matpoints(xtime,obs.scaled.s, type='p', pch=1:ncol(obs.scaled.s), col=1:ncol(obs.scaled.s), cex=0.75)
    }
  }
  if(plot2pdf) dev.off()

}#eof


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helpers for fn.ecosim_plot_fits ---------------------------------------------------------

#' @keywords internal
#' @noRd
.find_year_skip <- function(x, timestep){
  L <- readLines(x, n = 30L)
  if(timestep == "annual")  n <- which(tolower(substr(L, 1, 4)) == "year")[1] - 1
  if(timestep == "monthly") n <- which(tolower(substr(L, 1, 8)) == "timestep")[1] - 1
  n
}

#' @keywords internal
#' @noRd
.read_pred_wide <- function(dir.out, pattern, timestep){
  files <- unique(unlist(lapply(dir.out, function(d)
    list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE))))
  if(length(files) == 0) return(NULL)
  data_list <- lapply(files, function(x){
    nskip <- .find_year_skip(x, timestep)
    utils::read.csv(x, skip = nskip, row.names = 1, check.names = FALSE)
  })
  yrs  <- rownames(data_list[[1]])
  grps <- names(data_list[[1]])
  arr <- array(NA_real_, dim = c(length(yrs), length(grps), length(files)),
               dimnames = list(yrs, grps, basename(dirname(files))))
  for(r in seq_along(data_list)) arr[,,r] <- as.matrix(data_list[[r]])
  arr
}

#' @keywords internal
#' @noRd
.read_pred_long_group <- function(dir.out, pattern, timestep, agg_fun = sum){
  files <- unique(unlist(lapply(dir.out, function(d)
    list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE))))
  if(length(files) == 0) return(NULL)
  data_list <- lapply(files, function(x){
    nskip <- .find_year_skip(x, timestep)
    d <- utils::read.csv(x, skip = nskip, stringsAsFactors = FALSE)
    mat <- tapply(d$value, list(year = d$year, group = d$group), agg_fun)
    mat[is.na(mat)] <- 0
    mat <- mat[order(as.integer(rownames(mat))), order(as.integer(colnames(mat))), drop = FALSE]
    mat
  })
  yrs  <- sort(unique(unlist(lapply(data_list, rownames))))
  grps <- as.character(sort(as.integer(unique(unlist(lapply(data_list, colnames))))))
  arr <- array(0, dim = c(length(yrs), length(grps), length(files)),
               dimnames = list(yrs, grps, basename(dirname(files))))
  for(r in seq_along(data_list))
    arr[rownames(data_list[[r]]), colnames(data_list[[r]]), r] <- data_list[[r]]
  arr
}

#' @keywords internal
#' @noRd
.read_pred_long_fleetgroup <- function(dir.out, pattern, timestep, agg_fun = sum){
  files <- unique(unlist(lapply(dir.out, function(d)
    list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE))))
  if(length(files) == 0) return(NULL)
  parsed <- lapply(files, function(x){
    nskip <- .find_year_skip(x, timestep)
    d <- utils::read.csv(x, skip = nskip, stringsAsFactors = FALSE)
    d$fg <- paste0("F", d$fleet, "_G", d$group)
    mat <- tapply(d$value, list(year = d$year, fg = d$fg), agg_fun)
    mat[is.na(mat)] <- 0
    meta <- unique(d[, c("fleet", "group", "fg")])
    list(mat = mat, meta = meta)
  })
  yrs <- sort(unique(unlist(lapply(parsed, function(p) rownames(p$mat)))))
  meta <- unique(do.call(rbind, lapply(parsed, function(p) p$meta)))
  meta <- meta[order(meta$fleet, meta$group), ]
  fgs <- meta$fg
  arr <- array(0, dim = c(length(yrs), length(fgs), length(files)),
               dimnames = list(yrs, fgs, basename(dirname(files))))
  for(r in seq_along(parsed))
    arr[rownames(parsed[[r]]$mat), colnames(parsed[[r]]$mat), r] <- parsed[[r]]$mat
  attr(arr, "fleets") <- setNames(meta$fleet, meta$fg)
  attr(arr, "groups") <- setNames(meta$group, meta$fg)
  arr
}

#' @keywords internal
#' @noRd
.read_pred_discards_group <- function(dir.out, timestep){
  catch_grp <- .read_pred_long_group(dir.out, paste0("^catch-fleet-group_", timestep, "\\.csv$"), timestep)
  land_grp  <- .read_pred_long_group(dir.out, paste0("^landings_",          timestep, "\\.csv$"), timestep)
  if(is.null(catch_grp) || is.null(land_grp)) return(NULL)
  yrs  <- sort(unique(c(dimnames(catch_grp)[[1]], dimnames(land_grp)[[1]])))
  grps <- as.character(sort(as.integer(unique(c(dimnames(catch_grp)[[2]], dimnames(land_grp)[[2]])))))
  runs <- unique(c(dimnames(catch_grp)[[3]], dimnames(land_grp)[[3]]))
  expand <- function(a){
    out <- array(0, dim = c(length(yrs), length(grps), length(runs)),
                 dimnames = list(yrs, grps, runs))
    out[dimnames(a)[[1]], dimnames(a)[[2]], dimnames(a)[[3]]] <- a
    out
  }
  arr <- expand(catch_grp) - expand(land_grp)
  arr[arr < 0] <- 0
  arr
}

#' @keywords internal
#' @noRd
.read_pred_discards_fleetgroup <- function(dir.out, timestep){
  catch_fg <- .read_pred_long_fleetgroup(dir.out, paste0("^catch-fleet-group_", timestep, "\\.csv$"), timestep)
  land_fg  <- .read_pred_long_fleetgroup(dir.out, paste0("^landings_",          timestep, "\\.csv$"), timestep)
  if(is.null(catch_fg) || is.null(land_fg)) return(NULL)
  fgs_union <- unique(c(dimnames(catch_fg)[[2]], dimnames(land_fg)[[2]]))
  parsed <- regmatches(fgs_union, regexec("^F(\\d+)_G(\\d+)$", fgs_union))
  fleets_n <- as.integer(sapply(parsed, `[`, 2))
  groups_n <- as.integer(sapply(parsed, `[`, 3))
  ord <- order(fleets_n, groups_n)
  fgs  <- fgs_union[ord]; fleets_n <- fleets_n[ord]; groups_n <- groups_n[ord]
  yrs  <- sort(unique(c(dimnames(catch_fg)[[1]], dimnames(land_fg)[[1]])))
  runs <- unique(c(dimnames(catch_fg)[[3]], dimnames(land_fg)[[3]]))
  expand <- function(a){
    out <- array(0, dim = c(length(yrs), length(fgs), length(runs)),
                 dimnames = list(yrs, fgs, runs))
    out[dimnames(a)[[1]], dimnames(a)[[2]], dimnames(a)[[3]]] <- a
    out
  }
  arr <- expand(catch_fg) - expand(land_fg)
  arr[arr < 0] <- 0
  attr(arr, "fleets") <- setNames(fleets_n, fgs)
  attr(arr, "groups") <- setNames(groups_n, fgs)
  arr
}

#' @keywords internal
#' @noRd
.read_pred_dead_surv_discards_fg <- function(dir.out, timestep){
  disc <- .read_pred_discards_fleetgroup(dir.out, timestep)
  dm   <- .read_pred_long_fleetgroup(dir.out, paste0("^discardmortalityfleetgroup_", timestep, "\\.csv$"), timestep)
  ds   <- .read_pred_long_fleetgroup(dir.out, paste0("^discardsurvivalfleetgroup_", timestep, "\\.csv$"), timestep)
  if(is.null(disc) || is.null(dm) || is.null(ds)) return(NULL)

  fgs_union <- unique(c(dimnames(disc)[[2]], dimnames(dm)[[2]], dimnames(ds)[[2]]))
  parsed <- regmatches(fgs_union, regexec("^F(\\d+)_G(\\d+)$", fgs_union))
  fleets_n <- as.integer(sapply(parsed, `[`, 2))
  groups_n <- as.integer(sapply(parsed, `[`, 3))
  ord <- order(fleets_n, groups_n)
  fgs <- fgs_union[ord]; fleets_n <- fleets_n[ord]; groups_n <- groups_n[ord]
  yrs  <- sort(unique(c(dimnames(disc)[[1]], dimnames(dm)[[1]], dimnames(ds)[[1]])))
  runs <- unique(c(dimnames(disc)[[3]], dimnames(dm)[[3]], dimnames(ds)[[3]]))

  expand <- function(a){
    out <- array(0, dim = c(length(yrs), length(fgs), length(runs)),
                 dimnames = list(yrs, fgs, runs))
    out[dimnames(a)[[1]], dimnames(a)[[2]], dimnames(a)[[3]]] <- a
    out
  }
  d_arr <- expand(disc); m_arr <- expand(dm); s_arr <- expand(ds)
  denom <- m_arr + s_arr
  prop  <- ifelse(denom > 0, m_arr / denom, 0)

  dead <- d_arr * prop
  surv <- d_arr * (1 - prop)

  for(arr_name in c("dead", "surv")){
    attr_target <- get(arr_name)
    attr(attr_target, "fleets") <- setNames(fleets_n, fgs)
    attr(attr_target, "groups") <- setNames(groups_n, fgs)
    assign(arr_name, attr_target)
  }
  list(dead = dead, surv = surv)
}

#' @keywords internal
#' @noRd
.compute_obs_dead_surv_per_group <- function(obs.ts){
  ts.head <- obs.ts$ts.head; ts.data <- obs.ts$ts
  d_rows <- which(ts.head$Type == 13)
  m_rows <- which(ts.head$Type == 11)
  if(length(d_rows) == 0 || length(m_rows) == 0)
    return(list(dead = list(), surv = list()))

  to_combine <- list()
  for(i in d_rows){
    fleet <- ts.head$Poolcode[i]; grp <- ts.head$Poolcode2[i]
    disc <- ts.data[, i]; disc[disc < 0] <- NA
    m_match <- m_rows[ts.head$Poolcode[m_rows] == fleet &
                      (ts.head$Poolcode2[m_rows] == grp | ts.head$Poolcode2[m_rows] == 0)]
    if(length(m_match) == 0) next
    mort <- ts.data[, m_match[1]]
    mort[mort < 0 | mort > 1] <- NA
    gkey <- as.character(grp)
    if(is.null(to_combine[[gkey]])) to_combine[[gkey]] <- list(dead = list(), surv = list())
    to_combine[[gkey]]$dead[[length(to_combine[[gkey]]$dead) + 1]] <- disc * mort
    to_combine[[gkey]]$surv[[length(to_combine[[gkey]]$surv) + 1]] <- disc * (1 - mort)
  }

  sum_with_na <- function(mat){
    row_any <- rowSums(!is.na(mat)) > 0
    s <- rowSums(mat, na.rm = TRUE); s[!row_any] <- NA
    s
  }
  group_dead <- list(); group_surv <- list()
  for(gkey in names(to_combine)){
    dm_mat <- do.call(cbind, to_combine[[gkey]]$dead)
    sm_mat <- do.call(cbind, to_combine[[gkey]]$surv)
    group_dead[[gkey]] <- list(vec = sum_with_na(dm_mat), label = paste0("obs_dead_grp", gkey))
    group_surv[[gkey]] <- list(vec = sum_with_na(sm_mat), label = paste0("obs_surv_grp", gkey))
  }
  list(dead = group_dead, surv = group_surv)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helpers for fn.ecospace_plot_fits ------------------------------------------------------

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

#' @keywords internal
#' @noRd
.read_pred_ecospace_wide <- function(dir.out, varname, timestep, regions, styear,
                                     aggregate_by_group = FALSE){
  prefix <- if(timestep == "annual") "Ecospace_Annual_Average_Region_" else "Ecospace_Average_Region_"

  # Collect data: [[run_key]][[region_id]] = data frame (years x groups [x fleet pre-agg])
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

  # Union of group columns across all run/region data frames
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

#' @keywords internal
#' @noRd
.align_by_group <- function(a, b){
  common <- intersect(dimnames(a)[[2]], dimnames(b)[[2]])
  if(length(common) == 0) return(NULL)
  list(a = a[, common, , , drop = FALSE], b = b[, common, , , drop = FALSE])
}

#' @keywords internal
#' @noRd
.read_pred_ecospace_discards <- function(dir.out, timestep, regions, styear){
  catch <- .read_pred_ecospace_wide(dir.out, "Catch",    timestep, regions, styear, aggregate_by_group = TRUE)
  land  <- .read_pred_ecospace_wide(dir.out, "Landings", timestep, regions, styear, aggregate_by_group = TRUE)
  if(is.null(catch) || is.null(land)) return(NULL)
  al <- .align_by_group(catch, land); if(is.null(al)) return(NULL)
  d <- al$a - al$b
  d[d < 0] <- 0
  d
}

#' @keywords internal
#' @noRd
.read_pred_ecospace_F <- function(dir.out, timestep, regions, styear){
  catch <- .read_pred_ecospace_wide(dir.out, "Catch",   timestep, regions, styear, aggregate_by_group = TRUE)
  bio   <- .read_pred_ecospace_wide(dir.out, "Biomass", timestep, regions, styear, aggregate_by_group = FALSE)
  if(is.null(catch) || is.null(bio)) return(NULL)
  al <- .align_by_group(catch, bio); if(is.null(al)) return(NULL)
  f <- al$a / al$b
  f[!is.finite(f)] <- 0
  f
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title Multipanel Ecospace prediction-vs-observation plots, single combined PDF.
#' @description Reads predicted biomass, catch, landings, discards (catch - landings), and an
#'   approximate fishing mortality (catch / biomass) from one or more Ecospace run folders for each
#'   region present, overlays matching observed timeseries from
#'   \code{\link{fn.read_ecosim_timeseries}} on a per-region basis (using the Region row in the obs
#'   header to assign each obs series to a region), and writes all panels to a single combined PDF.
#'   Each panel shows one group with one predicted line per region; observed points are colored to
#'   match the prediction line of the region the obs is labeled with. For relative (non-absolute)
#'   series, obs are scaled to the prediction of their own region by
#'   \code{q = mean(obs[overlap]) / mean(pred[overlap, obs.region])}.
#' @param dir.out Path to an Ecospace run folder containing the
#'   \code{Ecospace_Annual_Average_Region_<n>_<Var>.csv} (or \code{Ecospace_Average_*} for monthly)
#'   outputs, or a parent folder containing multiple run subfolders. May be a vector of folders.
#' @param obs.ts Output of \code{\link{fn.read_ecosim_timeseries}}. Must include the optional
#'   \code{Region} column (auto-populated as \code{NA} when the ts file has no Region header row);
#'   obs series without a matching region are dropped from the overlays.
#' @param vars Which panels to draw. Any subset of
#'   \code{c("biomass","catch","landings","discards","F")}. Default is all five. \code{discards} is
#'   computed as \code{catch - landings}; \code{F} is \code{catch / biomass} per region.
#' @param timestep \code{"annual"} or \code{"monthly"}.
#' @param regions Integer vector of region IDs to plot. \code{NULL} (default) auto-detects all
#'   \code{Region_<n>_Biomass.csv} files found under \code{dir.out}.
#' @param groups One of \code{"with_obs"} (default), \code{"all"}, or an integer vector of group
#'   pool codes.
#' @param styear Calendar start year. \code{NULL} (default) uses the earliest year in
#'   \code{rownames(obs.ts$ts)}.
#' @param scale2run For multi-run output (\code{dir.out} is a parent of multiple Ecospace run
#'   subfolders), which run column is the scaling baseline / sole pred displayed. Default 1.
#' @param plt.dims Panel grid \code{c(rows, cols)}.
#' @param region.colors Vector of line/point colors per region (in detected-region order). Default
#'   \code{seq_along(regions)}.
#' @param region.names Optional character vector parallel to detected regions (sorted by region ID)
#'   used in the bottom-of-page legend. Default \code{"Region <n>"}.
#' @param sim.labels Multi-run labels. Default = run subfolder basenames.
#' @param group.names Character vector indexed by group pool code. Required to map obs (which uses
#'   pool codes) to Ecospace output columns (which use group names). If \code{NULL}, obs overlay is
#'   skipped (predictions still plotted).
#' @param plot2pdf Logical. If TRUE, writes \code{file.path(dir.plts, "ecospace fits_<run.label>.pdf")}.
#' @param dir.plts Directory to write the PDF.
#' @param run.label String appended to the PDF filename. Defaults to today's date.
#' @return Invisibly returns the list of prediction arrays used (each is
#'   \code{[years, groups, regions, runs]}).
#' @examples
#' \dontrun{
#' ts <- fn.read_ecosim_timeseries("ts_mice_v4_discards_ecospace_regions.csv")
#' fn.ecospace_plot_fits(dir.out = "sp00_5min_init", obs.ts = ts,
#'                       group.names = mygroupnames, plot2pdf = TRUE)
#' }
#' @export
fn.ecospace_plot_fits <- function(dir.out, obs.ts,
                                  vars          = c("biomass", "catch", "landings", "discards", "F"),
                                  timestep      = "annual",
                                  regions       = NULL,
                                  groups        = "with_obs",
                                  styear        = NULL,
                                  scale2run     = 1,
                                  plt.dims      = c(3, 3),
                                  region.colors = NULL,
                                  region.names  = NULL,
                                  sim.labels    = NULL,
                                  group.names   = NULL,
                                  plot2pdf      = FALSE,
                                  dir.plts      = dir.out[1],
                                  run.label     = Sys.Date()){

  vars <- tolower(vars)

  # 1. detect regions and styear -------------------------------------------------------------
  if(is.null(regions)) regions <- .detect_ecospace_regions(dir.out, timestep)
  if(length(regions) == 0)
    stop("No Ecospace per-region output files found under: ", paste(dir.out, collapse = ", "))
  if(is.null(styear)){
    yrs_obs <- suppressWarnings(as.numeric(rownames(obs.ts$ts)))
    styear  <- min(yrs_obs[is.finite(yrs_obs)])
    if(!is.finite(styear)) stop("Could not infer styear from obs.ts; please pass styear = <year>.")
  }

  # 2. read prediction arrays ----------------------------------------------------------------
  pred <- list()
  if("biomass"  %in% vars) pred$biomass  <- .read_pred_ecospace_wide   (dir.out, "Biomass",  timestep, regions, styear, aggregate_by_group = FALSE)
  if("catch"    %in% vars) pred$catch    <- .read_pred_ecospace_wide   (dir.out, "Catch",    timestep, regions, styear, aggregate_by_group = TRUE)
  if("landings" %in% vars) pred$landings <- .read_pred_ecospace_wide   (dir.out, "Landings", timestep, regions, styear, aggregate_by_group = TRUE)
  if("discards" %in% vars) pred$discards <- .read_pred_ecospace_discards(dir.out,            timestep, regions, styear)
  if("f"        %in% vars) pred$fmort    <- .read_pred_ecospace_F      (dir.out,            timestep, regions, styear)

  pred <- pred[!vapply(pred, is.null, logical(1))]
  if(length(pred) == 0)
    stop("No prediction files found under: ", paste(dir.out, collapse = ", "))

  # 3. defaults -------------------------------------------------------------------------------
  nruns <- dim(pred[[1]])[4]
  if(is.null(region.colors)) region.colors <- seq_along(regions)
  if(is.null(region.names))  region.names  <- paste0("Region ", regions)
  if(is.null(sim.labels))    sim.labels    <- dimnames(pred[[1]])[[4]]

  # group name -> pool code lookup
  norm_name <- function(x) gsub("\\s+", "_", trimws(gsub('"', "", as.character(x))))
  norm_gn <- if(!is.null(group.names)) norm_name(group.names) else character(0)

  # 4. open combined PDF ----------------------------------------------------------------------
  if(plot2pdf){
    grDevices::pdf(file.path(dir.plts, paste0("ecospace fits_", run.label, ".pdf")),
                   onefile = TRUE, width = 11, height = 8.5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  n_per_page <- prod(plt.dims)

  .pad_page <- function(n_drawn){
    rem <- (n_per_page - (n_drawn %% n_per_page)) %% n_per_page
    if(rem > 0) for(k in seq_len(rem)) graphics::plot.new()
  }

  # 5. per-variable panel routine -------------------------------------------------------------
  panel_by_group_regions <- function(arr, var_label, obs_type_codes, y_lab,
                                     obs_group_via = "Poolcode"){
    if(is.null(arr)) return(invisible(NULL))
    graphics::par(mfrow = plt.dims, mar = c(2, 4, 2, 1), oma = c(4, 0, 2, 1), xpd = FALSE)

    xtime         <- as.numeric(dimnames(arr)[[1]])
    eco_grp_names <- dimnames(arr)[[2]]
    region_labels <- dimnames(arr)[[3]]
    region_ids    <- as.integer(sub("R", "", region_labels))

    # map each Ecospace column to its pool code (NA if group.names not provided / no match)
    col_to_pc <- if(length(norm_gn))
                   match(norm_name(eco_grp_names), norm_gn)
                 else
                   rep(NA_integer_, length(eco_grp_names))

    obs_match_mask <- obs.ts$ts.head$Type %in% obs_type_codes
    obs_groups_col <- if(obs_group_via == "Poolcode2") obs.ts$ts.head$Poolcode2 else obs.ts$ts.head$Poolcode

    # determine columns to plot
    arr_cols <- seq_along(eco_grp_names)
    plot_cols <-
      if(identical(groups, "all")) arr_cols
      else if(identical(groups, "with_obs")){
        if(all(is.na(col_to_pc))) integer(0)  # no mapping, no obs
        else arr_cols[!is.na(col_to_pc) & col_to_pc %in% unique(obs_groups_col[obs_match_mask])]
      } else arr_cols[!is.na(col_to_pc) & col_to_pc %in% as.integer(groups)]

    if(length(plot_cols) == 0){
      graphics::plot.new()
      graphics::title(main = paste("No groups with obs for", var_label))
      .pad_page(1); return(invisible(NULL))
    }

    obs_yrs <- suppressWarnings(as.numeric(rownames(obs.ts$ts)))
    common_idx <- match(obs_yrs, xtime)

    for(idx in seq_along(plot_cols)){
      col <- plot_cols[idx]
      pc  <- col_to_pc[col]

      # pred matrix: rows = years, cols = regions (single run via scale2run)
      pred_g <- arr[, col, , scale2run, drop = FALSE]
      pred_g <- matrix(pred_g, nrow = length(xtime), ncol = length(region_ids),
                       dimnames = list(dimnames(arr)[[1]], region_labels))

      ylim_top <- max(pred_g, na.rm = TRUE) * 1.2
      obs_per_region <- vector("list", length(region_ids))
      obs_labels_per_region <- vector("list", length(region_ids))
      has_obs_any <- FALSE

      if(!is.na(pc)){
        for(reg_idx in seq_along(region_ids)){
          r <- region_ids[reg_idx]
          obs_rows <- which(obs_match_mask & obs_groups_col == pc &
                            !is.na(obs.ts$ts.head$Region) & obs.ts$ts.head$Region == r)
          if(length(obs_rows) == 0) next

          if(obs_group_via == "Poolcode2" && length(obs_rows) > 1){
            om <- as.matrix(obs.ts$ts[, obs_rows, drop = FALSE]); om[om < 0] <- NA
            row_any <- apply(!is.na(om), 1, any)
            osum <- rowSums(om, na.rm = TRUE); osum[!row_any] <- NA
            obs_s <- matrix(osum, ncol = 1,
                            dimnames = list(NULL, paste0("obs_grp", pc, "_R", r)))
            abs_vec <- all(obs.ts$ts.head$Absolute[obs_rows], na.rm = TRUE)
          } else {
            obs_s <- as.matrix(obs.ts$ts[, obs_rows, drop = FALSE]); obs_s[obs_s < 0] <- NA
            colnames(obs_s) <- gsub(" ", "_", obs.ts$ts.head$Title[obs_rows])
            abs_vec <- obs.ts$ts.head$Absolute[obs_rows]
          }

          pred_r <- pred_g[, reg_idx]
          obs_scaled <- obs_s
          for(j in seq_len(ncol(obs_s))){
            q <- 1
            if(!isTRUE(abs_vec[j])){
              obs_nonna <- !is.na(obs_s[, j])
              pi <- common_idx[obs_nonna]; pi <- pi[!is.na(pi)]
              if(length(pi) > 0){
                mo <- mean(obs_s[obs_nonna, j], na.rm = TRUE)
                mp <- mean(pred_r[pi], na.rm = TRUE)
                if(is.finite(mo) && is.finite(mp) && mp > 0) q <- mo / mp
              }
            }
            obs_scaled[, j] <- obs_s[, j] / q
          }
          ylim_top <- max(ylim_top, max(obs_scaled, na.rm = TRUE) * 1.2, na.rm = TRUE)
          obs_per_region[[reg_idx]]        <- obs_scaled
          obs_labels_per_region[[reg_idx]] <- colnames(obs_s)
          has_obs_any <- TRUE
        }
      }
      if(!is.finite(ylim_top) || ylim_top == 0) ylim_top <- 1

      g_name <- if(!is.null(group.names) && !is.na(pc) && pc <= length(group.names))
                  group.names[pc] else eco_grp_names[col]

      graphics::par(xpd = FALSE)
      graphics::matplot(xtime, pred_g, type = "l", lty = 1, lwd = 2,
                        ylim = c(0, ylim_top), main = g_name,
                        col = region.colors, xlab = "", ylab = y_lab)

      if(has_obs_any){
        filled_pch <- c(21, 22, 23, 24, 25)  # circle, square, diamond, tri-up, tri-down — all with fill via bg
        topleft_labels <- character(0); topleft_pch <- integer(0); topleft_col <- character(0)
        for(reg_idx in seq_along(region_ids)){
          if(is.null(obs_per_region[[reg_idx]])) next
          obs_scaled <- obs_per_region[[reg_idx]]
          col_reg <- region.colors[reg_idx]
          for(j in seq_len(ncol(obs_scaled))){
            sym <- filled_pch[((j - 1) %% length(filled_pch)) + 1]
            graphics::points(obs_yrs, obs_scaled[, j], pch = sym,
                             col = "black", bg = col_reg, cex = 1.1, lwd = 0.5)
            topleft_labels <- c(topleft_labels, obs_labels_per_region[[reg_idx]][j])
            topleft_pch    <- c(topleft_pch, sym)
            topleft_col    <- c(topleft_col, col_reg)
          }
        }
        if(length(topleft_labels))
          graphics::legend("topleft", legend = topleft_labels,
                           pch = topleft_pch, pt.bg = topleft_col, col = "black",
                           pt.cex = 1.0, cex = 0.65, bty = "n")
      }

      if(idx == 1 || ((idx - 1) %% n_per_page == 0))
        graphics::mtext(var_label, side = 3, line = 0.5, outer = TRUE, cex = 1.2, font = 2)

      if(idx %% n_per_page == 0 || idx == length(plot_cols)){
        graphics::par(xpd = NA)
        x_dev <- graphics::grconvertX(0.5, from = "ndc", to = "user")
        y_dev <- graphics::grconvertY(0.02, from = "ndc", to = "user")
        graphics::legend(x = x_dev, y = y_dev, legend = region.names, lty = 1, lwd = 2,
                         col = region.colors, bty = "n", xpd = NA,
                         xjust = 0.5, yjust = 0, ncol = min(4, length(region.names)))
      }
    }
    .pad_page(length(plot_cols))
  }

  # 6. draw each requested variable ----------------------------------------------------------
  if("biomass"  %in% vars) panel_by_group_regions(pred$biomass,  "Biomass",                        c(0, 1),    "biomass")
  if("catch"    %in% vars) panel_by_group_regions(pred$catch,    "Catch",                          c(6, 61, -6), "catch")
  if("landings" %in% vars) panel_by_group_regions(pred$landings, "Landings (group, fleets summed)", c(12),     "landings", obs_group_via = "Poolcode2")
  if("discards" %in% vars) panel_by_group_regions(pred$discards, "Discards (group)",               c(19, 20),  "discards")
  if("f"        %in% vars) panel_by_group_regions(pred$fmort,    "Fishing mortality (catch/biomass)", c(4, 104), "F")

  invisible(pred)
}#eof


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title Multipanel Ecosim prediction-vs-observation plots, single combined PDF.
#' @description Reads predicted biomass, catch, landings (group-aggregated and fleet x group),
#'   and fishing mortality from one or more Ecosim output folders, overlays the matching observed
#'   timeseries from \code{\link{fn.read_ecosim_timeseries}}, and writes all panels to a single PDF
#'   (or the active device). For relative (non-absolute) observation series a scalar
#'   \code{q = mean(obs[overlap]) / mean(pred[overlap, scale2run])} is applied so the obs lie on
#'   the prediction scale; absolute series are plotted as-is.
#' @param dir.out Path to an Ecosim run output folder (containing the required csvs directly), or
#'   a parent folder containing one subfolder per run. Files are located with
#'   \code{list.files(..., recursive = TRUE)} so subfolders are picked up automatically and become
#'   colored lines per panel. May be a vector of folders.
#' @param obs.ts Output of \code{\link{fn.read_ecosim_timeseries}}.
#' @param vars Which panels to draw. Any subset of
#'   \code{c("biomass","catch","landings","discards","F")}. Default is all five. When
#'   \code{"discards"} is included and the run folder contains
#'   \code{discardmortalityfleetgroup_*.csv} and \code{discardsurvivalfleetgroup_*.csv}, three
#'   additional sections are produced: dead-discards (group, fleets stacked), surviving-discards
#'   (group, fleets stacked), and a combined dead-vs-surviving view. Observed dead and surviving
#'   are derived from the type-13 (Discards) and type-11 (DiscardMortality) timeseries in
#'   \code{obs.ts}.
#' @param timestep \code{"annual"} or \code{"monthly"}.
#' @param groups One of \code{"with_obs"} (default; only groups with at least one matching obs
#'   series), \code{"all"} (every group in the pred file), or an integer vector of group pool
#'   codes.
#' @param scale2run Which run column is used as the reference for computing the obs scalar
#'   \code{q}. Default 1.
#' @param plt.dims Panel grid as \code{c(rows, cols)}.
#' @param plt.cols Vector of line colors for each run. Default seq_len(nruns).
#' @param sim.labels Legend labels for runs. Default = run subfolder basenames.
#' @param group.names Optional character vector indexed by group pool code; used in panel titles
#'   when provided. Falls back to \code{"Group <n>"}.
#' @param fleet.names Optional character vector indexed by fleet pool code; used in fleet x group
#'   panel titles and in the fleet legend on the F stacked-bar panels. Falls back to
#'   \code{"F<f>"}.
#' @param plot2pdf Logical. If TRUE, opens a PDF at
#'   \code{file.path(dir.plts, "ecosim fits_<run.label>.pdf")} and writes all panels to that single
#'   file.
#' @param dir.plts Directory to write the PDF. Defaults to the first entry of \code{dir.out}.
#' @param run.label String appended to the PDF filename. Defaults to today's date.
#' @return Invisibly returns a list of the predicted arrays used (one entry per variable plotted).
#' @examples
#' \dontrun{
#' ts <- fn.read_ecosim_timeseries("ts_mice_v4_discards.csv")
#' fn.ecosim_plot_fits(dir.out = "ecosim_sim00_init", obs.ts = ts, plot2pdf = TRUE)
#' }
#' @export
fn.ecosim_plot_fits <- function(dir.out, obs.ts,
                                vars        = c("biomass", "catch", "landings", "discards", "F"),
                                timestep    = "annual",
                                groups      = "with_obs",
                                scale2run   = 1,
                                plt.dims    = c(3, 3),
                                plt.cols    = NULL,
                                sim.labels  = NULL,
                                group.names = NULL,
                                fleet.names = NULL,
                                plot2pdf    = FALSE,
                                dir.plts    = dir.out[1],
                                run.label   = Sys.Date()){

  vars <- tolower(vars)

  # 1. read prediction arrays ----------------------------------------------------------------
  pred <- list()
  if("biomass" %in% vars)
    pred$biomass <- .read_pred_wide(dir.out, paste0("^biomass_", timestep, "\\.csv$"), timestep)
  if("catch" %in% vars)
    pred$catch <- .read_pred_wide(dir.out, paste0("^catch_", timestep, "\\.csv$"), timestep)
  if("landings" %in% vars){
    pred$landings_group <- .read_pred_long_group     (dir.out, paste0("^landings_", timestep, "\\.csv$"), timestep)
    pred$landings_fg    <- .read_pred_long_fleetgroup(dir.out, paste0("^landings_", timestep, "\\.csv$"), timestep)
  }
  if("discards" %in% vars){
    pred$discards_group <- .read_pred_discards_group     (dir.out, timestep)
    pred$discards_fg    <- .read_pred_discards_fleetgroup(dir.out, timestep)
    ds_split <- .read_pred_dead_surv_discards_fg(dir.out, timestep)
    if(!is.null(ds_split)){
      pred$dead_disc_fg <- ds_split$dead
      pred$surv_disc_fg <- ds_split$surv
    }
  }
  if("f" %in% vars){
    pred$fmort    <- .read_pred_long_group     (dir.out, paste0("^mort-fleet-group_", timestep, "\\.csv$"), timestep)
    pred$fmort_fg <- .read_pred_long_fleetgroup(dir.out, paste0("^mort-fleet-group_", timestep, "\\.csv$"), timestep)
  }

  pred <- pred[!vapply(pred, is.null, logical(1))]
  if(length(pred) == 0)
    stop("No prediction files found under: ", paste(dir.out, collapse = ", "))

  # 2. defaults derived from the arrays ------------------------------------------------------
  nruns <- dim(pred[[1]])[3]
  if(is.null(plt.cols))   plt.cols   <- seq_len(nruns)
  if(is.null(sim.labels)) sim.labels <- dimnames(pred[[1]])[[3]]

  # 3. open single combined PDF --------------------------------------------------------------
  if(plot2pdf){
    grDevices::pdf(file.path(dir.plts, paste0("ecosim fits_", run.label, ".pdf")),
                   onefile = TRUE, width = 11, height = 8.5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  n_per_page <- prod(plt.dims)

  .compute_q <- function(obs_vec, pred_vec, obs_year_idx_in_pred, is_absolute){
    if(isTRUE(is_absolute)) return(1)
    obs_nonna <- !is.na(obs_vec)
    pi <- obs_year_idx_in_pred[obs_nonna]
    pi <- pi[!is.na(pi)]
    if(length(pi) == 0) return(1)
    mo <- mean(obs_vec[obs_nonna], na.rm = TRUE)
    mp <- mean(pred_vec[pi], na.rm = TRUE)
    if(!is.finite(mo) || !is.finite(mp) || mp <= 0) return(1)
    mo / mp
  }

  .pad_page <- function(n_drawn){
    rem <- (n_per_page - (n_drawn %% n_per_page)) %% n_per_page
    if(rem > 0) for(k in seq_len(rem)) graphics::plot.new()
  }

  # 4. core panel routine for "by group" variables -------------------------------------------
  panel_by_group <- function(arr, var_label, obs_type_codes, y_lab,
                             group_via = "Poolcode"){
    if(is.null(arr)) return(invisible(NULL))
    graphics::par(mfrow = plt.dims, mar = c(2, 4, 2, 1), oma = c(4, 0, 2, 1), xpd = FALSE)
    xtime <- as.numeric(dimnames(arr)[[1]])
    arr_grp_ids <- as.integer(dimnames(arr)[[2]])

    obs_match_mask <- obs.ts$ts.head$Type %in% obs_type_codes
    obs_groups <- if(group_via == "Poolcode2") obs.ts$ts.head$Poolcode2 else obs.ts$ts.head$Poolcode

    plot_ids <- if(identical(groups, "all"))      arr_grp_ids
                else if(identical(groups, "with_obs")) intersect(arr_grp_ids, unique(obs_groups[obs_match_mask]))
                else as.integer(groups)
    plot_ids <- plot_ids[!is.na(plot_ids)]

    if(length(plot_ids) == 0){
      graphics::plot.new()
      graphics::title(main = paste("No groups with obs for", var_label))
      .pad_page(1)
      return(invisible(NULL))
    }

    obs_yrs <- suppressWarnings(as.numeric(rownames(obs.ts$ts)))
    common_idx <- match(obs_yrs, xtime)

    for(idx in seq_along(plot_ids)){
      g <- plot_ids[idx]
      g_col <- match(as.character(g), dimnames(arr)[[2]])
      if(is.na(g_col)) next

      pred_g <- matrix(arr[, g_col, ], nrow = dim(arr)[1], ncol = dim(arr)[3],
                       dimnames = list(dimnames(arr)[[1]], dimnames(arr)[[3]]))
      pred_base <- pred_g[, scale2run]

      obs_rows <- which(obs_match_mask & obs_groups == g)
      has_obs <- length(obs_rows) > 0
      ylims <- c(0, max(pred_g, na.rm = TRUE) * 1.2)
      obs_plot <- NULL; obs_titles <- character(0)

      if(has_obs){
        obs_s <- as.matrix(obs.ts$ts[, obs_rows, drop = FALSE])
        obs_s[obs_s < 0] <- NA
        colnames(obs_s) <- gsub(" ", "_", obs.ts$ts.head$Title[obs_rows])
        abs_vec <- obs.ts$ts.head$Absolute[obs_rows]

        if(group_via == "Poolcode2" && length(obs_rows) > 1){
          row_any <- apply(!is.na(obs_s), 1, any)
          obs_sum <- rowSums(obs_s, na.rm = TRUE)
          obs_sum[!row_any] <- NA
          obs_s <- matrix(obs_sum, ncol = 1,
                          dimnames = list(NULL, paste0("sum_fleets_grp", g)))
          abs_vec <- all(abs_vec, na.rm = TRUE)
        }

        obs_scaled <- obs_s
        for(j in seq_len(ncol(obs_s))){
          q <- .compute_q(obs_s[, j], pred_base, common_idx, abs_vec[j])
          obs_scaled[, j] <- obs_s[, j] / q
        }
        ylims <- c(0, max(c(pred_g, obs_scaled), na.rm = TRUE) * 1.2)
        if(!is.finite(ylims[2])) ylims <- c(0, 1)
        obs_plot <- obs_scaled
        obs_titles <- colnames(obs_s)
      }

      grp_title <- if(!is.null(group.names) && g <= length(group.names))
                     group.names[g] else paste("Group", g)
      graphics::par(xpd = FALSE)
      graphics::matplot(xtime, pred_g, type = "l", lty = 1, lwd = 2, ylim = ylims,
                        main = grp_title, col = plt.cols, xlab = "", ylab = y_lab)
      graphics::lines(xtime, pred_base, lwd = 2, lty = 1, col = "black")
      if(has_obs){
        graphics::legend("topleft", legend = obs_titles,
                         pch = seq_along(obs_titles), col = seq_along(obs_titles),
                         cex = 0.7, bty = "n")
        graphics::matpoints(obs_yrs, obs_plot, type = "p",
                            pch = seq_len(ncol(obs_plot)), col = seq_len(ncol(obs_plot)),
                            cex = 0.75)
      }
      if(idx == 1 || ((idx - 1) %% n_per_page == 0))
        graphics::mtext(var_label, side = 3, line = 0.5, outer = TRUE, cex = 1.2, font = 2)
      if(idx %% n_per_page == 0 || idx == length(plot_ids)){
        graphics::par(xpd = NA)
        x_dev <- graphics::grconvertX(0.5, from = "ndc", to = "user")
        y_dev <- graphics::grconvertY(0.02, from = "ndc", to = "user")
        graphics::legend(x = x_dev, y = y_dev, legend = sim.labels, lty = 1, col = plt.cols,
                         bty = "n", xpd = NA, lwd = 2, xjust = 0.5, yjust = 0,
                         ncol = min(4, length(sim.labels)))
      }
    }
    .pad_page(length(plot_ids))
  }

  # 5. panel routine for fleet x group variables ---------------------------------------------
  panel_by_fleetgroup <- function(arr, var_label, obs_type_codes, y_lab){
    if(is.null(arr)) return(invisible(NULL))
    graphics::par(mfrow = plt.dims, mar = c(2, 4, 2, 1), oma = c(4, 0, 2, 1), xpd = FALSE)
    xtime <- as.numeric(dimnames(arr)[[1]])
    fleets_attr <- attr(arr, "fleets")
    groups_attr <- attr(arr, "groups")
    fg_all <- dimnames(arr)[[2]]

    obs_match_mask <- obs.ts$ts.head$Type %in% obs_type_codes

    plot_fg <-
      if(identical(groups, "all")) fg_all
      else if(identical(groups, "with_obs")){
        if(length(obs_type_codes) == 0 || !any(obs_match_mask)){
          # no obs at fleet x group resolution; fall back to fg pairs with non-zero pred
          keep <- vapply(fg_all, function(fg) any(arr[, fg, ] != 0, na.rm = TRUE), logical(1))
          fg_all[keep]
        } else {
          keep <- vapply(fg_all, function(fg){
            f <- fleets_attr[fg]; g <- groups_attr[fg]
            any(obs_match_mask & obs.ts$ts.head$Poolcode == f & obs.ts$ts.head$Poolcode2 == g)
          }, logical(1))
          fg_all[keep]
        }
      } else fg_all[as.integer(groups_attr) %in% as.integer(groups)]

    if(length(plot_fg) == 0){
      graphics::plot.new()
      graphics::title(main = paste("No fleet x group panels with obs for", var_label))
      .pad_page(1)
      return(invisible(NULL))
    }

    obs_yrs <- suppressWarnings(as.numeric(rownames(obs.ts$ts)))
    common_idx <- match(obs_yrs, xtime)

    for(idx in seq_along(plot_fg)){
      fg <- plot_fg[idx]
      f <- fleets_attr[fg]; g <- groups_attr[fg]
      pred_g <- matrix(arr[, fg, ], nrow = dim(arr)[1], ncol = dim(arr)[3],
                       dimnames = list(dimnames(arr)[[1]], dimnames(arr)[[3]]))
      pred_base <- pred_g[, scale2run]

      obs_rows <- which(obs_match_mask &
                        obs.ts$ts.head$Poolcode  == f &
                        obs.ts$ts.head$Poolcode2 == g)
      has_obs <- length(obs_rows) > 0
      ylims <- c(0, max(pred_g, na.rm = TRUE) * 1.2)
      obs_plot <- NULL; obs_titles <- character(0)

      if(has_obs){
        obs_s <- as.matrix(obs.ts$ts[, obs_rows, drop = FALSE])
        obs_s[obs_s < 0] <- NA
        colnames(obs_s) <- gsub(" ", "_", obs.ts$ts.head$Title[obs_rows])
        abs_vec <- obs.ts$ts.head$Absolute[obs_rows]
        obs_scaled <- obs_s
        for(j in seq_len(ncol(obs_s))){
          q <- .compute_q(obs_s[, j], pred_base, common_idx, abs_vec[j])
          obs_scaled[, j] <- obs_s[, j] / q
        }
        ylims <- c(0, max(c(pred_g, obs_scaled), na.rm = TRUE) * 1.2)
        if(!is.finite(ylims[2])) ylims <- c(0, 1)
        obs_plot <- obs_scaled
        obs_titles <- colnames(obs_s)
      }

      f_lab <- if(!is.null(fleet.names) && f <= length(fleet.names)) fleet.names[f] else paste0("F", f)
      g_lab <- if(!is.null(group.names) && g <= length(group.names)) group.names[g] else paste0("Grp", g)
      grp_title <- paste0(f_lab, ": ", g_lab)
      graphics::par(xpd = FALSE)
      graphics::matplot(xtime, pred_g, type = "l", lty = 1, lwd = 2, ylim = ylims,
                        main = grp_title, col = plt.cols, xlab = "", ylab = y_lab)
      graphics::lines(xtime, pred_base, lwd = 2, lty = 1, col = "black")
      if(has_obs){
        graphics::legend("topleft", legend = obs_titles,
                         pch = seq_along(obs_titles), col = seq_along(obs_titles),
                         cex = 0.7, bty = "n")
        graphics::matpoints(obs_yrs, obs_plot, type = "p",
                            pch = seq_len(ncol(obs_plot)), col = seq_len(ncol(obs_plot)),
                            cex = 0.75)
      }
      if(idx == 1 || ((idx - 1) %% n_per_page == 0))
        graphics::mtext(var_label, side = 3, line = 0.5, outer = TRUE, cex = 1.2, font = 2)
      if(idx %% n_per_page == 0 || idx == length(plot_fg)){
        graphics::par(xpd = NA)
        x_dev <- graphics::grconvertX(0.5, from = "ndc", to = "user")
        y_dev <- graphics::grconvertY(0.02, from = "ndc", to = "user")
        graphics::legend(x = x_dev, y = y_dev, legend = sim.labels, lty = 1, col = plt.cols,
                         bty = "n", xpd = NA, lwd = 2, xjust = 0.5, yjust = 0,
                         ncol = min(4, length(sim.labels)))
      }
    }
    .pad_page(length(plot_fg))
  }

  # 5b. stacked-bar panel: one bar plot per group, fleet contributions stacked, obs overlaid ----
  panel_stacked_by_group <- function(arr_fg, var_label, y_lab,
                                     obs_type_codes = NULL,
                                     obs_group_via  = "Poolcode",
                                     always_scale   = FALSE,
                                     custom_obs     = NULL){
    if(is.null(arr_fg)) return(invisible(NULL))
    use_custom <- !is.null(custom_obs)
    graphics::par(mfrow = plt.dims, mar = c(3, 4, 2, 1), oma = c(4, 0, 2, 1), xpd = FALSE)
    xtime <- as.numeric(dimnames(arr_fg)[[1]])
    fleets_attr <- attr(arr_fg, "fleets")
    groups_attr <- attr(arr_fg, "groups")
    arr_run <- arr_fg[, , scale2run]  # years x fg, single run for bar composition

    all_groups <- sort(unique(groups_attr))
    if(use_custom){
      obs_match_mask <- rep(FALSE, nrow(obs.ts$ts.head))
      obs_groups_col <- rep(NA_real_, nrow(obs.ts$ts.head))
    } else {
      obs_match_mask <- obs.ts$ts.head$Type %in% obs_type_codes
      obs_groups_col <- if(obs_group_via == "Poolcode2") obs.ts$ts.head$Poolcode2 else obs.ts$ts.head$Poolcode
    }

    plot_groups <-
      if(identical(groups, "all")) all_groups
      else if(identical(groups, "with_obs")){
        keep <- vapply(all_groups, function(g){
          has_obs_g <- if(use_custom)
                         !is.null(custom_obs[[as.character(g)]]) &&
                           any(!is.na(custom_obs[[as.character(g)]]$vec))
                       else
                         any(obs_match_mask & obs_groups_col == g)
          has_obs_g || any(arr_run[, groups_attr == g, drop = FALSE] != 0, na.rm = TRUE)
        }, logical(1))
        all_groups[keep]
      } else intersect(all_groups, as.integer(groups))

    if(length(plot_groups) == 0){
      graphics::plot.new()
      graphics::title(main = paste("No groups to plot for", var_label))
      .pad_page(1)
      return(invisible(NULL))
    }

    # union of fleets that contribute non-zero F to any plotted group
    active_fleets <- sort(unique(fleets_attr[apply(arr_run, 2, function(x) any(x != 0, na.rm = TRUE))]))
    if(length(active_fleets) == 0) active_fleets <- sort(unique(fleets_attr))
    fleet_pal <- grDevices::hcl.colors(max(3L, length(active_fleets)), palette = "Set 2")[seq_along(active_fleets)]

    obs_yrs <- suppressWarnings(as.numeric(rownames(obs.ts$ts)))
    common_idx <- match(obs_yrs, xtime)

    for(idx in seq_along(plot_groups)){
      g <- plot_groups[idx]
      fg_cols <- which(groups_attr == g)
      fleets_g <- fleets_attr[fg_cols]

      stack_mat <- matrix(0, nrow = length(active_fleets), ncol = length(xtime),
                          dimnames = list(paste0("F", active_fleets), as.character(xtime)))
      for(j in seq_along(fg_cols)){
        f_idx <- match(fleets_g[j], active_fleets)
        if(!is.na(f_idx)) stack_mat[f_idx, ] <- stack_mat[f_idx, ] + arr_run[, fg_cols[j]]
      }
      total_pred <- colSums(stack_mat)

      if(use_custom){
        co <- custom_obs[[as.character(g)]]
        if(!is.null(co) && any(!is.na(co$vec))){
          has_obs <- TRUE
          obs_series <- list(list(vec = co$vec, abs = TRUE, label = co$label))
        } else { has_obs <- FALSE; obs_series <- list() }
      } else {
        obs_rows <- which(obs_match_mask & obs_groups_col == g)
        has_obs <- length(obs_rows) > 0
        obs_series <- list()
        if(has_obs){
          if(obs_group_via == "Poolcode2" && length(obs_rows) > 1){
            om <- as.matrix(obs.ts$ts[, obs_rows, drop = FALSE])
            om[om < 0] <- NA
            row_any <- apply(!is.na(om), 1, any)
            osum <- rowSums(om, na.rm = TRUE); osum[!row_any] <- NA
            obs_series <- list(list(vec   = osum,
                                    abs   = all(obs.ts$ts.head$Absolute[obs_rows], na.rm = TRUE),
                                    label = paste0("obs_total_grp", g)))
          } else {
            obs_series <- lapply(obs_rows, function(i){
              ov <- obs.ts$ts[, i]; ov[ov < 0] <- NA
              list(vec = ov, abs = isTRUE(obs.ts$ts.head$Absolute[i]),
                   label = obs.ts$ts.head$Title[i])
            })
          }
        }
      }
      ylim_top <- max(total_pred, na.rm = TRUE) * 1.2
      obs_plot_list <- list()
      if(has_obs){
        for(k in seq_along(obs_series)){
          ovec <- obs_series[[k]]$vec
          q <- 1
          if(always_scale || !isTRUE(obs_series[[k]]$abs)){
            obs_nonna <- !is.na(ovec)
            pi <- common_idx[obs_nonna]; pi <- pi[!is.na(pi)]
            if(length(pi) > 0){
              mo <- mean(ovec[obs_nonna], na.rm = TRUE)
              mp <- mean(total_pred[pi], na.rm = TRUE)
              if(is.finite(mo) && is.finite(mp) && mp > 0) q <- mo / mp
            }
          }
          obs_scaled <- ovec / q
          obs_plot_list[[obs_series[[k]]$label]] <- obs_scaled
          ylim_top <- max(ylim_top, max(obs_scaled, na.rm = TRUE) * 1.2, na.rm = TRUE)
        }
      }
      if(!is.finite(ylim_top) || ylim_top == 0) ylim_top <- 1

      g_name <- if(!is.null(group.names) && g <= length(group.names)) group.names[g] else paste("Group", g)

      graphics::par(xpd = FALSE)
      mids <- graphics::barplot(stack_mat, names.arg = xtime, col = fleet_pal, border = NA,
                                ylim = c(0, ylim_top), main = g_name, ylab = y_lab,
                                las = 2, cex.names = 0.55, space = 0.1)

      if(has_obs && length(obs_plot_list)){
        mid_for_year <- setNames(as.numeric(mids), as.character(xtime))
        ox_all <- mid_for_year[as.character(obs_yrs)]
        pch_seq <- seq_along(obs_plot_list)
        col_seq <- seq_along(obs_plot_list)
        for(k in seq_along(obs_plot_list)){
          ovec <- obs_plot_list[[k]]
          valid <- !is.na(ovec) & !is.na(ox_all)
          if(any(valid)){
            graphics::lines (ox_all[valid], ovec[valid], lwd = 1.5, col = col_seq[k])
            graphics::points(ox_all[valid], ovec[valid], pch = pch_seq[k], col = col_seq[k], cex = 0.7)
          }
        }
        graphics::legend("topleft", legend = names(obs_plot_list),
                         pch = pch_seq, col = col_seq, lwd = 1.5,
                         cex = 0.65, bty = "n")
      }

      if(idx == 1 || ((idx - 1) %% n_per_page == 0))
        graphics::mtext(var_label, side = 3, line = 0.5, outer = TRUE, cex = 1.2, font = 2)

      if(idx %% n_per_page == 0 || idx == length(plot_groups)){
        fleet_labs <- if(!is.null(fleet.names))
                        vapply(active_fleets,
                               function(f) if(f <= length(fleet.names)) fleet.names[f] else paste0("F", f),
                               character(1))
                      else paste0("F", active_fleets)
        graphics::par(xpd = NA)
        x_dev <- graphics::grconvertX(0.5, from = "ndc", to = "user")
        y_dev <- graphics::grconvertY(0.02, from = "ndc", to = "user")
        graphics::legend(x = x_dev, y = y_dev, legend = fleet_labs, fill = fleet_pal,
                         border = NA, bty = "n", xpd = NA, xjust = 0.5, yjust = 0,
                         ncol = min(6, length(fleet_labs)))
      }
    }
    .pad_page(length(plot_groups))
  }

  # 5c. combined dead-vs-surviving discards bar (2 segments per year, no fleet decomposition) ----
  panel_disc_combined <- function(dead_fg, surv_fg, var_label, y_lab){
    if(is.null(dead_fg) || is.null(surv_fg)) return(invisible(NULL))
    graphics::par(mfrow = plt.dims, mar = c(3, 4, 2, 1), oma = c(4, 0, 2, 1), xpd = FALSE)
    xtime <- as.numeric(dimnames(dead_fg)[[1]])
    groups_attr <- attr(dead_fg, "groups")
    dead_run <- dead_fg[, , scale2run]
    surv_run <- surv_fg[, , scale2run]
    all_groups <- sort(unique(groups_attr))

    obs_t13_mask <- obs.ts$ts.head$Type == 13

    plot_groups <-
      if(identical(groups, "all")) all_groups
      else if(identical(groups, "with_obs")){
        keep <- vapply(all_groups, function(g){
          any(obs_t13_mask & obs.ts$ts.head$Poolcode2 == g) ||
            any(dead_run[, groups_attr == g, drop = FALSE] != 0, na.rm = TRUE) ||
            any(surv_run[, groups_attr == g, drop = FALSE] != 0, na.rm = TRUE)
        }, logical(1))
        all_groups[keep]
      } else intersect(all_groups, as.integer(groups))

    if(length(plot_groups) == 0){
      graphics::plot.new(); graphics::title(main = paste("No groups to plot for", var_label))
      .pad_page(1); return(invisible(NULL))
    }

    pal2 <- c("dead" = grDevices::adjustcolor("firebrick", alpha.f = 0.8),
              "surv" = grDevices::adjustcolor("steelblue", alpha.f = 0.8))
    obs_yrs <- suppressWarnings(as.numeric(rownames(obs.ts$ts)))

    for(idx in seq_along(plot_groups)){
      g <- plot_groups[idx]
      fg_cols <- which(groups_attr == g)
      tot_dead <- rowSums(dead_run[, fg_cols, drop = FALSE], na.rm = TRUE)
      tot_surv <- rowSums(surv_run[, fg_cols, drop = FALSE], na.rm = TRUE)
      stack_mat <- rbind(dead = tot_dead, surv = tot_surv)
      colnames(stack_mat) <- as.character(xtime)

      obs_rows <- which(obs_t13_mask & obs.ts$ts.head$Poolcode2 == g)
      obs_total <- NULL
      if(length(obs_rows) > 0){
        om <- as.matrix(obs.ts$ts[, obs_rows, drop = FALSE]); om[om < 0] <- NA
        row_any <- apply(!is.na(om), 1, any)
        obs_total <- rowSums(om, na.rm = TRUE); obs_total[!row_any] <- NA
      }

      ylim_top <- max(colSums(stack_mat), na.rm = TRUE) * 1.2
      if(!is.null(obs_total))
        ylim_top <- max(ylim_top, max(obs_total, na.rm = TRUE) * 1.2, na.rm = TRUE)
      if(!is.finite(ylim_top) || ylim_top == 0) ylim_top <- 1

      g_name <- if(!is.null(group.names) && g <= length(group.names)) group.names[g] else paste("Group", g)

      graphics::par(xpd = FALSE)
      mids <- graphics::barplot(stack_mat, names.arg = xtime, col = pal2, border = NA,
                                ylim = c(0, ylim_top), main = g_name, ylab = y_lab,
                                las = 2, cex.names = 0.55, space = 0.1)

      if(!is.null(obs_total)){
        mid_for_year <- setNames(as.numeric(mids), as.character(xtime))
        ox <- mid_for_year[as.character(obs_yrs)]
        valid <- !is.na(obs_total) & !is.na(ox)
        if(any(valid)){
          graphics::lines (ox[valid], obs_total[valid], lwd = 1.5, col = "black")
          graphics::points(ox[valid], obs_total[valid], pch = 1, col = "black", cex = 0.7)
        }
        graphics::legend("topleft", legend = "obs total discards",
                         pch = 1, col = "black", lwd = 1.5, cex = 0.65, bty = "n")
      }

      if(idx == 1 || ((idx - 1) %% n_per_page == 0))
        graphics::mtext(var_label, side = 3, line = 0.5, outer = TRUE, cex = 1.2, font = 2)
      if(idx %% n_per_page == 0 || idx == length(plot_groups)){
        graphics::par(xpd = NA)
        x_dev <- graphics::grconvertX(0.5, from = "ndc", to = "user")
        y_dev <- graphics::grconvertY(0.02, from = "ndc", to = "user")
        graphics::legend(x = x_dev, y = y_dev,
                         legend = c("dead", "surviving"),
                         fill = pal2, border = NA, bty = "n",
                         xpd = NA, xjust = 0.5, yjust = 0, ncol = 2)
      }
    }
    .pad_page(length(plot_groups))
  }

  # 6. draw each requested variable ----------------------------------------------------------
  if("biomass" %in% vars)
    panel_by_group(pred$biomass, "Biomass", c(0, 1), "biomass")
  if("catch" %in% vars)
    panel_by_group(pred$catch, "Catch", c(6, 61, -6), "catch")
  if("landings" %in% vars){
    panel_stacked_by_group(pred$landings_fg, "Landings (group, fleets stacked)",
                           "landings", obs_type_codes = c(12),
                           obs_group_via = "Poolcode2", always_scale = FALSE)
    panel_by_fleetgroup   (pred$landings_fg, "Landings by fleet x group",
                           c(12), "landings")
  }
  if("discards" %in% vars){
    panel_stacked_by_group(pred$discards_fg, "Discards (group, fleets stacked)",
                           "discards", obs_type_codes = c(19, 20),
                           obs_group_via = "Poolcode", always_scale = FALSE)
    panel_by_fleetgroup   (pred$discards_fg, "Discards by fleet x group",
                           c(13), "discards")
    if(!is.null(pred$dead_disc_fg) && !is.null(pred$surv_disc_fg)){
      ds_obs <- .compute_obs_dead_surv_per_group(obs.ts)
      panel_stacked_by_group(pred$dead_disc_fg, "Dead discards (group, fleets stacked)",
                             "dead discards", custom_obs = ds_obs$dead)
      panel_stacked_by_group(pred$surv_disc_fg, "Surviving discards (group, fleets stacked)",
                             "surviving discards", custom_obs = ds_obs$surv)
      panel_disc_combined(pred$dead_disc_fg, pred$surv_disc_fg,
                          "Discards: dead vs surviving", "discards")
    }
  }
  if("f" %in% vars){
    panel_by_group        (pred$fmort, "Fishing mortality (group)", c(4, 104), "F")
    panel_stacked_by_group(pred$fmort_fg, "Fishing mortality stacked by fleet",
                           "F", obs_type_codes = c(4, 104),
                           obs_group_via = "Poolcode", always_scale = TRUE)
  }

  invisible(pred)
}#eof


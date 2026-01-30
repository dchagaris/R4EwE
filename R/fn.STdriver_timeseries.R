#' @title Get ST driver timeseries average.
#' @description Creates a movie (gif) file from ST driver ascii rasters.
#' @param dirs.stdriver A character vector containing the full directory paths to 1 or more ST driver variable folders.
#' @param grid.subset A numeric vector of lenght 4 (xmin,xmax,ymin,ymax) indicating which grid cells subset and then average over. These are row and column indexes. NULL, the default, averages over the entire spatial grid.  A value of c(1,3,1,3) would get the top left 3x3 corner. Need to switch to a grid mask so cells subset can be any shape/cells.
#' @param datestamps Optional character vector of timesteps for each file.  Setting to NULL (default) will pull the timestep from the ascii file suffix that is ...YYYYMM.asc. 
#' @return Saves a STdriver_timeseries.csv file in the parent directoryr  
#' @examples
#' #example code:
#' fn.STdriver_timeseries(dirs.stdriver=list.dirs("./ST drivers/30min")[1:2], grid.subset=NULL, datestamps=NULL)
#' @export
fn.STdriver_timeseries <- function(dirs.stdriver=dirs.stdriver, datestamps=NULL, grid.subset=NULL){
  
  # dirs.stdriver=list.dirs("C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\ST drivers\\30min")[c(4,28)]
  # datestamps=NULL
  # grid.subset = NULL

  #package setup
  required <- c("terra")
  #apply(1:length(required), FUN=function(i)library(required[i]))
  missing  <- required[!vapply(required, requireNamespace, FUN.VALUE = TRUE, quietly = TRUE)]
  if (length(missing)>0) {
    stop(
      sprintf(
        "Missing required packages: %s\nInstall with: install.packages(c(%s))",
        paste(missing, collapse = ", "),
        paste(sprintf('"%s"', missing), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # attach once
  invisible(lapply(required, function(p) {
    if (!p %in% .packages()) library(p, character.only = TRUE)
  }))
  
  for(i in 1:length(dirs.stdriver)){
    #i=1
    var.name=basename(dirs.stdriver[i])
    message(paste(i,'of',length(dirs.stdriver),":",var.name))
    files.st = list.files(dirs.stdriver[i],full.names = T,recursive=T)
    st.rast = rast(files.st)
    names(st.rast) = substr(basename(files.st),nchar(basename(files.st))-9,nchar(basename(files.st))-4)
    
    #unit conversions
    #npp: molC/m2/s to gC/m2/yr
    if(var.name=='wc_vert_int_npp_cefi') st.rast <- st.rast*(365*24*60*60)*12
    
    #grid subset
    if(!is.null(grid.subset)) st.rast = st.rast[grid.subset[1]:grid.subset[2],grid.subset[3]:grid.subset[4],1:nlyr(st.rast),drop=FALSE]
    
    #calculate mean
    st.rast.mean <- global(st.rast, "mean",na.rm=TRUE)
    names(st.rast.mean) <- var.name
    
    #store output
    if(i==1){
      mean.out <- st.rast.mean
    } else {
      mean.out <- merge(mean.out,st.rast.mean,by=0,all.y=T)
    }
    rm(st.rast);gc()
  }
  names(mean.out)[1] <- 'time'
  mean.out$time = as.numeric(substr(mean.out$time,1,4))+(as.numeric(substr(mean.out$time,5,6))-1)/12
  write.csv(mean.out,file.path(dirname(dirs.stdriver[1]),'STdriver_timeseries_mean.csv'),row.names=F)
} #eof 

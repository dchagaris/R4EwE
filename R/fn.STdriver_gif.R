#' @title Make ST driver gifs
#' @description Creates a movie (gif) file from ST driver ascii rasters.
#' @param dir.stdriver Full directory path to the ST driver folder.
#' @param do.files A numeric vector indexing which of the ascii files in dir.stdriver to plot.  The default NULL plots all files in the folder.
#' @param datestamps Optional character vector of timesteps for each file.  Setting to NULL (default) will pull the timestep from the ascii file suffix that is ...YYYYMM.asc. 
#' @return Place a gif file in a folder one level up from dri.stdriver  
#' @examples
#' #example code:
#' fn.STdriver_gif(dir.stdriver="./ST drivers/30min/wc_vert_int_npp_cefi", do.files=NULL, datestamps=NULL)
#' @export
fn.STdriver_gif <- function(dir.stdriver=dir.stdriver, do.files=NULL, datestamps=NULL){
  
  # dir.stdriver=("C:\\Users\\dchagaris\\OneDrive - University of Florida\\Atlantic menhaden\\NWACS-MICE\\Ecospace\\ST drivers\\30min\\wc_vert_int_npp_cefi")
  # do.files=1:12
  # datestamps=NULL

  #package setup
  required <- c("terra", "viridis", "maps", "gifski", "fields")
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
  
  # ---- Input validation ----
  if (missing(dir.stdriver) || is.null(dir.stdriver) || !dir.exists(dir.stdriver)) {
    stop("Provide a valid 'dir.stdriver' folder that contains your raster files.", call. = FALSE)
  }
  
  
  
  dir.gif=file.path(dirname(dir.stdriver),'gifs')
  dir.png = file.path(dirname(dir.stdriver),'png_tmp')
  var.name=basename(dir.stdriver)
  
  if(!dir.exists(dir.gif))dir.create(dir.gif)
  if(!dir.exists(dir.png))dir.create(dir.png)
  
  files.st = list.files(dir.stdriver,full.names = T,recursive=T)
  if(!is.null(do.files))   files.st = files.st[do.files]
  
  #read them into terra stack
  st.rast = rast(files.st)
  names(st.rast) = substr(basename(files.st),nchar(basename(files.st))-9,nchar(basename(files.st))-4)
  
  #unit conversions
  #npp: molC/m2/s to gC/m2/yr
  st.units<-""
  if(var.name=='wc_vert_int_npp_cefi') st.rast <- st.rast*(365*24*60*60)*12
  if(var.name=='wc_vert_int_npp_cefi') st.units <- "gC/m2/yr"
  
  #get limits
  st.rast.min <- min(global(st.rast, "min",na.rm=TRUE))
  st.rast.max <- max(global(st.rast, "max", na.rm=TRUE))
  
  #plot-------------------------------------------------------------------------
  ##get zlim and colors----
  plt.zlims <- c(st.rast.min, st.rast.max)
  plt.brks    = seq(plt.zlims[1],plt.zlims[2],length.out=100)
  plt.col = viridis(n=length(plt.brks)-1)
  
  # Choose fixed legend ticks you want to see on every plot
  tick_at <- pretty(plt.zlims, n = 6)  # ~5–6 ticks; adjust n to taste
  tick_lab <- format(tick_at, digits = 3, scientific = FALSE)  # optional formatting
  
  ##plot pngs
  #graphics.off();rm(.SavedPlots);windows(record=T)
  for(i in 1:dim(st.rast)[3]){
    #i=1
    
    png(file.path(dir.png, paste0(var.name, "_", names(st.rast)[i], ".png")),
        width = 1400, height = 1000, res = 150)
    
    # Create a two-column layout: big map + skinny legend
    layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.15))   # adjust 0.07 to be thinner or wider
    
    ## Panel 1: Map (tight margins)
    par(mar = c(2, .1, 2, 0.1)) 
    plot(st.rast[[i]], col=plt.col, zlim=plt.zlims, main=paste(var.name,names(st.rast)[i]),breaks=plt.brks, legend=FALSE,
         mar=c(c(2, 0, 2, 0)))
    map(database='world', add=T, fill=T, col='gray')
    map(database='state',add=T,fill=F)
    
    ## Panel 2: Legend (continuous bar; tiny margins)
    par(mar = c(.1, 0.1, .1, 1))
    # IMPORTANT: supply zlim for robustness
    fields::image.plot(legend.only = TRUE,
                       zlim       = plt.zlims,
                       col        = plt.col,
                       axis.args  = list(at = tick_at, labels = tick_lab, cex.axis = 0.85),
                       legend.lab = paste0("NPP (",st.units,")"),
                       legend.cex = 0.9,
                       legend.mar = 0.1,
                       horizontal = FALSE)
    dev.off()
  }
  
  #make gif----
  png_files <- list.files(dir.png, pattern = "png$", full.names = TRUE)
  gifski(png_files, gif_file = file.path(dir.gif,paste0(var.name,".gif")), width = 1400, height = 1000, delay = .1)
  unlink(dir.png, recursive=T, force=T)
} #eof 

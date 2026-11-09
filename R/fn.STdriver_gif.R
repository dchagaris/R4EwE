#' @title Make ST driver gifs (deprecated)
#' @description Deprecated. Use \code{\link{fn.ecospace_gif}}, which animates ST driver
#'   folders along with Ecospace output and fixes two defects in this function: frames
#'   were dated by a fixed character offset, which misreads the \code{_YYYY_MM} names of
#'   a monthly climatology and silently interleaves it with the real-dated series, and
#'   the legend label was drawn on top of the legend's tick labels. This wrapper is kept
#'   so existing scripts keep working; it forwards to \code{fn.ecospace_gif()} with
#'   \code{source = 'stdriver'}.
#' @param dir.stdriver Full directory path to the ST driver folder.
#' @param do.files A numeric vector indexing which ascii files in dir.stdriver to plot.
#'   The default NULL plots all of them. Forwarded as \code{do.steps}, so it now indexes
#'   the dated grids in date order rather than whatever \code{list.files()} returned;
#'   an explicit index may therefore select different frames than it used to.
#' @param datestamps Optional character vector of timesteps for each file. Setting to
#'   NULL (default) pulls the timestep from the ascii file name.
#' @return Invisibly, the gif file written, in a 'gifs' folder one level up from
#'   dir.stdriver. (This function previously returned nothing.)
#' @examples
#' # example code:
#' \dontrun{
#' # deprecated:
#' fn.STdriver_gif(dir.stdriver = "./ST drivers/30min/wc_vert_int_npp_cefi")
#' # use instead:
#' fn.ecospace_gif(dir.out = "./ST drivers/30min/wc_vert_int_npp_cefi",
#'                 inset = TRUE)
#' }
#' @seealso \code{\link{fn.ecospace_gif}}
#' @export
fn.STdriver_gif <- function(dir.stdriver=dir.stdriver, do.files=NULL, datestamps=NULL){

  .Deprecated("fn.ecospace_gif", package = "R4EwE",
              msg = paste("fn.STdriver_gif() is deprecated; use fn.ecospace_gif()",
                          "with an ST driver folder as 'dir.out'."))

  # ---- Input validation ----
  if (missing(dir.stdriver) || is.null(dir.stdriver) || !dir.exists(dir.stdriver)) {
    stop("Provide a valid 'dir.stdriver' folder that contains your raster files.",
         call. = FALSE)
  }

  # The old function converted NPP units for one hard-coded folder name. Preserved so
  # existing scripts get the same numbers; fn.ecospace_gif()'s 'transform' argument is
  # the general way to do this for any variable.
  var.name  <- basename(dir.stdriver)
  transform <- NULL
  legend.lab <- NULL
  if (var.name == "wc_vert_int_npp_cefi") {
    transform  <- function(x) x * (365 * 24 * 60 * 60) * 14.007
    legend.lab <- "NPP (gC/m2/yr)"
    message("Applying the legacy NPP unit conversion for '", var.name,
            "'. In fn.ecospace_gif(), pass 'transform' explicitly instead.")
  }

  fn.ecospace_gif(dir.out    = dir.stdriver,
                  source     = "stdriver",
                  do.steps   = do.files,
                  datestamps = datestamps,
                  transform  = transform,
                  legend.lab = legend.lab,
                  width      = 1400,
                  height     = 1000,
                  res        = 150,
                  delay      = 0.1)
} #eof

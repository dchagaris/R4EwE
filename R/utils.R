#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.longvuls ---------------------------------------------------------------------------------------
#' @title Melt vulnerability matrix
#' @description melt a vulnerability matrix.
#' @param vuls A vulnerability matrix
#' @keywords internal
#' @noRd
fn.longvuls <- function(vuls=vuls){
  #vuls=vuls.base
  longvuls <- reshape2::melt(vuls, id.vars=1) ## Melt to long
  names(longvuls) <- c('prey','pred','basevul')  ## Rename
  longvuls$prey <- match(longvuls$prey,vuls$prey)
  longvuls$pred <- match(longvuls$pred,names(vuls)[-1])
  longvuls <- longvuls[!is.na(longvuls$basevul),]
  return(longvuls)
}
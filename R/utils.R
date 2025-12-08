#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.longvuls ---------------------------------------------------------------------------------------
#' @title Melt vulnerability matrix
#' @description melt a vulnerability matrix.
#' @param vuls A vulnerability matrix
#' @keywords internal
#' @noRd
fn.longvuls <- function(vuls=vuls){
  longvuls <- reshape2::melt(vuls, id.vars <- 1) ## Melt to long
  names(longvuls)[2:3] <- c('pred','baseval')  ## Rename
  longvuls$prey <- as.factor(longvuls$prey)    
  longvuls$grp.pred <- as.integer(longvuls$pred)  ## Add group number for pred. Note this doesn't work for prey
  longvuls <- merge(longvuls, df.names, by.x="prey", by.y ="group.names") ## Merge in prey group number from df.names
  names(longvuls)[which(names(longvuls)=='num')] = 'grp.prey'
  #longvuls <- rename(longvuls, c(grp.prey = "num"))
  longvuls <- longvuls[!is.na(longvuls$baseval),]
  return(longvuls)
}
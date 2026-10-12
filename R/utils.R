#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.longvuls ---------------------------------------------------------------------------------------
#' @title Melt vulnerability matrix
#' @description melt a vulnerability matrix.
#' @param vuls A vulnerability matrix
#' @param by \code{"pair"} (default): one row per predator-prey pair.
#'   \code{"pred"}: one row per predator with \code{prey = NA}, which the
#'   console applies to every prey in that predator's column. The base value
#'   is the column's value when uniform, else its geometric mean (with a
#'   warning, since the base model is then not reproduced exactly).
#' @export
fn.longvuls <- function(vuls=vuls, by=c("pair","pred")){
  #vuls=vuls.base
  by <- match.arg(by)
  longvuls <- reshape2::melt(vuls, id.vars=1) ## Melt to long
  names(longvuls) <- c('prey','pred','basevul')  ## Rename
  longvuls$prey <- match(longvuls$prey,vuls$prey)
  longvuls$pred <- match(longvuls$pred,names(vuls)[-1])
  longvuls <- longvuls[!is.na(longvuls$basevul),]
  if(by == "pred"){
    cols <- split(longvuls$basevul, longvuls$pred)
    mixed <- names(cols)[vapply(cols, function(x) length(unique(x)) > 1, logical(1))]
    if(length(mixed))
      warning("Predator column(s) ", paste(mixed, collapse=", "), " hold more than one vulnerability; ",
              "using the geometric mean as the per-predator base value.")
    longvuls <- data.frame(prey = NA_integer_, pred = as.integer(names(cols)),
                           basevul = vapply(cols, function(x) exp(mean(log(x))), numeric(1)))
    rownames(longvuls) <- NULL
  }
  return(longvuls)
}

# Internal: index part of an <ECOSIM_VULNERABILITIES_INDEXED> tag -- "pred prey"
# for a pair, "pred" alone (whole predator column) when prey is NA.
.vul_tag_index <- function(pred, prey) ifelse(is.na(prey), as.character(pred), paste(pred, prey))
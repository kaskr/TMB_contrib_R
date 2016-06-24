
#' Calculate marginal AIC for a fitted model
#'
#' \code{TMBAIC} calculates AIC for a given model fit
#'
#' @param opt, the output from \code{nlminb} or \code{optim}
#' @param k the penalty on additional fixed effects (default=2, for AIC)
#'
#' @return AIC, where a parsimonious model has a AIC relative to other candidate models

#' @export
TMBAIC=function(opt, k=2){
  if( all(c("par","objective") %in% names(opt)) ) Return = 2*length(opt[["par"]]) + k*opt[["objective"]]
  if( all(c("par","value") %in% names(opt)) ) Return = 2*length(opt[["par"]]) + k*opt[["value"]]
  return( Return )
}


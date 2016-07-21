
#' Optimize a TMB model
#'
#' \code{Optimize} runs a TMB model and generates standard diagnostics
#'
#' @param obj The compiled TMB object
#' @param startpar Starting values for fixed effects
#' @param lower lower bounds on fixed effects
#' @param upper upper bounds on fixed effects
#' @param getsd Boolean whether to run standard error calculation
#' @param control list of options to pass to \code{nlminb}
#' @param ... list of settings to pass to \code{sdreport}
#'
#' @return the standard output from \code{nlminb}, except with additional diagnostics and timing info, and a new slot containing the output from \code{sdreport}

#' @examples
#' TMBhelper::Optimize( Obj ) # where Obj is a compiled TMB object

#' @export
Optimize = function( obj, startpar=obj$par, lower=rep(-Inf,length(startpar)), upper=rep(Inf,length(startpar)), getsd=TRUE, control=list(eval.max=1e4, iter.max=1e4, trace=TRUE), ... ){
  # Run
  start_time = Sys.time()
  opt = nlminb( start=startpar, objective=obj$fn, gradient=obj$gr, control=control )

  # Add diagnostics
  opt[["run_time"]] = Sys.time() - start_time
  opt[["diagnostics"]] = data.frame( "Param"=names(obj$par), "starting_value"=startpar, "Lower"=lower, "MLE"=opt$par, "Upper"=upper, "final_gradient"=obj$gr(opt$par) )

  # Get standard deviations
  if(getsd==TRUE) opt[["SD"]] = sdreport( obj, ... )

  # Return stuff
  return( opt )
}


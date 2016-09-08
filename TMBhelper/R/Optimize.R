
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
#' @param savedir directory to save results (if \code{savedir=NULL}, then results aren't saved)
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param ... list of settings to pass to \code{sdreport}
#'
#' @return the standard output from \code{nlminb}, except with additional diagnostics and timing info, and a new slot containing the output from \code{sdreport}

#' @examples
#' TMBhelper::Optimize( Obj ) # where Obj is a compiled TMB object

#' @export
Optimize = function( obj, startpar=obj$par, lower=rep(-Inf,length(startpar)), upper=rep(Inf,length(startpar)), getsd=TRUE, control=list(eval.max=1e4, iter.max=1e4, trace=TRUE),
  savedir=NULL, loopnum=3, ... ){

  # Run first time
  start_time = Sys.time()
  opt = nlminb( start=startpar, objective=obj$fn, gradient=obj$gr, control=control, lower=lower, upper=upper )

  # Re-run to further decrease final gradient
  for( i in seq(2,loopnum,length=max(0,loopnum-1)) ){
    opt = nlminb( start=opt$par, objective=obj$fn, gradient=obj$gr, control=control, lower=lower, upper=upper )
  }

  # Add diagnostics
  opt[["run_time"]] = Sys.time() - start_time
  opt[["number_of_coefficients"]] = c("Total"=length(unlist(obj$env$parameters)), "Fixed"=length(obj$par), "Random"=length(unlist(obj$env$parameters))-length(obj$par) )
  opt[["AIC"]] = 2*opt$objective + 2*length(opt$par)
  opt[["diagnostics"]] = data.frame( "Param"=names(obj$par), "starting_value"=startpar, "Lower"=lower, "MLE"=opt$par, "Upper"=upper, "final_gradient"=obj$gr(opt$par) )

  # Get standard deviations
  if(getsd==TRUE) opt[["SD"]] = sdreport( obj, ... )

  # Save results (excluding 'env' which is too big to routinely archive)
  if( !is.null(savedir) ){
  }

  # Return stuff
  return( opt )
}



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
#' @param newtonsteps number of extra newton steps to take after optimization (alternative to \code{loopnum})
#' @param n sample sizes (if \code{n!=Inf} then \code{n} is used to calculate BIC and AICc)
#' @param ... list of settings to pass to \code{sdreport}
#'
#' @return the standard output from \code{nlminb}, except with additional diagnostics and timing info, and a new slot containing the output from \code{sdreport}

#' @examples
#' TMBhelper::Optimize( Obj ) # where Obj is a compiled TMB object

#' @export
Optimize = function( obj, fn=obj$fn, gr=obj$gr, startpar=obj$par, lower=rep(-Inf,length(startpar)), upper=rep(Inf,length(startpar)),
  getsd=TRUE, control=list(eval.max=1e4, iter.max=1e4, trace=0), bias.correct=FALSE,
  bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct=NULL),
  savedir=NULL, loopnum=3, newtonsteps=0, n=Inf, ... ){

  # Run first time
  start_time = Sys.time()
  opt = nlminb( start=startpar, objective=obj$fn, gradient=obj$gr, control=control, lower=lower, upper=upper )

  # Re-run to further decrease final gradient
  for( i in seq(2,loopnum,length=max(0,loopnum-1)) ){
    opt = nlminb( start=opt$par, objective=obj$fn, gradient=obj$gr, control=control, lower=lower, upper=upper )
  }

  ## Run some Newton steps
  for(i in seq_len(newtonsteps)) {
    g <- as.numeric( obj$gr(opt$par) )
    h <- optimHess(opt$par, obj$fn, obj$gr)
    opt$par <- opt$par - solve(h, g)
    opt$objective <- obj$fn(opt$par)
  }

  # Add diagnostics
  opt[["run_time"]] = Sys.time() - start_time
  opt[["max_gradient"]] = max(abs(obj$gr(opt$par)))
  opt[["number_of_coefficients"]] = c("Total"=length(unlist(obj$env$parameters)), "Fixed"=length(obj$par), "Random"=length(unlist(obj$env$parameters))-length(obj$par) )
  opt[["AIC"]] = TMBhelper::TMBAIC( opt=opt )
  if( n!=Inf ){
    opt[["AICc"]] = TMBhelper::TMBAIC( opt=opt, n=n )
    opt[["BIC"]] = TMBhelper::TMBAIC( opt=opt, p=log(n) )
  }
  opt[["diagnostics"]] = data.frame( "Param"=names(obj$par), "starting_value"=startpar, "Lower"=lower, "MLE"=opt$par, "Upper"=upper, "final_gradient"=as.vector(obj$gr(opt$par)) )

  # Get standard deviations
  if(getsd==TRUE){
    # Compute hessian
    h <- optimHess(opt$par, obj$fn, obj$gr)
    # Check for problems
    if( is.character(try(chol(h),silent=TRUE)) ){
      warning("Hessian is not positive definite, so standard errors are not available")
      return( list("opt"=opt, "h"=h) )
    }
    # Compute standard errors
    if( bias.correct==FALSE | is.null(bias.correct.control[["vars_to_correct"]]) ){
      opt[["SD"]] = sdreport( obj=obj, par.fixed=opt$par, hessian.fixed=h, bias.correct=bias.correct, bias.correct.control=bias.correct.control[c("sd","split","nsplit")], ... )
    }else{
      if( "ADreportIndex" %in% names(obj$env) ){
        Which = as.vector(unlist( Obj$env$ADreportIndex()[ bias.correct.control[["vars_to_correct"]] ] ))
      }else{
        # Run first time to get indices
        opt[["SD"]] = sdreport( obj=obj, par.fixed=opt$par, hessian.fixed=h, bias.correct=FALSE, ... )
        # Determine indices
        Which = which( rownames(summary(opt[["SD"]],"report")) %in% bias.correct.control[["vars_to_correct"]] )
      }
      # Split up indices
      if(bias.correct.control[["nsplit"]]>1) Which = split( Which, cut(seq_along(Which), bias.correct.control[["nsplit"]]) )
      Which = Which[sapply(Which,FUN=length)>0]
      if(length(Which)==0) Which = NULL
      # Repeat SD with indexing
      message( paste0("Bias correcting ", length(Which), " derived quantities") )
      opt[["SD"]] = sdreport( obj=obj, par.fixed=opt$par, hessian.fixed=h, bias.correct=TRUE, bias.correct.control=list(sd=bias.correct.control[["sd"]], split=Which, nsplit=NULL), ... )
    }
  }

  # Save results
  if( !is.null(savedir) ){
    parameter_estimates = opt
    #parameter_estimates$SD = parameter_estimates$SD[ setdiff(names(parameter_estimates$SD),"env") ]
    save( parameter_estimates, file=file.path(savedir,"parameter_estimates.RData"))
    capture.output( parameter_estimates, file=file.path(savedir,"parameter_estimates.txt"))
  }

  # Return stuff
  return( opt )
}


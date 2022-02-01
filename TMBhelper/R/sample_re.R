
#' Sample random effects to correct for re-transformation bias
#'
#' \code{sample_re} calculates MCMC samples of random effects conditional upon
#' estimated MLE for fixed effects, and then uses each sample to calculate objects in the report.
#' This is useful e.g., in correcting for re-transformation bias (by calculating the
#' posterior mean of a nonlinear transformation of random effects) or visualizing
#' random-effect variance (which often can be time-consuming using the delta-method
#' in models with many random effects)
#'
#' @inheritParams rstan::sampling
#' @param obj the TMB object after parameter estimation
#' @param report_names which elements of \code{obj$report()} should be recorded;
#'        default \code{report_names=NULL} uses \code{report_names=names(obj$report())}
#' @param ... adding arguments to pass to \code{\link[tmbstan]{tmbstan}}
#'
#' @return A tagged list containing:
#' \describe{
#'   \item{\code{stan_out}}{output from \code{\link[tmbstan]{tmbstan}}}
#'   \item{\code{report_full}}{A list of output from \code{obj$report()[report_names], except with extra dimension for each MCMC sample}}
#'   \item{\code{run_time}}{total run time}
#' }
#'
#' @references For a discussion of the epsilon-estimator as alternative method
#'             to correct for re-transformation bias see \url{https://doi.org/10.1016/j.fishres.2015.11.016}
#' @export
sample_re <-
function( obj,
          warmup = 50,
          iter = 150,
          report_names = NULL,
          dat = obj$env$data,
          ... ){

  # unload stuff
  start_time = Sys.time()
  init.map = obj$env$map
  MLE = obj$env$parList()
  # dat is passed as argument in case there's some special types involved, e.g., for explicit units

  # set defaults
  report = obj$report()
  if(is.null(report_names)){
    report_names = names(report)
  }else{
    report_names = intersect( report_names, names(report) )
  }
  if(length(report_names)==0) stop("Check `report_names`")

  # process fixed effects
  FE <- obj$env$last.par[-obj$env$random]
  FE <- split( unname(FE), names(FE) )

  # re-make map
  #map <- lapply(names(FE), function(x) factor(FE[[x]]*NA))
  map <- lapply( names(FE), function(x) factor(MLE[[x]]*NA) )
  names(map) <- names(FE)
  map <- c(map, init.map[!(names(init.map)%in%names(map))])

  # check for issues
  ok = ( sapply(MLE[names(map)],length) == sapply(map,length) )
  if(!all(ok)) stop("Check problem with `mcmc_detransformation`")

  # Rebuild with original FE mapped off and RE as FE
  obj_random <- MakeADFun( data = dat,
                           parameters = MLE,
                           map = map )

  # run STAN
  stan_out <- tmbstan::tmbstan( obj_random,
                                chains = 1,
                                warmup = warmup,
                                iter = iter,
                                refresh = -1,
                                ... )
  stan_samples = as.array(stan_out)[,1,1:length(obj_random$par)]

  # Compile output
  report_full = objmle$report(samples[1,])[report_names]
  for( i in 2:nrow(samples)){
    report_sample = objmle$report(samples[i,])[report_names]
    for( z in seq_along(report_full) ){
      if(is.vector(report_sample[[z]])){
        report_full[[z]] = cbind( report_full[[z]], report_sample[[z]] )
      }else{
        report_full[[z]] = abind::abind( report_full[[z]], report_sample[[z]], along=length(dim(report_sample[[z]]))+1 )
      }
    }
  }
  run_time = Sys.time() - start_time

  #
  out = list( "stan_out" = stan_out,
              "report_full" = report_full,
              "run_time" = run_time )
  return( out )
}

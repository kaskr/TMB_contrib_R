

#' Extract fixed effects
#'
#' \code{extract_fixed} extracts the best previous value of fixed effects, in a way that works for both mixed and fixed effect models
#'
#' @param obj, The compiled object
#'
#' @return A vector of fixed-effect estimates

extract_fixed = function( obj ){
  if( length(obj$env$random)==0 ){
    Return = obj$env$last.par.best
  }else{
    Return = obj$env$last.par.best[-c(obj$env$random)]
  }
  return( Return )
}


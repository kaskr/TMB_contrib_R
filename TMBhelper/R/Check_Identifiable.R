
#' Check for identifiability of fixed effects
#'
#' \code{Check_Identifiable} calculates the matrix of second-derivatives of the marginal likelihood w.r.t. fixed effects, to see if any linear combinations are unidentifiable
#'
#' @param obj, The compiled object
#'
#' @return A tagged list of the hessian and the message

#' @export
Check_Identifiable = function( obj ){
  # Finite-different hessian
  List = NULL
  List[["Hess"]] = optimHess( par=obj$env$last.par.best, fn=obj$fn, gr=obj$gr )

  # Check eigendecomposition
  List[["Eigen"]] = eigen( List[["Hess"]] )
  List[["WhichBad"]] = which( List[["Eigen"]]$values < sqrt(.Machine$double.eps) )

  # Check for parameters
  RowMax = apply( List[["Eigen"]]$vectors[,List[["WhichBad"]]], MARGIN=1, FUN=function(vec){max(abs(vec))} )
  List[["BadParams"]] = data.frame("Param"=names(obj$par), "MLE"=obj$env$last.par.best, ifelse(RowMax>0.1, "Bad","OK"))

  # Message
  if( length(List[["WhichBad"]])==0 ){
    print( "All are identifiable" )
  }else{
    print( List[["BadParams"]] )
  }

  # Return
  return( invisible(List) )
}

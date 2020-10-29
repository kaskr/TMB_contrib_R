
#' Check for identifiability of fixed effects
#'
#' \code{check_estimability} calculates the matrix of second-derivatives of the marginal likelihood
#' w.r.t. fixed effects, to see if any linear combinations are not estimable (i.e. cannot be
#' uniquely estimated conditional upon model structure and available data, e.g., resulting
#' in a likelihood ridge and singular, non-invertable Hessian matrix)
#'
#' @param obj The compiled object
#' @param h optional argument containing pre-computed Hessian matrix
#'
#' @return A tagged list of the hessian and the message

#' @export
check_estimability = function( obj, h ){

  # Extract fixed effects
  ParHat = TMBhelper:::extract_fixed( obj )

  # Check for problems
  Gr = obj$gr( ParHat )
  if( any(Gr>0.01) ) stop("Some gradients are high, please improve optimization and only then use `Check_Identifiable`")

  # Finite-different hessian
  List = NULL
  if(missing(h)){
    List[["Hess"]] = optimHess( par=ParHat, fn=obj$fn, gr=obj$gr )
  }else{
    List[["Hess"]] = h
  }

  # Check eigendecomposition
  List[["Eigen"]] = eigen( List[["Hess"]] )
  List[["WhichBad"]] = which( List[["Eigen"]]$values < sqrt(.Machine$double.eps) )

  # Check result
  if( length(List[["WhichBad"]])==0 ){
    # print message
    message( "All parameters are estimable" )
  }else{
    # Check for parameters
    RowMax = apply( List[["Eigen"]]$vectors[,List[["WhichBad"]],drop=FALSE], MARGIN=1, FUN=function(vec){max(abs(vec))} )
    List[["BadParams"]] = data.frame("Param"=names(obj$par), "MLE"=ParHat, "Param_check"=ifelse(RowMax>0.1, "Bad","OK"))
    # print message
    print( List[["BadParams"]] )
  }

  # Return
  return( invisible(List) )
}

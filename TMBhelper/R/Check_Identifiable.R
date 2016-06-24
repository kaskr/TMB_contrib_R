
#' Check for identifiability of fixed effects
#'
#' \code{Check_Identifiable} calculates the matrix of second-derivatives of the marginal likelihood w.r.t. fixed effects, to see if any linear combinations are unidentifiable
#'
#' @param obj, The compiled object
#'
#' @return A tagged list of the hessian and the message

Check_Identifiable = function( obj ){
  # Finite-different hessian
  Hess = optimHess( par=obj$env$last.par.best, fn=obj$fn, gr=obj$gr )

  # Check eigendecomposition
  Eigen = eigen( Hess )
  Which = which( Eigen$values < sqrt(.Machine$double.eps) )

  # Return
  List = NULL
  if( length(Which)==0 ){
    List[["Message"]] = "All are identifiable"
  }else{
    List[["Message"]] = data.frame( names(obj$par), Eigen$vectors[,Which] )
  }
  List[["Hess"]] = H

  # Message and return
  print( List[["Message"]] )
  return( List )
}

#' Copy of fit_tmb
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{?fit_tmb} to see list of arguments and usage
#' @export
Optimize = function( ... ){
  .Deprecated( new="fit_tmb" )
  fit_tmb( ... )
}

#' Copy of check_estimability
#'
#' Included for continuity when using old scripts
#'
#' Please use \code{\link{check_estimability}} to see list of arguments and usage
#' @export
Check_Identifiable = function( ... ){
  .Deprecated( new="check_estimability" )
  check_estimability( ... )
}

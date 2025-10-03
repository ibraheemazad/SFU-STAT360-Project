#' Print method for the MARS object
#'
#' Prints the function call of MARS object, then each coefficient of selected variables.
#'
#'
#' @param x A `mars` object
#' @param ... Additional arguments
#'
#' @family method
#'
#' @return Function call and coefficients of the MARS object
#' @export
#'
#' @examples
#' mm <- mars(y~., data=mars::marstestdata)
#' print(mm)
print.mars <- function(x, ...){
  print("Call:")
  print(x$call)
  print("Coefficients:")
  print(x$coefficients)
}

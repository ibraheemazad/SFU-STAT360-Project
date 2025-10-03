#' Plot method for the mars object
#'
#' Plots the fitted basis functions that depend on
#' explanatory variable(main effects) or two explanatory variables
#' (two-way interactions).
#'
#'
#' @param x A `mars` object
#' @param ... Additional arguments for plotting
#'
#' @family method
#'
#' @return Four diagnose plots:
#' \itemize{
#' \item Residuals vs Fitted value
#' \item Normal Q-Q plot
#' \item Squared standardized residuals vs. Fitted value
#' \item Cook's distance.
#' }
#' @export
#'
#' @import stats
#' @import grDevices
#' @import graphics
#'
#' @examples
#' mm <- mars(y~., data=mars::marstestdata)
#' plot(mm)
plot.mars <- function(x, ...){
  data <- eval(x$call$data)
  tt <- terms(x$formula, data=data)
  tt <- delete.response(tt)
  mf <- model.frame(tt, data)
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)[,-1] # remove intercept
  Bf <- x$Bfuncs

  # opar <- par(mfrow=c(2,2))
  # par(mfrow=c(2,2))
  # par(mar=c(1,1,1,1))

  # plot 1: fitted vs residual
  yh <- x$fitted.values[order(x$fitted.values)]
  r <- x$residuals[order(yh)]
  ylim <- range(r, na.rm=TRUE)
  ylim <- extendrange(r = ylim, f = 0.08)
  plot(yh, r, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs. Fitted",
       ylim = ylim, type = "b", ...)

  abline(h = 0, lty = 3, col = "gray")

  # plot 2: qq plot
  rs <- scale(r)
  qqnorm(rs); qqline(rs, col="red")

  # plot 3: sqrt. std. residual vs fitted
  sqrtabsr <- sqrt(abs(rs))
  ylim <- range(sqrtabsr, na.rm=TRUE)
  ylim <- extendrange(r = ylim, f = 0.08)
  plot(yh, sqrtabsr, xlab = "Fitted values", ylab = "Squared Std. Residuals", main = "Squared Std. Residuals vs. Fitted",
       ylim = ylim, type = "b", ...)

  abline(h = 0, lty = 3, col = "gray")

  # plot 4: cook distance
  cook <- cooks.distance(x)
  ymx <- max(cook) * 1.075
  plot(cook, type = "h", ylim = c(0, ymx), main = "Cook's Distance",
       xlab = "Obs. number", ylab = "Cook's distance", ...)
  # on.exit(par(opar),add=TRUE)

}



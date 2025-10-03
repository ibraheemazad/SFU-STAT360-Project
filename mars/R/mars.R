#' Fit a MARS model
#'
#' @description Build a regression model using the techniques in Friedman's papers "Multivariate Adaptive Regression Splines (MARS)".
#'
#' @details The algorithm will search for, and discover, nonlinearities in the data that help maximize predictive accuracy.
#' Multivariate adaptive regression splines (MARS), an algorithm that essentially creates a piecewise linear model
#' which provides an intuitive stepping block into nonlinearity after grasping the concept of linear regression
#' and other intrinsically linear models.
#'
#' @author Tauseef Kashtwari, Promit Chowdhury, Ibraheem Azad
#'
#' @usage mars(formula, data, control)
#'
#' @param formula An R formula
#' @param data A data frame containing the data for the model
#' @param control An object of class `mars.control`
#'
#' @family method
#'
#' @return An object of class `mars` that includes the final regression and a description
#' of the basis functions. There are plot, predict, summary and print methods for mars object.
#'
#' @export
#'
#' @seealso
#'  [mars.control()] for constructing the control object \cr
#'  [plot.mars()] for plotting the results \cr
#'  [predict.mars()] for predictions \cr
#'  [summary.mars()] for summarizing mars objects \cr
#'  [print.mars()] for printing mars objects
#'
#'
#' @examples
#' mm <- mars(y~., data=mars::marstestdata)
#' @import stats
#' @import ISLR
#' @import grDevices
#' @import graphics
#'
#' @references Jerome H. Friedman. “Multivariate Adaptive Regression Splines.”
#' Ann. Statist. 19(1) 1 - 67, March, 1991. https://doi.org/10.1214/aos/1176347963.
#'
mars <- function(formula,data,control=NULL) {
  cc <- match.call() # save the call
  mf <- model.frame(formula,data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)[,-1, drop=FALSE]
  x_names <- colnames(x)
  if(is.null(control)) control <- mars.control()
  control <- validate_mars.control(control)
  fwd <- fwd_stepwise(y,x,control) # Note change from earlier implementation
  bwd <- bwd_stepwise(fwd,control)
  final_model <- lm(y~.-1, data=data.frame(y=y,bwd$B))
  output <- c(list(call=cc, formula=formula, y=y, B=bwd$B, Bfuncs=bwd$Bfuncs,
                   x_names=x_names), final_model)
  class(output) <- c("mars", "lm")
  return(output)
}

#' Mars Control Object
#' @description Constructor for `mars.control` objects
#'
#' This function constructs a `mars.control` object that specifies
#' parameters used in the model fitting procedure
#'
#' @param Mmax Maximum number of basis functions. Should be an even integer. Default is 2.
#' @param d The coefficient in the penalty term of the generalized cross validation measure. Default is 3.
#' @param trace Should we print status information about the fitting? Default is `FALSE`.
#'
#' @return a `mars.control` object
#' @export
#'
#' @examples
#' mc <- mars.control(Mmax=10)
mars.control <- function(Mmax=2, d=3, trace=FALSE){
  if(Mmax < 2) {
    Mmax = 2
    warning("Mmax is smaller than 2, and therefore, set to be 2")
  }
  R <- list(Mmax=Mmax, d=d, trace=trace)
  R <- validate_mars.control(R)
  new_mars.control(R)
}


new_mars.control <- function(R){
  structure(R,class="mars.control")
}

validate_mars.control <- function(R){
  if(R$Mmax %% 2 == 1 | (!is.numeric(R$Mmax)) | R$Mmax < 2) stop("Mmax is not an even integer >= 2")
  if(!is.numeric(R$d)) stop("d is not a numeric value")
  if(!is.logical(R$trace)) stop("trace is not a logical value")
  R
}


# Forward stepwise function takes no arguments and return an empty list
fwd_stepwise <- function(y,x,control=mars.control()){

  Mmax <- control$Mmax

  N <- length(y) # sample size
  n <- ncol(x) # number of predictors
  B <- init_B(N,Mmax) # N-by-(Mmax+1) matrix
  # splits <- data.frame(m=rep(NA,Mmax),v=rep(NA,Mmax),t=rep(NA,Mmax))
  Bfuncs <- vector("list", Mmax+1)

  for(i in 1:(Mmax/2)) {
    M <- 2*i - 1
    lof_best <- Inf
    for(m in 1:M) { # choose a basis function to split
      if (is.null(Bfuncs[[m]][,"v"])) {
        svars <- 1:n
      } else {
        svars <- setdiff(1:n, Bfuncs[[m]][,"v"])
      }
      for(v in svars){ # select a variable to split on
        tt <- split_points(x[,v],B[,m])
        for(t in tt) {
          Bnew <- data.frame(B[,1:M],
                             Btem1=B[,m]*h(x[,v], +1, t),
                             Btem2=B[,m]*h(x[,v], -1, t))
          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(y~.,gdat, control)
          if(lof < lof_best) {
            lof_best <- lof
            split_best <- c(m=m, v=v, t=t)
          } # end if
        } # end loop over splits
      } # end loop over variables
    } # end loop over basis functions to split
    m <- split_best["m"]; v <- split_best["v"]; t <- split_best["t"]
    # B[,M+1] <- cbind(B[,m]*h(x[,v], -1, t))
    # B[,M+2] <- cbind(B[,m]*h(x[,v], +1, t))
    B[,M+1] <- B[,m]*h(x[,v], -1, t)
    B[,M+2] <- B[,m]*h(x[,v], +1, t)
    Bfuncs[[M+1]] <- rbind(Bfuncs[[m]], c(s=-1, v, t))
    Bfuncs[[M+2]] <- rbind(Bfuncs[[m]], c(s=+1, v, t))
  } # end loop over i

  colnames(B) <- paste0("B", (0:(ncol(B)-1)))
  return(list(y=y, B=B, Bfuncs=Bfuncs))
}

# Backward stepwise takes the output of fwd_stepwise() and for now return an empty list
bwd_stepwise <- function(fwd, control){
  y <- fwd$y; B <- fwd$B; Bfuncs <- fwd$Bfuncs
  Mmax <- ncol(B)-1
  J_star <- 2:(Mmax+1) # All non-intercept basis functions
  K_star <- J_star
  data <- data.frame(y=y, B)
  lof_star <- LOF(y~.-1, data, control)
  for (M in (Mmax+1):2){
    L <- K_star
    b <- Inf
    for (m in L){
      K <- setdiff(L,m)
      gdat <- data.frame(y=y,B=B[,K])
      lof <- LOF(y~.,gdat,control)
      if (b > lof){
        b <- lof
        K_star <- K
      }
    }
    if (b < lof_star){
      lof_star <- b
      J_star <- K_star
    }
  }
  J_star <- c(1,J_star)
  return(list(y=y, B=B[,J_star], Bfuncs=Bfuncs[J_star]))
}

init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}
split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}

# hinge function
h <- function(x,s,t){
  return(pmax(0, s*(x-t)))
}

LOF <- function(formula, data, control){
  ff <- lm(formula,data)
  RSS <- sum(residuals(ff)^2)
  C_tilda <- sum(hatvalues(ff)) + control$d * (length(coefficients(ff))-1)
  return(RSS * (nrow(data))/(nrow(data)-C_tilda)^2)
}

#' Predict method for the MARS object
#'
#' Provides predictions from the results of the MARS object's basis functions.
#'
#' @param object A `mars` object
#' @param newdata New data **optional*
#' @param ... Additional arguments
#'
#' @family method
#'
#' @return If `newdata` is missing, fitted values are returned;
#' otherwise, predicted values on the new data are returned
#' @export
#'
#' @examples
#' mm <- mars(y~x1+x2, data=mars::marstestdata)
#' pred <- predict(mm, newdata=data.frame(x1=rnorm(100), x2=rnorm(100)))
predict.mars <- function(object,newdata, ...) {
  if(missing(newdata) || is.null(newdata)) {
    # extract fitted model
    B <- as.matrix(object$B)
  } else {
    # predict on new data
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1, drop=F] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  pred <- rep(NA, nrow(B))
  for (n in 1:nrow(B)){
    pred[n] <- sum(beta * B[n,])
  }
  pred
  # drop(B %*% beta)
}

make_B <- function(data, Bfuncs){
  new_B <- matrix(1, nrow=nrow(data), ncol=length(Bfuncs))
  # new_B[,1] <- rep(1, nrow(data)) # intercept
  for (i in 2:length(Bfuncs)){
    new_B[,i] <- create_basis_function(data, Bfuncs[[i]])
  }
  new_B
}
create_basis_function <- function(data, Child_Bfuncs){
  B <- rep(1, nrow(data))  # initialize
  for (j in 1:nrow(Child_Bfuncs)){
    # "ind" selects the points that are within the child region
    ind <- sign(data[,Child_Bfuncs[j,2]]-Child_Bfuncs[j,3]) == Child_Bfuncs[j,1]
    B <- B*Child_Bfuncs[j,1]*(data[,Child_Bfuncs[j,2]]-Child_Bfuncs[j,3])*ind
  }
  as.matrix(B, ncol=1)
}

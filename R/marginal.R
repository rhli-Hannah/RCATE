#' Marginal treatment effect plot from machine learning algorithms.
#'
#' \code{marginal.rcate.ml} Returns the variable importance level from "rcate.ml" model.
#'
#' @param object "rcate.ml" object.
#' @param variable.col the column number of interested variable. Default is 1.
#' @param ... other.
#' @rdname marginal.rcate.ml
#' @export
marginal.rcate.ml <- function(object, variable.name=NULL,...){
  model <- object$model
  algorithm <- object$algorithm
  x <- object$x
  x.mean <- object$param$x.mean
  x.sd <- object$param$x.sd

  if (is.null(variable.name)) {
    variable.name <- colnames(object$param$x.scaled)[1]
  }

  x.select <- seq(min(x[,variable.name]),max(x[,variable.name]),length.out = 200)
  x.cond.scaled <- (x.select-x.mean[variable.name])/x.sd[variable.name]
  x.cond <- data.frame(matrix(0,ncol = ncol(x),nrow = 200))
  colnames(x.cond) <- colnames(object$param$x.scaled)
  x.cond[,variable.name] <- x.cond.scaled
  if (algorithm =='GBM') {
    pred <- predict(model,x.cond)
    graphics::plot(x.select,pred,type = 'l',xlab=variable.name)
  } else if (algorithm == 'NN') {
    pred <- rowMeans(predict(model,as.matrix(x.cond)))
    graphics::plot(x.select,pred,type = 'l',xlab=variable.name)
  }
}


#' Marginal treatment effect plot from additive model.
#'
#' \code{marginal.rcate.am} Returns the variable importance level from "rcate.am" model.
#'
#' @param object "rcate.am" object.
#' @param variable.col the column number of interested variable. Default is 1.
#' @param ... other.
#' @rdname marginal.rcate.am
#' @export
marginal.rcate.am <- function(object, variable.col=1,...){
  algorithm <- object$algorithm
  model <- object$model
  colnum <- object$colnum

  x <- object$x
  x.colmean <- colMeans(object$x)
  x.cond <- matrix(rep(colMeans(object$x), times=200),
                   nrow = 200,ncol = ncol(object$x), byrow=TRUE)
  x.cond[,variable.col] <- seq(min(object$x[,variable.col]),max(object$x[,variable.col]),length.out = 200)
  x.select <- seq(min(object$x[,variable.col]),max(object$x[,variable.col]),length.out = 200)

  center.x <- apply(x, 2, function(x) (x - object$mean.x)/object$sd.x)
  center.xval <- apply(x.cond, 2, function(x) (x - object$mean.x)/object$sd.x)

  if (algorithm == "SAM") {
    center.xval <- center.xval[, colnum]
    center.x <- center.x[, colnum]

    if (is.vector(center.xval) & is.vector(center.x)) {
      center.xval <- matrix(center.xval, ncol = 1)
      center.x <- matrix(center.x, ncol = 1)
    }

    Btilde.val <- B_R(center.xval, center.x, object$lambda.smooth, object$nknots, colnum)
    predict <- cbind(rep(1, nrow(Btilde.val)), Btilde.val) %*% object$coef
  }

  graphics::plot(x.select,predict,type = 'l')
}

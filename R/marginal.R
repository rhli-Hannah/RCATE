#' Marginal treatment effect plot.
#'
#' \code{marginal.rcate} Returns the variable importance level from "rcate.ml" model.
#'
#' @param object "rcate.ml" or "rcate.am" object.
#' @param variable.name the name of interested variable. Default is the name of the
#' first continuous variable.
#' @param ... other.
#' @rdname marginal.rcate
#' @export
marginal.rcate <- function(object, variable.name=NULL,...){
  model <- object$model
  algorithm <- object$algorithm
  x <- object$x
  x.mean <- object$param$x.mean
  x.sd <- object$param$x.sd

  if (is.null(variable.name)) {
    variable.name <- colnames(object$param$x.scaled)[1]
  }

  x.select <- seq(min(x[,variable.name]),max(x[,variable.name]),length.out = 100)
  x.cond.scaled <- (x.select-x.mean[variable.name])/x.sd[variable.name]
  x.cond <- data.frame(matrix(0,ncol = ncol(x),nrow = 100))
  colnames(x.cond) <- colnames(object$param$x.scaled)
  x.cond[,variable.name] <- x.cond.scaled

  if (algorithm =='GBM') {
    pred <- predict(model,x.cond, ntrees=object$n.trees.gbm)
    graphics::plot(x.select,pred,type = 'l',xlab=variable.name)
  } else if (algorithm == 'NN') {
    pred <- rowMeans(predict(model,as.matrix(x.cond)))
    graphics::plot(x.select,pred,type = 'l',xlab=variable.name)
  } else if (algorithm == 'SAM') {
    center.xval <- x.cond[, object$colnum]
    center.xval <- center.xval+matrix(stats::rnorm(100*ncol(center.xval),0,0.001),nrow = 100)
    center.x <- object$param$x.scaled[, object$colnum]

    if (is.vector(center.xval) & is.vector(center.x)) {
      center.xval <- matrix(center.xval, ncol = 1)
      center.x <- matrix(center.x, ncol = 1)
    }

    Btilde.val <- B_R(center.xval, center.x, object$lambda.smooth, object$nknots, object$colnum)
    predict <- as.matrix(cbind(rep(1, nrow(Btilde.val)),
                               Btilde.val, x.cond[, !names(x.cond) %in% object$param$name.num])) %*% object$coef
    graphics::plot(x.select,predict,type = 'l')
  }
}

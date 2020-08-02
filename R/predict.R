#' Predict treatment effect.
#'
#' \code{predict} returns the robust treatment effect estimation result.
#'
#' @param object "RCATE" object.
#' @param x matrix or a data frame of predictors.
#' @return a list of components
#' \itemize{
#'  \item predict - the robust estimation result of CATE.
#'  \item x - matrix of predictors.
#'  \item algorithm - fitting algorithm.
#'  \item model - "RCATE" object.
#'  \item method - estimation method.
#'  }
#' @examples
#' n <- 1000; p <- 3
#' X <- matrix(runif(n*p,-3,3),nrow=n,ncol=p)
#' tau = sin(X[,1])
#' p = 1/(1+exp(-X[,1]))
#' d = rbinom(n,1,p)
#' t = 2*d-1
#' y = 1+tau*t/2 + rnorm(n,0,0.5)
#' fit <- rcate.ml(X,y,d)
#' y_pred <- predict.rcate.ml(fit,X)$predict
predict.rcate.ml <- function(object, x) {
  algorithm <- object$algorithm
    model <- object$model

    if (algorithm == "GBM") {
      predict <- predict(model, data.frame(x), n.trees = object$n.trees.gbm)
    } else if (algorithm == "NN") {
      predict <- model %>% predict(x)
    }

  return(list(predict = predict, x = x, algorithm = object$algorithm,
              model = model, method = object$method))
}

predict.rcate.am <- function(object, x) {
    algorithm <- object$algorithm
    model <- object$model
    colnum <- object$colnum

    center.x <- apply(object$x, 2, function(x) (x - object$mean.x)/object$sd.x)
    center.xval <- apply(x, 2, function(x) (x - object$mean.x)/object$sd.x)

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


  return(list(predict = predict, x = x, algorithm = object$algorithm,
              model = model, method = object$method))
}

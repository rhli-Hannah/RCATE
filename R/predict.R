#' Predict treatment effect from machine learning algorithms.
#'
#' \code{predict.rcate.ml} Returns the predicted treatment effect from "rcate.ml" model.
#'
#' @param object "rcate.ml" object.
#' @param x matrix or a data frame of predictors.
#' @param ... other.
#' @return a list of components
#' \itemize{
#'  \item predict - a robust estimation result of CATE.
#'  \item x - matrix of predictors.
#'  \item algorithm - fitting algorithm.
#'  \item model - "rcate.ml" object.
#'  \item method - estimation method.
#'  }
#' @rdname predict.rcate.ml
#' @export
predict.rcate.ml <- function(object, x, ...) {
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

#' Prediction of treatment effect from an \eqn{L_1} b-spline additive regression penalized model
#'
#' \code{predict.rcate.am} Returns predicted treatment effect from "rcate.am" model.
#'
#' @param object "rcate.am" object.
#' @param x matrix or a data frame of predictors.
#' @param ... other.
#' @return a list of components
#' \itemize{
#'  \item predict - a robust estimation result of CATE.
#'  \item x - matrix of predictors.
#'  \item algorithm - fitting algorithm.
#'  \item model - "rcate.am" object.
#'  \item method - estimation method.
#'  }
#' @rdname predict.rcate.am
#' @export
predict.rcate.am <- function(object, x,...) {
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



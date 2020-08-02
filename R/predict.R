#' Predict treatment effect from machine learning algorithms.
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
#' n <- 1000; p <- 10
#' X <- matrix(rnorm(n*p,0,1),nrow=n,ncol=p)
#' tau = 6*sin(2*X[,1])+3*(X[,2]+3)*X[,3]+9*tanh(0.5*X[,4])+3*X[,5]*(2*I(X[,4]<1)-1)
#' p = 1/(1+exp(-X[,1]+X[,2]))
#' d = rbinom(n,1,p)
#' t = 2*d-1
#' y = 100+4*X[,1]+X[,2]-3*X[,3]+tau*t/2 + rnorm(n,0,1)
#' x_val = matrix(rnorm(200*10,0,1),nrow=200,ncol=10)
#' tau_val = 6*sin(2*x_val[,1])+3*(x_val[,2]+3)*x_val[,3]+9*tanh(0.5*x_val[,4])+3*x_val[,5]*(2*I(x_val[,4]<1)-1)
#' # Use MCM-EA transformation and GBM to estimate CATE
#' fit <- rcate.ml(X,y,d,method='RL')
#' y_pred <- predict.rcate.ml(fit,x_val)$predict
#' plot(tau_val,y_pred);abline(0,1)
#' @export
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

#' Predict treatment effect from additive model.
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
#' n <- 1000; p <- 10
#' X <- matrix(rnorm(n*p,0,1),nrow=n,ncol=p)
#' tau = 6*sin(2*X[,1])+3*(X[,2])+X[,3]+9*tanh(0.5*X[,4])+3*X[,5]
#' p = 1/(1+exp(-X[,1]+X[,2]))
#' d = rbinom(n,1,p)
#' t = 2*d-1
#' y = 100+4*X[,1]+X[,2]-3*X[,3]+tau*t/2 + rnorm(n,0,1)
#' x_val = matrix(rnorm(200*10,0,1),nrow=200,ncol=10)
#' tau_val = 6*sin(2*x_val[,1])+3*(x_val[,2])+x_val[,3]+9*tanh(0.5*x_val[,4])+3*x_val[,5]
#'
#' fit <- rcate.am(X,y,d)
#' y_pred <- predict.rcate.am(fit,x_val)$pred
#' plot(tau_val,y_pred);abline(0,1)
#' @export
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

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
    param <- object$param
    x.mean <- param$x.mean
    x.sd <- param$x.sd
    name.num <- param$name.num

    x.num <- dplyr::select_if(x, is.numeric)
    scaled <- NULL
    for (i in 1:ncol(x.num)) {
      scaled1 <- (x.num[,i]-x.mean[i])/x.sd[i]
      scaled <- cbind(scaled,scaled1)
    }
    x.num.scaled <- scaled
    x.other <- data.frame(x[ , -which(names(x) %in% name.num)])
    if (ncol(x.other)==0) {
      x.scaled <- x.num.scaled
    } else {
      x.other <- apply(x.other, 2, function(x) as.numeric(as.character(x)))
      x.scaled <- cbind(x.num.scaled,x.other)
    }
    colnames(x.scaled) <- colnames(object$param$x.scaled)

    if (algorithm == "GBM") {
      predict <- predict(model, data.frame(x.scaled), n.trees = object$n.trees.gbm)
    } else if (algorithm == "NN") {
      predict <- rowMeans(predict(model,as.matrix(x.scaled)))
      model <- NULL
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
    name.num <- object$param$name.num

    center.x <- object$param$x.num.scaled
    center.xval <- NULL
    for (i in 1:ncol(center.x)) {
      xval1 <- (x[,name.num][,i]-object$param$x.mean[i])/object$param$x.sd[i]
      center.xval <- cbind(center.xval,xval1)
    }
    center.xval <- center.xval+matrix(stats::rnorm(200*ncol(center.xval),0,0.001),nrow = 200)

    if (algorithm == "SAM") {
      if (is.vector(center.xval) & is.vector(center.x)) {
        center.xval <- matrix(center.xval, ncol = 1)
        center.x <- matrix(center.x, ncol = 1)
      }


      Btilde.val <- B_R(center.xval, center.x, object$lambda.smooth, object$nknots, seq(1:ncol(center.xval)))
      mat.val <- as.matrix(cbind(rep(1, nrow(Btilde.val)), Btilde.val,
                                 x[, !names(x) %in% name.num]))
      predict <- mat.val %*% object$coef
    }


  return(list(predict = predict, x = x, algorithm = object$algorithm,
              model = model, method = object$method))
}


#' Predict treatment effect from robust random forests.
#'
#' \code{predict.rcate.rf} predicts treatment effect from robust random forests.
#'
#' @param object a robust random forests.
#' @param newdata dataframe contains covariates.
#' @param ... other.
#' @return a list of components
#' \itemize{
#'  \item pred - prediction of newdata.
#'  \item newdata - a test data frame.
#'  }
#' @rdname predict.rcate.rf
#' @export
predict.rcate.rf <- function(object, newdata,...) {
  tree_info <- object$tree
  pred.mat <- matrix(NA, nrow = length(tree_info), ncol = nrow(newdata))
  nodepreds <- object$nodepreds
  param <- object$param
  x.mean <- param$x.mean
  x.sd <- param$x.sd
  name.num <- param$name.num

  x.num <- dplyr::select_if(newdata, is.numeric)
  scaled <- NULL
  for (i in 1:ncol(x.num)) {
    scaled1 <- (x.num[,i]-x.mean[i])/x.sd[i]
    scaled <- cbind(scaled,scaled1)
  }
  x.num.scaled <- scaled
  x.other <- data.frame(newdata[ , -which(names(newdata) %in% name.num)])
  if (ncol(x.other)==0) {
    x.scaled <- x.num.scaled
  } else {
    x.other <- apply(x.other, 2, function(x) as.numeric(as.character(x)))
    x.scaled <- cbind(x.num.scaled,x.other)
  }
  colnames(x.scaled) <- colnames(object$param$x.scaled)

  for (j in 1:length(tree_info)) {
    leafs <- tree_info[[j]][tree_info[[j]]$TERMINAL=='LEAF',]
    nodepred <- nodepreds[[j]]
    predicted <- c()
    for (i in seq_len(nrow(leafs))) {
      # extract index
      ind <- as.numeric(rownames(subset(as.data.frame(x.scaled),
                                        eval(parse(text = leafs[i, "FILTER"])))))
      # estimator is the median y value of the leaf
      predicted[ind] <- nodepred[i]
    }
    pred.mat[j,] <- predicted
  }
  pred <- colMeans(pred.mat)
  return(list(pred = pred, newdata=newdata))
}




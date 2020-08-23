#' Marginal treatment effect plot.
#'
#' \code{marginal.rcate} Returns the marginal treatment effect plot.
#'
#' @param object "rcate.ml", "rcate.rf" or "rcate.am" object.
#' @param variable.name the name of interested variable. Default is the name of the
#' first continuous variable.
#' @param ... other.
#' @rdname marginal.rcate
#' @export
marginal.rcate <- function(object, variable.name=NULL,...){
  x.mean <- object$param$x.mean
  x.sd <- object$param$x.sd

  if (is.null(variable.name)) {
    variable.name <- colnames(object$param$x.scaled)[1]
  }

  algorithm <- object$algorithm

  if (algorithm =='GBM') {
    model <- object$model
    x <- object$x

    x.select <- seq(min(x[,variable.name]),max(x[,variable.name]),length.out = 100)
    x.cond.scaled <- (x.select-x.mean[variable.name])/x.sd[variable.name]
    x.cond <- data.frame(matrix(0,ncol = ncol(x),nrow = 100))
    colnames(x.cond) <- colnames(object$param$x.scaled)
    x.cond[,variable.name] <- x.cond.scaled

    pred <- predict(model,x.cond, ntrees=object$n.trees.gbm)
    graphics::plot(x.select,pred,type = 'l',xlab=variable.name, ylab='Predict')
  } else if (algorithm == 'NN') {
    model <- object$model
    x <- object$x

    x.select <- seq(min(x[,variable.name]),max(x[,variable.name]),length.out = 100)
    x.cond.scaled <- (x.select-x.mean[variable.name])/x.sd[variable.name]
    x.cond <- data.frame(matrix(0,ncol = ncol(x),nrow = 100))
    colnames(x.cond) <- colnames(object$param$x.scaled)
    x.cond[,variable.name] <- x.cond.scaled

    pred <- rowMeans(predict(model,as.matrix(x.cond)))
    graphics::plot(x.select,pred,type = 'l',xlab=variable.name, ylab='Predict')
  } else if (algorithm == 'SAM') {
    model <- object$model
    x <- object$x

    x.select <- seq(min(x[,variable.name]),max(x[,variable.name]),length.out = 100)
    x.cond.scaled <- (x.select-x.mean[variable.name])/x.sd[variable.name]
    x.cond <- data.frame(matrix(0,ncol = ncol(x),nrow = 100))
    colnames(x.cond) <- colnames(object$param$x.scaled)
    x.cond[,variable.name] <- x.cond.scaled

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
    graphics::plot(x.select,predict,type = 'l',xlab=variable.name, ylab='Predict')
  } else if (algorithm == 'RF') {
    tree_info <- object$tree
    pred.mat <- matrix(NA, nrow = length(tree_info), ncol = 100)
    nodepreds <- object$nodepreds
    param <- object$param
    x <- object$x

    x.select <- seq(min(x[,variable.name]),max(x[,variable.name]),length.out = 100)
    x.cond.scaled <- (x.select-x.mean[variable.name])/x.sd[variable.name]
    x.cond <- data.frame(matrix(0,ncol = ncol(x),nrow = 100))
    colnames(x.cond) <- colnames(object$param$x.scaled)
    x.cond[,variable.name] <- x.cond.scaled

    for (j in 1:length(tree_info)) {
      leafs <- tree_info[[j]][tree_info[[j]]$TERMINAL=='LEAF',]
      nodepred <- nodepreds[[j]]
      predicted <- c()
      for (i in seq_len(nrow(leafs))) {
        # extract index
        ind <- as.numeric(rownames(subset(as.data.frame(x.cond),
                                          eval(parse(text = leafs[i, "FILTER"])))))
        # estimator is the median y value of the leaf
        predicted[ind] <- nodepred[i]
      }
      pred.mat[j,] <- predicted
    }
    pred <- colMeans(pred.mat)
    graphics::plot(x.select,pred,type = 'l',xlab=variable.name, ylab='Predict')
  }
}

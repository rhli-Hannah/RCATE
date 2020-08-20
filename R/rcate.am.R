# Functions for smoothness-sparsity penalty
# library(rqPen);library(splines)

Btilde.fun <- function(x, lambda2, knots2) {
  B <- splines::bs(x, degree = 3, knots = knots2, Boundary.knots = c(-4, 4))
  D <- diff(diag(ncol(B)), differences = 2)
  Omega <- crossprod(D)
  n <- length(x)
  M <- 1/n * crossprod(B) + lambda2 * Omega
  R <- chol(M)
  R1 <- solve(R)
  Btilde.mat <- B %*% R1
  return(Btilde.mat)
}

B_R <- function(x2, x3, lambda2br, knots.op1,colnum) {
  Btilde0 <- NULL
  result <- knot(x3, knots.op1)
  for (i in 1:length(colnum)) {
    Btilde1 <- Btilde.fun(x2[, i], lambda2br, result[, i])
    Btilde0 <- cbind(Btilde0, Btilde1)
  }
  return(Btilde0)
}

# Genrate knots
knot <- function(x, knots.op) {
  knot.mat = apply(x, 2, function(y)
    stats::quantile(y, probs = utils::head(seq(0, 1, length.out = as.numeric(knots.op + 2)), -1)[-1]))
  return(knot.mat)
}


#' Robust estimation of treatment effect using additive b-spline model.
#'
#' \code{rcate.am} fit additive model for robust treatment effect estimation.
#'
#' Fits a \eqn{L_1} additive regression model with the group SCAD penalty.
#'
#' @param x matrix or a data frame of predictors.
#' @param y vector of response values.
#' @param d vector of binary treatment assignment (0 or 1).
#' @param method character string of CATE estimation method: "MCMEA" - modified co-variate
#' method with efficiency augmentation or "RL" - R-learning.
#' @param  NIS logical of using non-parametric independent screening or not. The default is TRUE
#' when \eqn{p < n/log(n)}.
#' @param  lambda.smooth scalar represent the smoothness penalty lambda2. The default is 2.
#' @param  nlambda number of lambda1s for cross-validation. The default is 30.
#' @param  nfolds number of folds for cross-validation. The default is 5.
#' @param n.trees.p tuning parameter the number of trees used for estimating propensity score
#' with GBM. the default value is 40000.
#' @param shrinkage.p tuning parameter the shrinkage level for estimating propensity score
#' with GBM. the default value is 0.005.
#' @param n.minobsinnode.p tuning parameter the minimum node size for estimating propensity score
#' with GBM. the default value is 10.
#' @param interaction.depth.p tuning parameter the number of interactions for estimating propensity
#' score with GBM. the default value is 1.
#' @param cv.p tuning parameter the number of folds in cross-validation for estimating propensity
#' score with GBM. the default value is 2.
#' @param n.trees.mu scalar or vector of the number of trees for estimating mean function with GBM.
#' The default is (1:50)*50.
#' @param shrinkage.mu tuning parameter the shrinkage level for estimating mean function
#' with GBM. the default value is 0.01.
#' @param n.minobsinnode.mu tuning parameter the minimum node size for estimating mean function
#' with GBM. the default value is 10.
#' @param interaction.depth.mu tuning parameter the number of interactions for estimating mean
#' function with GBM. the default value is 5.
#' @param cv.mu tuning parameter the folds for cross-validation for estimating mean function with
#' GBM. The default value is 5.
#' @return a list of components
#' \itemize{
#'  \item model - the robust estimation model of CATE.
#'  \item method - estimation method.
#'  \item algorithm - fitting algorithm.
#'  \item lambda.smooth
#'  \item fitted.values - vector of fitted values.
#'  \item x - matrix of predictors.
#'  \item y - vector of response values.
#'  \item d - vector of treatment assignment.
#'  \item y.tr - transformed outcome.
#'  \item w.tr - transformed weight.
#'  \item coef -coefficients.
#'  \item colnum - column number.
#'  \item nknots - number of knots of cubic spline.
#'  \item param - required parameters for utility functions.
#'  }
#' @examples
#' n <- 1000; p <- 2; set.seed(2223)
#' X <- as.data.frame(matrix(runif(n*p,-3,3),nrow=n,ncol=p))
#' tau = 3*X[,1]-2*X[,2]
#' p = 1/(1+exp(-X[,1]+X[,2]))
#' d = rbinom(n,1,p)
#' t = 2*d-1
#' y = 100+tau*t/2 + rnorm(n,0,1); set.seed(2223)
#' x_val = as.data.frame(matrix(rnorm(200*2,0,1),nrow=200,ncol=2))
#' tau_val = 3*x_val[,1]-2*x_val[,2]
#'
#' fit <- rcate.am(X,y,d,lambda.smooth = 4, method = 'RL')
#' y_pred <- predict(fit,x_val)$pred
#' plot(tau_val,y_pred);abline(0,1)
#'
#' @export
rcate.am <- function(x, y, d, method = "MCMEA", nknots = NA,
                     lambda.smooth = 1, nlambda = 30, nfolds = 5, n.trees.p = 40000,
                     shrinkage.p = 0.005, n.minobsinnode.p = 10,
                     interaction.depth.p = 1, cv.p = 2, n.trees.mu = c(1:50) * 50,
                     shrinkage.mu = 0.01,
                     n.minobsinnode.mu = 5, interaction.depth.mu = 5, cv.mu = 5) {
  options(warn=-1)

  # Calculate T=2D-1
  t <- 2 * d - 1
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }
  # Number of rows and columns of X
  row <- nrow(x)
  col <- ncol(x)

  # Calculate number of knots
  if (is.na(nknots)) {
     nknots <- floor(sqrt(row)/2)
  }

  # Standardize X
  x.num <- dplyr::select_if(x, is.numeric)
  x.mean <- apply(x.num, 2, mean)
  x.sd <- apply(x.num, 2, sd)
  x.num.scaled <- scale(x.num)
  name.num <- colnames(x.num)
  x.other <- data.frame(x[ , -which(names(x) %in% name.num)])
  if (ncol(x.other)==0) {
    x.scaled <- x.num.scaled
    group2 <- rep(1:ncol(x.num.scaled), each = nknots + 3)
    group <- group2
  } else {
    x.other <- apply(x.other, 2, function(x) as.numeric(as.character(x)))
    x.scaled <- cbind(x.num.scaled,x.other)
    group2 <- rep(1:ncol(x.num.scaled), each = nknots + 3)
    group1 <- seq(1:ncol(x.other))
    group <- c(group2,group1+max(group2))
  }

  colnum <- seq(1:ncol(x.num.scaled))
  x.tr <- x.scaled


  # If X has only one dimension, transform it into a matrix with one column
  if (is.vector(x.tr)) {
    x.tr <- matrix(x.tr, ncol = 1)
  }

  # Estimate mu0(x), mu1(x) and p(x)
  data.p <- data.frame(cbind(factor(d), x))
  colnames(data.p)[1] <- c("d")
  data.p$d <- as.factor(data.p$d)
  gbmGrid.p <- expand.grid(interaction.depth = interaction.depth.p,
                           n.trees = n.trees.p, shrinkage = shrinkage.p,
                           n.minobsinnode = n.minobsinnode.p)
  gbmFit.p <- caret::train(d ~ ., data = data.p, method = "gbm",
                           verbose = FALSE, trControl = caret::trainControl(method = "cv", number = cv.p),
                           tuneGrid = gbmGrid.p)
  pscore.hat <- caret::predict.train(gbmFit.p, newdata = data.p, type = "prob")[, 2]

  data00 <- data.frame(cbind(y, x, d))
  colnames(data00)[c(1,ncol(data00))] <- c("y", "d")
  # data020 <- data00[data00$d == 0, ]
  # data021 <- data00[data00$d == 1, ]

  gbmGrid.mu <- expand.grid(interaction.depth = interaction.depth.mu,
                            n.trees = n.trees.mu, shrinkage = shrinkage.mu,
                            n.minobsinnode = n.minobsinnode.mu)
  gbmFit.mu <- caret::train(y ~ ., data = data00[, -ncol(data00)],
                            method = "gbm", verbose = FALSE, trControl = caret::trainControl(method = "cv",
                            number = cv.mu), tuneGrid = gbmGrid.mu, metric = "MAE")
  # gbmFit.mu1 <- caret::train(y ~ ., data = data021[, -ncol(data021)],
  #                            method = "gbm", verbose = FALSE, trControl = caret::trainControl(method = "cv",
  #                            number = cv.mu), tuneGrid = gbmGrid.mu, metric = "MAE")
  # gbmFit.mu0 <- caret::train(y ~ ., data = data020[, -ncol(data020)],
  #                            method = "gbm", verbose = FALSE, trControl = caret::trainControl(method = "cv",
  #                            number = cv.mu), tuneGrid = gbmGrid.mu, metric = "MAE")
  #
  # mu0 <- caret::predict.train(gbmFit.mu0, newdata = data00)
  # mu1 <- caret::predict.train(gbmFit.mu1, newdata = data00)
  mu.ea <- caret::predict.train(gbmFit.mu, newdata = data00)

  # Do transformation
  if (method == "MCMEA") {
    w.tr <- 1/(t * pscore.hat + (1 - t)/2)
    y.tr <- (y - mu.ea) * w.tr
    wmat.tr <- matrix(0, row, row)
    diag(wmat.tr) <- w.tr * t/2

    Btilde <- B_R(x.num.scaled, x.num.scaled, lambda.smooth, nknots, colnum)
    x.tr1 <- wmat.tr %*% as.matrix(cbind(rep(1, row), Btilde, x.other))
  } else if (method == "RL") {
    w.tr <- 1
    y.tr <- (y - mu.ea) * w.tr
    wmat.tr <- matrix(0, row, row)
    diag(wmat.tr) <- w.tr * (t - 2 * pscore.hat + 1)/2

    Btilde <- B_R(x.num.scaled, x.num.scaled, lambda.smooth, nknots, colnum)
    x.tr1 <- wmat.tr %*% as.matrix(cbind(rep(1, row), Btilde, x.other))
  }

  # Fit the weighted LAD model with group SCAD
  model <- rqPen::cv.rq.group.pen(x = x.tr1, y = y.tr, groups = as.factor(c(max(group)+1, group)),
                                  tau = 0.5, intercept = FALSE, penalty = 'SCAD',
                                  nfolds = nfolds, criteria = "BIC", nlambda = nlambda)
  # Get coefficients and fitted value
  coef <- coef(model)
  fitted.values <- as.matrix(cbind(rep(1, nrow(Btilde)), Btilde, x.other)) %*% coef

  result <- list(model = model, method = method, algorithm = "SAM",
                 lambda.smooth = lambda.smooth, fitted.values = fitted.values,
                 x = x, y = y, d = d, y.tr = y.tr, w.tr = w.tr,
                 coef = coef, colnum = colnum, nknots = nknots,
                 param = list(x.scaled=x.scaled, name.num=name.num,
                              x.mean=x.mean, x.sd=x.sd, x.num.scaled = x.num.scaled))
  class(result) <- "rcate.am"
  return(result)
}

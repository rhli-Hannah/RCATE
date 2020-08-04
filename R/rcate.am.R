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

# Adaptive INIS
adaptINIS <- function(x, y, testdata = NULL, lambda.pen.list = NULL,
                      folds = NULL, quant = NULL, kfold = NULL, knots = NULL, eps0 = 1e-06,
                      DOISIS = TRUE, maxloop = 10, trace = FALSE, detailed = FALSE) {
  t0 = proc.time()[1]
  cat("starting adaptINIS, NIS algorithm, adatively choose number of variables\n")
  n <- nrow(x)
  p <- ncol(x)


  # if(is.null(nsis)) nsis=min(floor(n/log(n)),p-1)
  if (is.null(knots)) {
    knots = ceiling(n^0.2)
  }
  if (is.null(folds)) {
    temp = sample(1:n, n, replace = FALSE)
    if (is.null(kfold))
      kfold = 5
    for (i in 1:kfold) {
      folds[[i]] = setdiff(1:n, temp[seq(i, n, kfold)])
    }
  }
  if (is.null(quant)) {
    quant = 1
  }
  df0 <- knots + 1

  xbs = matrix(0, n, df0 * p)


  for (i in 1:p) {
    xbs[, (i - 1) * (df0) + (1:df0)] = splines::ns(x[, i], df = df0)
  }


  tempresi <- rep(0, p)

  curloop = 1
  for (i in 1:p) {
    tempfit <- stats::lm.fit(x = cbind(1, xbs[, (i - 1) * df0 + 1:df0]), y = y)
    tempresi[i] <- sum(tempfit$residuals^2)
  }


  used.list <- tempresi

  used.sort <- sort(used.list, method = "sh", index = TRUE, decreasing = FALSE)
  initRANKorder <- used.sort$ix


  mindex <- sample(1:n)
  mresi = NULL
  for (i in 1:p) {
    tempfit <- stats::lm.fit(x = cbind(1, xbs[, (i - 1) * df0 + 1:df0]), y = y[mindex])
    mresi[i] <- sum(tempfit$residuals^2)
  }
  resi.thres = stats::quantile(mresi, 1 - quant)
  nsis <- max(min(sum(used.list < resi.thres), floor(n/df0/3)), 2)


  SISind <- sort(initRANKorder[1:nsis])
  if (!DOISIS)
    return(list(initRANKorder = initRANKorder, SISind = SISind, nsis = nsis))

  cat("loop ", curloop, "...SISind ", SISind, "\n")
  pick.ind = initRANKorder[1:nsis]

  return(pick.ind)
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
#' @param method character string of CATE estimation method: "MCMEA" - modified covariate
#' method with efficiency augmentation, "RL" - R-learning, or "DR" - doubly robust method.
#' @param  NIS logical of using non-parametric independice screening or not. The default is TRUE
#' when \eqn{p < n/log(n)}.
#' @param  nknots number of knots for cubic spline. The default is \eqn{floor(sqrt(n)/2)}.
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
#'  \item mean.x - column means of predictors.
#'  \item sd.x -column standard deviation of predictors.
#'  \item coef -coefficients.
#'  \item colnum - column number.
#'  \item nknots - number of knots of cubic spline.
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
#' y_pred <- predict(fit,x_val)$pred
#' plot(tau_val,y_pred);abline(0,1)
#' @export
rcate.am <- function(x, y, d, method = "MCMEA", NIS = TRUE, nknots = NA,
                     lambda.smooth = 2, nlambda = 30, nfolds = 5, n.trees.p = 40000,
                     shrinkage.p = 0.005, n.minobsinnode.p = 10,
                     interaction.depth.p = 1, cv.p = 2, n.trees.mu = c(1:50) * 50,
                     shrinkage.mu = 0.01,
                     n.minobsinnode.mu = 5, interaction.depth.mu = 5, cv.mu = 5) {
  # Calculate T=2D-1
  t <- 2 * d - 1
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }
  # Number of rows and columns of X
  row <- nrow(x)
  col <- ncol(x)

  # Standardize X
  mean.x <- apply(x, 2, mean)
  sd.x <- apply(x, 2, sd)
  center.x <- (x - mean.x)/sd.x

  # Calculate number of knots
  if (is.na(nknots)) {
    nknots <- floor(sqrt(row)/2)
  }

  # NIS
  if (col < floor(row/log(row)) | NIS == FALSE) {
    colnum <- 1:col
  } else {
    colnum_d <- adaptINIS(x, d)
    colnum_y1 <- adaptINIS(x[d == 1, ], y[d == 1])
    colnum_y0 <- adaptINIS(x[d == 0, ], y[d == 0])
    colnum <- sort(unique(union(union(colnum_d, colnum_y0), colnum_y1)))
  }
  group <- rep(1:length(colnum), each = nknots + 3)
  x.tr <- center.x[, colnum]

  # If X has only one dimension, transform it into a matrix with one column
  if (is.vector(x.tr)) {
    x.tr <- matrix(x.tr, ncol = 1)
  }

  # Estimate mu0(x), mu1(x) and p(x)
  data.p <- data.frame(cbind(factor(d), x))
  colnames(data.p) <- c("d", paste0("X", 1:ncol(x)))
  data.p$d <- as.factor(data.p$d)
  gbmGrid.p <- expand.grid(interaction.depth = interaction.depth.p,
                           n.trees = n.trees.p, shrinkage = shrinkage.p,
                           n.minobsinnode = n.minobsinnode.p)
  gbmFit.p <- caret::train(d ~ ., data = data.p, method = "gbm",
                           verbose = FALSE, trControl = caret::trainControl(method = "cv", number = cv.p),
                           tuneGrid = gbmGrid.p)
  pscore.hat <- caret::predict.train(gbmFit.p, newdata = data.p, type = "prob")[, 2]

  data00 <- data.frame(cbind(y, x, d))
  colnames(data00) <- c("y", paste0("X", 1:ncol(x)), "d")
  data020 <- data00[data00$d == 0, ]
  data021 <- data00[data00$d == 1, ]

  gbmGrid.mu <- expand.grid(interaction.depth = interaction.depth.mu,
                            n.trees = n.trees.mu, shrinkage = shrinkage.mu,
                            n.minobsinnode = n.minobsinnode.mu)
  gbmFit.mu <- caret::train(y ~ ., data = data00[, -ncol(data00)],
                            method = "gbm", verbose = FALSE, trControl = caret::trainControl(method = "cv",
                            number = cv.mu), tuneGrid = gbmGrid.mu, metric = "MAE")
  gbmFit.mu1 <- caret::train(y ~ ., data = data021[, -ncol(data021)],
                             method = "gbm", verbose = FALSE, trControl = caret::trainControl(method = "cv",
                             number = cv.mu), tuneGrid = gbmGrid.mu, metric = "MAE")
  gbmFit.mu0 <- caret::train(y ~ ., data = data020[, -ncol(data020)],
                             method = "gbm", verbose = FALSE, trControl = caret::trainControl(method = "cv",
                             number = cv.mu), tuneGrid = gbmGrid.mu, metric = "MAE")

  mu0 <- caret::predict.train(gbmFit.mu0, newdata = data00)
  mu1 <- caret::predict.train(gbmFit.mu1, newdata = data00)
  mu.ea <- caret::predict.train(gbmFit.mu, newdata = data00)

  # Do transformation
  if (method == "MCMEA") {
    w.tr <- 1/(t * pscore.hat + (1 - t)/2)
    y.tr <- (y - mu.ea) * w.tr
    wmat.tr <- matrix(0, row, row)
    diag(wmat.tr) <- w.tr * t/2

    Btilde <- B_R(x.tr, x.tr, lambda.smooth, nknots, colnum)
    x.tr1 <- wmat.tr %*% cbind(rep(1, row), Btilde)
  } else if (method == "RL") {
    w.tr <- 1
    y.tr <- (y_tr - mu.ea) * w.tr
    wmat.tr <- matrix(0, row, row)
    diag(wmat.tr) <- w.tr * (t - 2 * pscore.hat + 1)/2

    Btilde <- B_R(x.tr, x.tr, lambda2, nknots, colnum)
    x.tr1 <- wmat.tr %*% cbind(rep(1, row), Btilde)
  } else if (method == "DR") {
    y.tr <- (t - 2 * pscore.hat + 1) * (y)/(2 * pscore.hat * (1 - pscore.hat)) +
      (pscore.hat - d)/pscore.hat * mu1 + (pscore.hat - d)/(1 - pscore.hat) * mu0
    w.tr <- abs((t - 2 * pscore.hat + 1)/(2 * pscore.hat * (1 - pscore.hat)))
  }

  # Fit the weighted LAD model with group SCAD
  model <- rqPen::cv.rq.group.pen(x = x.tr1, y = y.tr, groups = as.factor(c(max(group)+1, group)),
                                  tau = 0.5, intercept = FALSE, penalty = 'SCAD',
                                  nfolds = nfolds, criteria = "BIC", nlambda = nlambda)
  # Get coefficients and fitted value
  coef <- coef(model)
  fitted.values <- cbind(rep(1, nrow(Btilde)), Btilde) %*% coef

  result <- list(model = model, method = method, algorithm = "SAM",
                 lambda.smooth = lambda.smooth, fitted.values = fitted.values,
                 x = x, y = y, d = d, mean.x = mean.x, sd.x = sd.x,
                 coef = coef, colnum = colnum, nknots = nknots)
  class(result) <- "rcate.am"
  return(result)
}

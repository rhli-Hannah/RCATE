#library(MASS);library(gbm);library(caret);library(plyr);library(dplyr);library(randomForest);library(rmutil);library(keras);library(tensorflow)


#' Robust estimation of treatment effect.
#'
#' \code{rcate.ml} fit ML algorithm for robust treatment effect estimation.
#'
#' Fit a GBM or NN to estimate treatment effect estimation that robust to outliers.
#'
#' @param x matrix or a data frame of predictors.
#' @param y vector of response values.
#' @param d vector of binary treatment assignment (0 or 1).
#' @param  algorithm character string of algorithm: "GBM" - gradient boosting machine or
#' "NN" - neural network. The random forests is available in rcate.rf.
#' @param method character string of CATE estimation method: "MCMEA" - modified co-variate
#' method with efficiency augmentation, "RL" - R-learning, or "DR" - doubly robust method.
#' @param n.trees.gbm tuning parameter the number of trees used in GBM for estimating treatment
#' effect function if algorithm="GBM". The default is 1000.
#' @param interaction.depth.gbm tuning parameter the number of interactions for estimating treatment
#' effect function if algorithm="GBM". The default value is 2.
#' @param n.cells.nn vector of the number of neurals in each hidden layer if algorithm='NN'.
#' The default is two layers with each layer the half size of previous layer.
#' @param dropout.nn vector of the dropout rate of each hidden layer if algorithm='NN'.
#' The default is no dropout.
#' @param epochs.nn scalar of the number of epochs for neural network if algorithm='NN'.
#' The defualt is 100.
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
#'  \item fitted.values - vector of fitted values.
#'  \item x - matrix of predictors.
#'  \item y - vector of response values.
#'  \item d - vector of treatment assignment.
#'  \item y.tr - vector of transformed outcome.
#'  \item w.tr - vector of transformed weight.
#'  \item n.trees.gbm - number of trees for estimating treatment effect function if algorithm='GBM'.
#'  \item history - model fitting history.
#'  \item importance - variable importance level.
#'  \item param - required parameters for utility functions.
#'  }
#' @examples
#' n <- 1000; p <- 3; set.seed(2223)
#' X <- as.data.frame(matrix(runif(n*p,-3,3),nrow=n,ncol=p))
#' tau = 6*sin(2*X[,1])+3*(X[,2]+3)*X[,3]
#' p = 1/(1+exp(-X[,1]+X[,2]))
#' d = rbinom(n,1,p)
#' t = 2*d-1
#' y = 100+4*X[,1]+X[,2]-3*X[,3]+tau*t/2 + rnorm(n,0,1); set.seed(2223)
#' x_val = as.data.frame(matrix(rnorm(200*3,0,1),nrow=200,ncol=3))
#' tau_val = 6*sin(2*x_val[,1])+3*(x_val[,2]+3)*x_val[,3]
#' # Use R-learning method and GBM to estimate CATE
#' fit <- rcate.ml(X,y,d,method='RL',algorithm='GBM')
#' y_pred <- predict(fit,x_val)$predict
#' plot(tau_val,y_pred);abline(0,1)
#'
#' # Use doubly robust method and neural network to estimate CATE
#' fit <- rcate.ml(X,y,d,method='DR',algorithm='NN',dropout.nn=c(0),n.cells.nn=c(3,3))
#' @importFrom stats predict
#' @export
rcate.ml <- function(x, y, d, method = "MCMEA", algorithm = "GBM",
                  n.trees.p = 40000, shrinkage.p = 0.005, n.minobsinnode.p = 10,
                  interaction.depth.p = 1, cv.p = 5, n.trees.mu = c(1:50) * 50,
                  shrinkage.mu = 0.01, n.minobsinnode.mu = 5,
                  interaction.depth.mu = 5, cv.mu = 5, n.trees.gbm = 1000,
                  interaction.depth.gbm = 2, n.cells.nn = NA, dropout.nn = NA,
                  epochs.nn = 100) {
  options(warn=-1)

  # Calculate T=2D-1
  t <- 2 * d - 1

  # Generate the number of cells used in NN
  if (is.na(n.cells.nn)) {
    n.cells.nn <- c(ncol(x), ceiling(max(ncol(x) * 0.5,1)))
  }

  # Generate the dropout rate of NN
  if (is.na(dropout.nn)) {
    dropout.nn <- c(0.5)
  }

  x.num <- dplyr::select_if(x, is.numeric)
  x.mean <- apply(x.num, 2, mean)
  x.sd <- apply(x.num, 2, sd)
  x.num.scaled <- scale(x.num)
  name.num <- colnames(x.num)
  x.other <- data.frame(x[ , -which(names(x) %in% name.num)])
  if (ncol(x.other)==0) {
    x.scaled <- x.num.scaled
  } else {
    x.other <- apply(x.other, 2, function(x) as.numeric(as.character(x)))
    x.scaled <- cbind(x.num.scaled,x.other)
  }

  # Estimate mu0(x), mu1(x) and p(x)
  data.p <- data.frame(cbind(d, x.scaled))
  #colnames(data.p) <- c("d", paste0("X", 1:ncol(x)))
  data.p$d <- as.factor(data.p$d)
  gbmGrid.p <- expand.grid(interaction.depth = interaction.depth.p,
                           n.trees = n.trees.p, shrinkage = shrinkage.p,
                           n.minobsinnode = n.minobsinnode.p)
  gbmFit.p <- caret::train(d ~ ., data = data.p, method = "gbm",
                           verbose = FALSE,
                           trControl = caret::trainControl(method = "cv", number = cv.p),
                           tuneGrid = gbmGrid.p)
  pscore.hat <- caret::predict.train(gbmFit.p, newdata = data.p, type = "prob")[, 2]

  data00 <- data.frame(cbind(y, x.scaled, d))
  colnames(data00)[c(1,ncol(data00))] <- c("y", "d")
  data020 <- data00[data00$d == 0, ]
  data021 <- data00[data00$d == 1, ]

  gbmGrid.mu <- expand.grid(interaction.depth = interaction.depth.mu,
                            n.trees = n.trees.mu, shrinkage = shrinkage.mu,
                            n.minobsinnode = n.minobsinnode.mu)
  gbmFit.mu <- caret::train(y ~ ., data = data00[, -ncol(data00)],
                            method = "gbm", verbose = FALSE,
                            trControl = caret::trainControl(method = "cv",
                            number = cv.mu), tuneGrid = gbmGrid.mu, metric = "MAE")
  gbmFit.mu1 <- caret::train(y ~ ., data = data021[, -ncol(data021)],
                             method = "gbm", verbose = FALSE,
                             trControl = caret::trainControl(method = "cv",
                             number = cv.mu), tuneGrid = gbmGrid.mu, metric = "MAE")
  gbmFit.mu0 <- caret::train(y ~ ., data = data020[, -ncol(data020)],
                             method = "gbm", verbose = FALSE,
                             trControl = caret::trainControl(method = "cv",
                             number = cv.mu), tuneGrid = gbmGrid.mu, metric = "MAE")

  mu0 <- caret::predict.train(gbmFit.mu0, newdata = data00)
  mu1 <- caret::predict.train(gbmFit.mu1, newdata = data00)
  mu.ea <- caret::predict.train(gbmFit.mu, newdata = data00)

  # Do transformation
  if (method == "MCMEA") {
    y.tr <- 2 * t * (y - mu.ea)
    w.tr <- 1/2 * (t * pscore.hat + (1 - t)/2)
  } else if (method == "RL") {
    y.tr <- 2 * (y - mu.ea)/(t - 2 * pscore.hat + 1)
    w.tr <- abs(t - 2 * pscore.hat + 1)/2
  } else if (method == "DR") {
    y.tr <- (t - 2 * pscore.hat + 1) * (y)/(2 * pscore.hat * (1 - pscore.hat)) +
      (pscore.hat - d)/pscore.hat * mu1 + (pscore.hat - d)/(1 - pscore.hat) * mu0
    w.tr <- abs((t - 2 * pscore.hat + 1)/(2 * pscore.hat * (1 - pscore.hat)))
  }

  # Fit GBM
  if (algorithm == "GBM") {
    data.gbm <- data.frame(cbind(y.tr, x.scaled))
    colnames(data.gbm)[1] <- c("y.tr")
    model <- gbm::gbm(y.tr ~ ., data = data.gbm, distribution = "laplace",
                      weights = w.tr, n.trees = n.trees.gbm,
                      interaction.depth = interaction.depth.gbm)
    df.x <- data.frame(x.scaled)
    #colnames(df.x) <- paste0("X", 1:ncol(x))
    fitted.values <- gbm::predict.gbm(model, df.x, n.trees = n.trees.gbm)

    history <- NULL

    importance <- gbm::relative.influence(model,n.trees = n.trees.gbm)
  } else if (algorithm == "NN") {
    x = as.matrix(x.scaled)
    y = as.matrix(y.tr)

    keras_model_simple_mlp <- function(use_bn = FALSE, use_dp = FALSE,
                                       name = NULL) {

      # define and return a custom model
      keras::keras_model_custom(name = name, function(self) {

        # create layers we'll need for the call (this code executes once)
        self$dense1 <- keras::layer_dense(units = n.cells.nn[1], activation = "relu")
        self$dense1.1 <- keras::layer_dense(units = n.cells.nn[1], activation = "relu")
        self$dense1.5 <- keras::layer_dropout(rate = dropout.nn[1])
        self$dense2 <- keras::layer_dense(units = n.cells.nn[2:length(n.cells.nn)], activation = "relu")
        self$dense2.5 <- keras::layer_dropout(rate = dropout.nn[1])
        self$dense3 <- keras::layer_dense(units = 1, activation = "linear")

        if (use_dp)
          self$dp <- keras::layer_dropout(rate = 0.5)
        if (use_bn)
          self$bn <- keras::layer_batch_normalization(axis = -1)

        # implement call (this code executes during training & inference)
        function(inputs, mask = NULL) {
          x <- self$dense1(inputs)
          if (use_dp)
            x <- self$dp(x)
          if (use_bn)
            x <- self$bn(x)
          self$dense3(x)
        }
      })
    }

    model <- keras_model_simple_mlp()
    keras::compile(model, loss = "mae", optimizer = "adam")
    history <- keras::fit(model,x, y, epochs = epochs.nn, verbose = 0, sample_weight = w.tr)

    fitted.values <- rowMeans(predict(model,x))

    pred_wrapper <- function(object, newdata) {
      as.vector(predict(object, x = as.matrix(newdata)))
    }
    p1 <- vip::vip(
      object = model,                     # fitted model
      method = "permute",                 # permutation-based VI scores
      num_features = ncol(x),       # default only plots top 10 features
      pred_wrapper = pred_wrapper,            # user-defined prediction function
      train = as.data.frame(x) ,    # training data
      target = y,                   # response values used for training
      metric = "rsquared",                # evaluation metric
      # progress = "text"                 # request a text-based progress bar
    )
    importance <- p1$data
  }

  result <- list(model = model, method = method, algorithm = algorithm,
                 fitted.values = fitted.values, x = x, y = y, d = d,
                 y.tr = y.tr, w.tr = w.tr,
                 n.trees.gbm = n.trees.gbm, history = history, importance = importance,
                 param = list(x.scaled=x.scaled, name.num=name.num,
                              x.mean=x.mean, x.sd=x.sd))
  class(result) <- "rcate.ml"
  return(result)
}

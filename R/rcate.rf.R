#Spliting rule.
wmed_var <- function(x, y, weight) {
    splits <- sort(unique(x))
    wmed <- c()
    for (i in seq_along(splits)) {
        sp <- splits[i]
        wmed[i] <- sum(abs(y[x < sp] - stats::median(y[x < sp])) * weight[x < sp]) +
            sum(abs(y[x >= sp] - stats::median(y[x >= sp])) * weight[x >= sp])
    }
    split_at <- splits[which.min(wmed)]
    return(c(wmed = min(wmed), split = split_at))
}

#' Robust regression tree.
#'
#' \code{reg_tree} is robust regression tree algorithm based on weighted median
#' splitting rule.
#'
#' @param formula an object of class "formula".
#' @param data a training data frame.
#' @param minsize minimal leaf size of tree. Default is 5.
#' @param newdata an optional test data frame. If NULL, the newdata=data.
#' @param weights an optional vector of weights.
#' @return a list of components
#' \itemize{
#'  \item tree - tree information.
#'  \item fit - estimation method.
#'  \item formula - an object of class "formula".
#'  \item importance - vector of importance level.
#'  \item data - a training data frame.
#'  \item pred - prediction of newdata.
#'  \item newdata - a test data frame.
#'  }
#' @examples
#' n <- 1000; p <- 3
#' X <- matrix(runif(n*p,-3,3),nrow=n,ncol=p)
#' y = 1+sin(X[,1]) + rnorm(n,0,0.5)
#' df <- data.frame(y,X)
#' tree <- reg_tree_imp(y~X1+X2+X3,data=df,minsize=3,newdata=df,weights=rep(1,1000))
#' y_pred <- tree$pred
#' plot(y,y_pred)
#' @export
reg_tree_imp <- function(formula, data, minsize, newdata, weights) {

    # coerce to data.frame
    data <- as.data.frame(data)
    row.names(data) <- seq(1:nrow(data))
    newdata <- as.data.frame(newdata)

    # handle formula
    formula <- stats::terms.formula(formula)

    # get the design matrix
    X <- stats::model.matrix(formula, data)

    # extract target
    y <- data[, as.character(formula)[2]]

    weight <- weights

    # initialize while loop
    do_splits <- TRUE

    # create output data.frame with splitting rules and observations
    tree_info <- data.frame(NODE = 1, NOBS = nrow(data), FILTER = NA, TERMINAL = "SPLIT",
                            IMP_GINI = NA, SPLIT = NA, stringsAsFactors = FALSE)

    # keep splitting until there are only leafs left
    while(do_splits) {

        # which parents have to be splitted
        to_calculate <- which(tree_info$TERMINAL == "SPLIT")

        for (j in to_calculate) {

            # handle root node
            if (!is.na(tree_info[j, "FILTER"])) {
                # subset data according to the filter
                this_data <- subset(data, eval(parse(text = tree_info[j, "FILTER"])))
                # get the design matrix
                X <- stats::model.matrix(formula, this_data)
                weight <- weights[as.numeric(row.names(X))]
                # print(X)
                # print(weight)
            } else {
                this_data <- data
            }

            # estimate splitting criteria
            splitting <- apply(X,  MARGIN = 2, FUN = wmed_var, y = this_data[, all.vars(formula)[1]], weight=weight)

            # get the min SSE
            tmp_splitter <- which.min(splitting[1,])

            # define maxnode
            mn <- max(tree_info$NODE)

            # paste filter rules
            current_filter <- c(paste(names(tmp_splitter), ">=",
                                      splitting[2,tmp_splitter]),
                                paste(names(tmp_splitter), "<",
                                      splitting[2,tmp_splitter]))

            # Error handling! check if the splitting rule has already been invoked
            split_here  <- !sapply(current_filter,
                                   FUN = function(x,y) any(grepl(x, x = y)),
                                   y = tree_info$FILTER)

            # append the splitting rules
            if (!is.na(tree_info[j, "FILTER"])) {
                current_filter  <- paste(tree_info[j, "FILTER"],
                                         current_filter, sep = " & ")
            }

            # calculate metrics within the children
            metr <- lapply(current_filter,
                           FUN = function(i, x, data, formula) {
                               df <- subset(x = x, subset = eval(parse(text = i)))
                               nobs <- nrow(df)
                               w <- nobs/nrow(data)
                               y <- df[, all.vars(formula)[1]]
                               imp <- mean(abs(y - stats::median(y, na.rm = TRUE)))
                               return(c(nobs, w*imp))
                           },
                           x = this_data, data = data, formula = formula)

            # extract relevant information
            current_nobs <- sapply(metr, function(x) x[[1]])
            imp_sum_child <- sum(sapply(metr, function(x) x[[2]]))
            current_y <- this_data[, all.vars(formula)[1]]
            imp_parent <- nrow(this_data)/nrow(data) * mean(abs(current_y-stats::median(current_y)))
            imp_gini <- imp_parent - imp_sum_child

            # insufficient minsize for split
            if (any(current_nobs <= minsize)) {
                split_here <- rep(FALSE, 2)
            }

            # create children data frame
            children <- data.frame(NODE = c(mn+1, mn+2),
                                   NOBS = current_nobs,
                                   FILTER = current_filter,
                                   TERMINAL = rep("SPLIT", 2),
                                   IMP_GINI = NA,
                                   SPLIT = NA,
                                   row.names = NULL)[split_here,]

            # overwrite state of current node, add gini importance and split variable
            tree_info[j, "TERMINAL"] <- ifelse(all(!split_here), "LEAF", "PARENT")
            tree_info[j, "IMP_GINI"] <- imp_gini
            if (tree_info[j, "TERMINAL"] == "PARENT") {
                tree_info[j, "SPLIT"] <- names(tmp_splitter)
            }

            # bind everything
            tree_info <- rbind(tree_info, children)

            # check if there are any open splits left
            do_splits <- !all(tree_info$TERMINAL != "SPLIT")
        } # end for
    } # end while

    # calculate fitted values
    leafs <- tree_info[tree_info$TERMINAL == "LEAF", ]
    fitted <- c()
    nodepred <- rep(NA,nrow(leafs))
    for (i in seq_len(nrow(leafs))) {
        # extract index
        ind <- as.numeric(rownames(subset(data, eval(parse(text = leafs[i, "FILTER"])))))
        # estimator is the median y value of the leaf
        fitted[ind] <- stats::median(y[ind])
        nodepred[i] <- stats::median(y[ind])
    }

    # calculate predicted values
    predicted <- c()
    for (i in seq_len(nrow(leafs))) {
        # extract index
        ind <- as.numeric(rownames(subset(newdata, eval(parse(text = leafs[i, "FILTER"])))))
        # estimator is the median y value of the leaf
        predicted[ind] <- nodepred[i]
    }

    # calculate feature importance
    imp <- tree_info[, c("SPLIT", "IMP_GINI")]

    if (!all(is.na(imp$SPLIT))) {
        imp <- stats::aggregate(IMP_GINI ~ SPLIT, FUN = function(x, all) sum(x, na.rm = T)/sum(all, na.rm = T),
                         data = imp, all = imp$IMP_GINI)
    }

    # rename to importance
    names(imp) <- c("FEATURES", "IMPORTANCE")
    imp <- imp[order(imp$IMPORTANCE, decreasing = TRUE),]

    # return everything
    return(list(tree = tree_info, fit = fitted, formula = formula,
                importance = imp, data = data, pred=predicted, newdata=newdata))
}

#' Robust random forests.
#'
#' \code{reg_rf} is robust random forest algorithm based on weighted median
#' splitting rule.
#'
#' @param formula an object of class "formula".
#' @param n_trees number of trees. Default is 50.
#' @param feature_frac fraction of features used in each split. Default is 1/2.
#' @param data a training data frame.
#' @param minnodes minimal leaf size of tree. Default is 5.
#' @param newdata an optional test data frame. If NULL, the newdata=data.
#' @param weights an optional vector of weights.
#' @return a list of components
#' \itemize{
#'  \item fit - estimation method.
#'  \item pred - prediction of newdata.
#'  }
#' @examples
#' n <- 1000; p <- 3
#' X <- matrix(runif(n*p,-3,3),nrow=n,ncol=p)
#' y = 1+sin(X[,1]) + rnorm(n,0,0.5)
#' df <- data.frame(y,X)
#' RF <- reg_rf(y~X1+X2+X3,data=df,newdata=df,weights=rep(1,1000))
#' y_pred <- RF$pred
#' plot(y,y_pred)
#' @export
reg_rf <- function(formula, n_trees=50, feature_frac=1/2, data, newdata, weights, minnodes=5) {

    # define function to sprout a single tree
    sprout_tree <- function(formula, feature_frac, data, newdata, weights) {
        # extract features
        features <- all.vars(formula)[-1]

        # extract target
        target <- all.vars(formula)[1]

        # bag the data
        # - randomly sample the data with replacement (duplicate are possible)
        train.id <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
        train <- data[train.id,]
        w.train <- weights[train.id]


        # randomly sample features
        # - only fit the regression tree with feature_frac * 100 % of the features
        features_sample <- sample(features,
                                  size = ceiling(length(features) * feature_frac),
                                  replace = FALSE)

        # create new formula
        formula_new <-
            stats::as.formula(paste0(target, " ~ -1 + ", paste0(features_sample,
                                                         collapse =  " + ")))

        # fit the regression tree
        tree <- reg_tree_imp(formula = formula_new,
                             data = train, newdata=newdata,
                             minsize = minnodes, weights = w.train)

        # save the fit and the importance
        return(list(tree$fit, tree$importance, tree$pred))
    }

    # apply the rf_tree function n_trees times with plyr::raply
    # - track the progress with a progress bar
    trees <- plyr::raply(
        n_trees,
        sprout_tree(
            formula = formula,
            feature_frac = feature_frac,
            data = data,
            newdata = newdata,
            weights = weights
        ),
        .progress = "text"
    )

    # extract fit
    fits <- do.call("cbind", trees[, 1])
    preds <- do.call("cbind", trees[, 3])

    # calculate the final fit as a mean of all regression trees
    rf_fit <- apply(fits, MARGIN = 1, mean, na.rm = TRUE)
    rf_pred <- apply(preds, MARGIN = 1, mean, na.rm = TRUE)

    # # extract the feature importance
    # imp_full <- do.call("rbind", trees[, 2])
    #
    # # build the mean feature importance between all trees
    # imp <- aggregate(IMPORTANCE ~ FEATURES, FUN = mean, imp_full)
    #
    # # build the ratio for interpretation purposes
    # imp$IMPORTANCE <- imp$IMPORTANCE / sum(imp$IMPORTANCE)

    # export
    return(list(fit = rf_fit, pred = rf_pred))
}

#' Robust estimation of treatment effect using random forests.
#'
#' \code{rcate.rf} returns the robust treatment effect estimation model.
#'
#' This is a generic function:
#'
#' @param x matrix or a data frame of predictors.
#' @param y vector of response values.
#' @param d vector of binary treatment assignment (0 or 1).
#' @param method character string of CATE estimation method: "MCMEA" - modified covariate
#' method with efficiency augmentation, "RL" - R-learning, or "DR" - doubly robust method.
#' @param  algorithm character string of algorithm: "GBM" - gradient boosting machine or
#' "NN" - neural network.
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
#' @param n.trees.rf tuning parameter the number of trees used in GBM for estimating treatment
#' effect function if algorithm="GBM". The default is 1000.
#' @param feature.frac tuning parameter the number of interactions for estimating treatment
#' effect function if algorithm="GBM". The default value is 2.
#' @param newdata vector of the number of neurals in each hidden layer if algorithm='NN'.
#' The default is two layers with each layer the half size of previous layer.
#' @param minnodes vector of the dropout rate of each hidden layer if algorithm='NN'.
#' The default is no dropout.
#' @return a list of components
#' \itemize{
#'  \item fit - estimation method.
#'  \item pred - prediction of newdata.
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
#' tau_val = 6*sin(2*x_val[,1])+3*(x_val[,2]+3)*x_val[,3]+9*tanh(0.5*x_val[,4])+
#' 3*x_val[,5]*(2*I(x_val[,4]<1)-1)
#'
#' # Use MCM-EA transformation and GBM to estimate CATE
#' fit <- rcate.rf(X,y,d,newdata=data.frame(x_val),method='RL')
#' y_pred <- fit$pred
#' plot(tau_val,y_pred);abline(0,1)
#' @export
rcate.rf <- function(x, y, d, method = "MCMEA", algorithm = "GBM",
                  n.trees.p = 40000, shrinkage.p = 0.005, n.minobsinnode.p = 10,
                  interaction.depth.p = 1, cv.p = 5, n.trees.mu = c(1:50) * 50,
                  shrinkage.mu = 0.01, n.minobsinnode.mu = 5,
                  interaction.depth.mu = 5, cv.mu = 5,
                  n.trees.rf = 50, feature.frac = 1/2, newdata = NULL, minnodes = 5) {
    # Calculate T=2D-1
    t <- 2 * d - 1

    # Estimate mu0(x), mu1(x) and p(x)
    data.p <- data.frame(cbind(d, x))
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

    data00 <- data.frame(cbind(y, x, d))
    colnames(data00) <- c("y", paste0("X", 1:ncol(x)), "d")
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

    data.rf <-  data.frame(y.tr,x)
    if (is.null(newdata)) {newdata.rf <- data.rf} else {newdata.rf <- data.frame(newdata)}
    formula.rf <- stats::as.formula(paste0('y.tr', " ~ ",paste0(colnames(data.frame(x)), collapse = " + ")))
    result <- reg_rf(formula=formula.rf, n_trees = n.trees.rf, feature_frac = feature.frac,
           data=data.rf, newdata=newdata.rf, weights = w.tr, minnodes = 5)

    return(result)
}






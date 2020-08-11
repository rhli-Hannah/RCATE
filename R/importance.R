#' Importance level from machine learning algorithms.
#'
#' \code{importance.rcate.ml} Returns the variable importance level from "rcate.ml" model.
#'
#' @param object "rcate.ml" object.
#' @param plotit whether plot the importance level.
#' @param ... other.
#' @return a list of components
#' \itemize{
#'  \item importance - vector of variable importance level.
#'  }
#' @rdname importance.rcate.ml
#' @export
importance.rcate.ml <- function(object, plotit=TRUE,...) {
  importance <- object$importance
  algorithm <- object$algorithm

  x <- object$x
  col <- ncol(x)

  if (plotit==TRUE) {
    if (algorithm == "GBM") {
    importance <- data.frame(Variable=paste0('X',seq(1:col)),Importance = importance)
    graphics::barplot(importance$Importance,
            horiz = TRUE,names.arg = importance$Variable,
            main = 'Variable Importance from GBM',
            xlab = 'Importance', ylab= 'Variable')
  } else if (algorithm == "NN") {
    importance <- data.frame(importance)
    importance[,1] <- stringr::str_replace(importance[,1],'V','X')
    graphics::barplot(importance$Importance,
            horiz = TRUE,names.arg = importance$Variable,
            main = 'Variable Importance from Neural Network',
            xlab = 'Importance', ylab= 'Variable')
  }
  }
  return(importance)
}

#' Importance level from machine learning algorithms.
#'
#' \code{importance.rcate.rf} Returns the variable importance level from "rcate.rf" model.
#'
#' @param object "rcate.rf" object.
#' @param plotit whether plot the importance level.
#' @param ... other.
#' @return a list of components
#' \itemize{
#'  \item importance - vector of variable importance level.
#'  }
#' @rdname importance.rcate.rf
#' @export
importance.rcate.rf <- function(object, plotit=TRUE,...) {
 importance <- object$importance

 if (plotit==TRUE) {
   graphics::barplot(importance$IMPORTANCE,
           horiz = TRUE,names.arg = importance$FEATURES,
           main = 'Variable Importance from Random Forests',
           xlab = 'Importance', ylab= 'Variable')
 }
 colnames(importance) <- c('Variable','Importance')
 return(importance)
}

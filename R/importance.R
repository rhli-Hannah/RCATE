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
    importance <- data.frame(Importance=importance,name=names(importance))
    par(mgp=c(5,1,0))
    par(mar=c(8,8,4,2)+0.1)
    graphics::barplot(importance$Importance,
                      horiz=TRUE,names.arg = importance$name,
            main = 'Variable Importance from GBM',
            xlab = 'Importance', ylab= 'Variable', las=2, cex.names=0.8)
  } else if (algorithm == "NN") {
    importance <- data.frame(importance)
    par(mgp=c(5,1,0))
    par(mar=c(8,8,4,2)+0.1)
    graphics::barplot(importance$Importance,
            horiz = TRUE,names.arg = importance$Variable,
            main = 'Variable Importance from Neural Network',
            xlab = 'Importance', ylab= 'Variable', las=2, cex.names=0.8)
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

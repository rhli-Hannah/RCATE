
#' Plot performance history of machine learning algorithm.
#'
#' \code{perm.plot.rcate.ml} Plot performance history of "rcate.ml" model.
#'
#' @param object "rcate.ml" object.
#' @param ... other.
#' @rdname perm.plot.rcate.ml
#' @export
perm.plot.rcate.ml <- function(object,...) {
  model <- object$model
  algorithm <- object$algorithm
  if (algorithm == 'GBM') {
    gbm::gbm.perf(model,method = 'OOB', oobag.curve = TRUE, plot.it = TRUE)
  } else if (algorithm == 'NN') {
    history <- data.frame(object$history)
    graphics::plot(history$epoch,history$value,type='l',xlab='epoch',ylab = 'Absolute loss')
  }
}

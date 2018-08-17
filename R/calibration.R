#' Calibration of the observed probability vs. prediction probability
#' @description https://rdrr.io/cran/caret/man/calibration.html
#' @param model The prediction model.
#' @param data A dataframe containing features of the test samples.
#' @param class Class of the test samples.
#' @keywords
#' @export
#' @examples
#' calPlot()

calPlot <- function(model, data, class) {
  require(caret)
  require(ggplot2)
  prob <- data.frame(dat = class,
                     mod = predict(model, data)$posterior[, 1])
  calData <- calibration(dat ~ mod, data = prob)
  p <- ggplot(calData)
  return(p)
}

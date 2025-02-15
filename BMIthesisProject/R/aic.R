#' Computes the Akaike information criterion of an ordinary least square regression model.
#' The AIC was corrected for synthetically oversampled data, for which the sample size is unrepresentable
#'
#' @param model An ordinary least square regression model
#' @param oversampled A logical value indicating whether the data set used to train the regression model was oversampled or not.
#' @param amount The number of synthetically generated observations added to the oversampled data set
#' @return A numeric AIC value calculated as $n log(SSE) - n log(n) + 2p$
#' @export
#'
aic <- function(model, oversampled=FALSE, amount=0) {
  n <- length(model$coefficients)+model$df.residual
  SSE <- sum(model$residuals**2)
  p <- length(model$coefficients)

  if (oversampled) {
    n <- n-amount   # substract the amount of synthetical observations from the sample size
    SSE <- SSE*n/(n+amount)   # rescale the SSE to the original sample size
  }

  n*log(SSE) - n*log(n) + 2*p
}

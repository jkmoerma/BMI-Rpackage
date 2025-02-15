#' Estimate the fraction of category 2 + 95pc CI from a sample with two-category outcome counts.
#' Return result in numeric format
#'
#' @param n1 Number of observed category 1 outcomes
#' @param n2 Number of observed category 2 outcomes
#' @param reps Number of bootstrap replicates used to estimate the 95pc CI, defaults to 1000
#' @return An estimate of the fraction + 95pc CI as a numeric vector
#' @examples
#' calculateFractionNumeric(n1=90, n2=10)
#'
calculateFractionNumeric <- function(n1, n2, reps=1000) {
  fracs <- replicate(n=reps, {
    repCounts <- table(sample(factor(1:2), size=n1+n2, prob=c(n1, n2)/(n1+n2), replace=TRUE))
    repCounts[2]/(repCounts[1]+repCounts[2])
  })
  quantile(fracs, prob=c(0.025, 0.5, 0.975))
}

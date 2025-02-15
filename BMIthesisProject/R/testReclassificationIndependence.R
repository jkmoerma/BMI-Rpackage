#' Test for reclassification independence between two tables of reclassification healthy/syndromatic counts by ethnicity.
#'
#' @param classCounts1 A tibble or data frame containing columns "pred.: Healthy" and "pred.: Syndr."
#' @param classCounts2 A second table with the same number of rows and columns "pred.: Healthy" and "pred.: Syndr."
#' @return A chi-square p-value for the test of independence in reclassification between the two reclassification tables
#' @examples
#'
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' normals1 <- data.frame(ObesityClass=rep("Normal weight", times=50),
#'                       Reclassified=c("Healthy", "Syndromatic")[1+rbinom(n=50, size=1, prob=0.2)])
#' overweights1 <- data.frame(ObesityClass=rep("Overweight", times=50),
#'                           Reclassified=c("Healthy", "Syndromatic")[1+rbinom(n=50, size=1, prob=0.5)])
#' obeses1 <- data.frame(ObesityClass=rep("Obese", times=50),
#'                           Reclassified=c("Healthy", "Syndromatic")[1+rbinom(n=50, size=1, prob=0.8)])
#' data1 <- bind_rows(normals1, overweights1, obeses1)
#' classCounts1 <- countClassesNumeric(data1, race="Not specified", strat="Unknown")
#'
#' normals2 <- data.frame(ObesityClass=rep("Normal weight", times=50),
#'                       Reclassified=c("Healthy", "Syndromatic")[1+rbinom(n=50, size=1, prob=0.1)])
#' overweights2 <- data.frame(ObesityClass=rep("Overweight", times=50),
#'                           Reclassified=c("Healthy", "Syndromatic")[1+rbinom(n=50, size=1, prob=0.4)])
#' obeses2 <- data.frame(ObesityClass=rep("Obese", times=50),
#'                           Reclassified=c("Healthy", "Syndromatic")[1+rbinom(n=50, size=1, prob=0.8)])
#' data2 <- bind_rows(normals2, overweights2, obeses2)
#' classCounts2 <- countClassesNumeric(data2, race="Not specified", strat="Also Unknown")
#'
#' testReclassificationIndependence(classCounts1, classCounts2)
#'
testReclassificationIndependence <- function(classCounts1, classCounts2) {

  stopifnot('classCounts1 must contain columns "pred.: Healthy" and "pred.: Syndr."' = all(c("pred.: Healthy", "pred.: Syndr.") %in% names(classCounts1)))
  stopifnot('classCounts2 must contain columns "pred.: Healthy" and "pred.: Syndr."' = all(c("pred.: Healthy", "pred.: Syndr.") %in% names(classCounts2)))
  stopifnot("classCounts1 and classCounts2 must contain the same number of rows" = nrow(classCounts1)==nrow(classCounts2))

  # By obesity class the total amounts of patients reclassified as healthy and
  # syndromatic, pooled for the two tables
  totalPredHealthy <- classCounts1$`pred.: Healthy`+classCounts2$`pred.: Healthy`
  totalPredSyndr <- classCounts1$`pred.: Syndr.`+classCounts2$`pred.: Syndr.`

  # By obesity class the marginal reclassification rates for healthy and syndromatic
  relHealthy <- rep(totalPredHealthy/(totalPredHealthy+totalPredSyndr), times=2)
  relSyndr <- rep(totalPredSyndr/(totalPredHealthy+totalPredSyndr), times=2)

  # the estimated class counts for classCounts1 and classCounts2 under the independence assumption
  classCountsAll <- bind_rows(classCounts1, classCounts2)
  indepHealthy <- relHealthy*(classCountsAll$`pred.: Healthy`+classCountsAll$`pred.: Syndr.`)
  indepSyndr <- relSyndr*(classCountsAll$`pred.: Healthy`+classCountsAll$`pred.: Syndr.`)

  # Test of independence for classCounts1 and classCounts2
  # The degrees of freedom are calculated as:
  #         2 * number of obesity classes * 2 (all entries of classCounts1 and classCounts2)
  #           - 2 * number of obesity classes (totals of patients in every obesity class for every table is fixed)
  #           - number of obesity classes (syndromatic rates under independence)
  #         = number of obesity classes
  df <- nrow(classCounts1)
  X2 <- sum((classCountsAll$`pred.: Healthy`-indepHealthy)**2/indepHealthy) +
    sum((classCountsAll$`pred.: Syndr.`-indepSyndr)**2/indepSyndr)
  1-pchisq(X2, df=df)
}

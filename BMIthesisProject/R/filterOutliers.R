#' This function splits the data set in outliers and regular observations
#'
#' @param df A data frame as given by the owner of the data with imputed values, it must contain a column with unique patient identifiers, called "ID".
#' @param outliers A data frame of logical values indicating whether values of a metabolite were flagged as outliers, it must contain a column with unique patient identifiers, called "ID".
#' @param mets The metabolites used to filter out outliers
#' @return A list containing the data frames of "regulars" and "outliers"
#' @examples
#' df <- data.frame(ID=c("a", "b", "c", "d", "e"),
#'                  met1=c(2, 1, 10, 4, 3),
#'                  met2=c(20, 11, 12, 8, 9)
#'                  met3=c(0.5, 6, 0.1, 0.4, 1.1))
#' outliers <- data.frame(ID=c("a", "b", "c", "d", "e"),
#'                        met1=c(FALSE, FALSE, TRUE, FALSE, FALSE),
#'                        met2=c(TRUE, FALSE, FALSE, FALSE, FALSE)
#'                        met3=c(FALSE, TRUE, FALSE, FALSE, FALSE))
#' filterOutliers(df, outliers, mets=c("met1", "met2"))
#' filterOutliers(df, outliers, mets=c("met1", "met3"))
#'
filterOutliers <- function(df, outliers, mets) {
  subsetOutlier <- with(outliers, eval(parse(text=paste(mets, collapse="|"))))
  outlierIDs <- outliers$ID[which(subsetOutlier)]
  w <- which(df$ID %in% outlierIDs)
  list(regulars=df[-w,], outliers=df[w,])
}

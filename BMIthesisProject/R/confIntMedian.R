#' Takes a wilcox.test object with option "conf.int=TRUE" and returns the pseudo-median
#' difference between two population in character format.
#'
#' @param wilcoxTest a wilcox.test object with option "conf.int=TRUE"
#' @param negative Logical adressing whether the roles of the two groups must be reversed
#' @return An estimate of the pseudo median difference + 95pc CI
#' @examples
#'
#' dataA <- data.frame(group=rep("A", times=20),
#'                     val=rnorm(n=20, mean=0, sd=1))
#' dataB <- data.frame(group=rep("B", times=20),
#'                     val=rnorm(n=20, mean=1, sd=1))
#' data <- bind_rows(dataA, dataB)
#' wilcoxTest <- wilcox.test(formula=val~group, data=data, conf.int=TRUE)
#' confIntMedian(wilcoxTest, negative=FALSE)
#' confIntMedian(wilcoxTest, negative=TRUE)
#'
confIntMedian <- function(wilcoxTest, negative=FALSE) {
  if (negative) {
    return(sprintf("%.2f[%.2f %.2f]",
                   -wilcoxTest$estimate,
                   -wilcoxTest$conf.int[2],
                   -wilcoxTest$conf.int[1]))
  } else {
    return(sprintf("%.2f[%.2f %.2f]",
                   wilcoxTest$estimate,
                   wilcoxTest$conf.int[1],
                   wilcoxTest$conf.int[2]))
  }

}

#' Assess log-normality of metabolite distributions + tagging outliers
#'
#' @param df A data frame containing columns "ID", "Race", "BMI" and metabolite levels.
#' @return A list object containing a data frame with outliers for every metabolite, a count of outliers + frequency of common outliers with other metabolites, a table of metabolites which were not log-normally distributed + p-value of the Shapiro-Wilk test
#'
logNormalityAndOutliers <- function(df) {

  outliers <- data.frame(ID=df$ID, Race=df$Race)
  logNormality <- vector(mode="numeric", length=length(metabolites)+1)
  names(logNormality) <- c("BMI", metabolites)

  for (met in c("BMI", metabolites)) {

    # for standardizing log concentrations of metabolites
    standardize <- paste0("(log(", met, ")-mean(log(", met,
                          "), na.rm=TRUE))/sd(log(", met, "), na.rm=TRUE)")

    # for calculating theoretical quantiles under log-normality
    Zscore <- paste0("qqnorm(", standardize, ", plot.it=FALSE)$x")

    # for retrieving the smooth distribution value of all observations
    quantSmooth <- "smoother <- loess(formula=standardized~normQuantile, span=0.4,
                                     degree=2, family='symmetric');
                    predict(smoother, normQuantile)"

    # applying the above defined formulas on the data
    investigation <- df %>% group_by(Race) %>%
      mutate(standardized=eval(parse(text=standardize)),
             normQuantile=eval(parse(text=Zscore)),
             quantSmooth=eval(parse(text=quantSmooth)),
             outliers=(standardized-quantSmooth)*normQuantile/abs(normQuantile)>0.7)

    # store detected outliers
    outliers <- merge(outliers, investigation[,c("ID", "outliers")], by="ID")
    outliers[[met]] <- outliers$outliers
    outliers$outliers <- NULL

    # test normality
    logNormality[met] <- shapiro.test(investigation$standardized[which(!outliers[[met]])])$p.value

  }

  nonNormals <- logNormality[which(logNormality<1e-4)]
  totable <- cbind("met."=names(nonNormals), "p-val."=signif(nonNormals, digits=2))

  # create table to investigate dependencies in outliers between metabolites
  countOutliers <- data.frame(met=metabolites,
                              count=vector(mode="numeric", length=length(metabolites)),
                              commons=vector(mode="character", length=length(metabolites)),
                              BMIperc=vector(mode="character", length=length(metabolites)))
  BMIpercentages <- df %>% group_by(Race) %>% mutate(percent=percent_rank(BMI)*100)
  BMIpercentages <- BMIpercentages$percent
  for (i in 1:length(metabolites)) {
    met1 <- countOutliers$met[i]
    countOutliers$count[i] <- sum(outliers[[met1]])
    if (sum(outliers[[met1]], na.rm=TRUE)>0) {
      countOutliers$BMIperc[i] <- paste(paste0(round(BMIpercentages, digits=1)[which(outliers[[met1]])], "\\%"), collapse=" ")
    }
    commonSet <- c()
    for (j in (1:length(metabolites))[-i]) {
      met2 <- countOutliers$met[j]
      if (length(which(outliers[[met1]]&outliers[[met2]]))>0) {
        commonSet <- c(commonSet, paste0(met2,":", length(which(outliers[[met1]]&outliers[[met2]]))))
      }
    }
    if (length(commonSet) > 0) {
      countOutliers$commons[i] <- paste0(commonSet, collapse=" ")
    }
  }

  list(outliers=outliers, countOutliers=countOutliers, nonNormals=totable)

}

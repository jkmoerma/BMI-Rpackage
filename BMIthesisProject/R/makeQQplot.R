#' Makes QQ-plots of a subselection of metabolites, stored in the file "QQplots.pdf"
#'
#' @param df A data frame containing columns "ID", "BMI", "Race" and metabolite levels.
#' @param mets A subselection of metabolites to make a QQ-plot of
#' @return NULL
#'
makeQQplot <- function(df, mets, pvals) {

  outliers <- data.frame(ID=df$ID, Race=df$Race)
  logNormality <- vector(mode="numeric", length=length(metabolites)+1)
  names(logNormality) <- c("BMI", metabolites)
  colors <- c("red","blue","green","violet","black")
  chars <- c("W", "B", "S", "E", "M")

  pdf(file="QQplots.pdf", width=8.27, height=11.69)
  layout(matrix(1:6, byrow=TRUE, ncol=2))
  par(mai=c(0.6, 0.6, 0.4, 0.2), mex=1)
  for (i in 1:length(mets)) {

    met <- mets[i]
    pval <- pvals[i]

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

    # plot QQ-plot of metabolite distribution
    normRangeX <- range(investigation$normQuantile, na.rm=TRUE)
    normRangeY <- range(investigation$standardized, na.rm=TRUE)
    plot(x=normRangeX, y=normRangeY, col="white", main=met,
         xlab="theoretical quantiles", ylab="standardized log concentrations")
    abline(a=0, b=1, lty="dashed", col="darkgrey")
    for (j in 1:length(ethnicities)) {
      race <- ethnicities[j]
      normPlot <- subset(investigation, subset = Race==race)
      plotchar <- rep(chars[j], times=nrow(normPlot))
      w <- which(normPlot$outliers)
      plotchar[w] <- "."
      points(x=normPlot$normQuantile, y=normPlot$standardized, col=colors[j], pch=plotchar)
      if (length(w)>0) {
        text(x=normPlot$normQuantile[w], y=normPlot$standardized[w], col=colors[j],
             labels=normPlot$ID[w], cex=0.7)
      }
    }

    # result normality test
    text(x=normRangeX[1], y=normRangeY[2], col="red", adj=c(0,1),
         labels=paste0("Shap.-Wilk p-val.: %.1e", pval))
    #labels=sprintf("Shap.-Wilk p-val.: %.1e", pval))

  }
  dev.off()


}

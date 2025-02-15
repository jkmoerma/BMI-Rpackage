#' Stores a plot to demonstrate outlier detection with metabolites and ethnicities of interest, stored in the file "outlierDemonstration.pdf"
#'
#' @param df A data frame containing columns "ID", "BMI", "Race" and metabolite levels.
#' @param mets Metabolite names of interest
#' @param ethnicity Ethnicities according to the metabolites of interest
#'
demonstrateOutlierDetection <- function(df, mets, ethnicity) {

  stopifnot("arguments 'mets' and 'ethnicity' must have the same length" =
              length(mets)==length(ethnicity))

  outliers <- data.frame(ID=df$ID, Race=df$Race)
  colors <- c("red","blue","green","violet","black")
  fillRegular <- c("indianred","cadetblue","lightgreen","lightpink","grey")
  chars <- c("W", "B", "S", "E", "M")
  if (length(mets)%%2 == 0) {
    height <- 11.69 * length(mets)/6
    plot_pos <- 1:length(mets)
  }
  if (length(mets)%%2 == 1) {
    height <- 11.69 * (length(mets)+1)/6
    plot_pos <- c(1:length(mets), 0)
  }
  pdf(file="outlierDemonstration.pdf", width=8.27, height=height)
  layout(matrix(plot_pos, byrow=TRUE, ncol=2))
  par(mai=c(0.6, 0.6, 0.4, 0.2), mex=0.5)
  for (i in 1:length(mets)) {
    met <- mets[i]
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
      if (race == ethnicity[i]) {
        lines(x=sort(normPlot$normQuantile), y=sort(normPlot$quantSmooth), col=colors[j])
        polygon(x=c(sort(normPlot$normQuantile, decreasing=FALSE), sort(normPlot$normQuantile, decreasing=TRUE)),
                y=c(sort(normPlot$quantSmooth, decreasing=FALSE)-0.7,
                    sort(normPlot$quantSmooth, decreasing=TRUE)+0.7),
                col=fillRegular[j], border="white")
      }
      points(x=normPlot$normQuantile, y=normPlot$standardized, col=colors[j], pch=plotchar)
      if (length(w)>0) {
        text(x=normPlot$normQuantile[w], y=normPlot$standardized[w], col=colors[j],
             labels=normPlot$ID[w], cex=0.7)
      }
    }

  }
  dev.off()
}

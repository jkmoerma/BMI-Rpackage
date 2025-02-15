#' Calculates the difference in metabolic levels between ON, NO, OO and NN.
#' A plot is made of the scaled metabolic differences for comparison.
#'
#' @param data Data with metabolite levels and assigned prediction group "predictionGroup"
#' @param ethnicity Ethnicity of the patients in the data set
#' @return A list of 2 ggplot objects for NO ("normal") and ON ("obese") outliers
#' @examples
#'
#' metabolites <- c("met1", "met2", "met_coll", "met_110")
#'
#' model <- trainRidgeLASSO(effect="main", type="Ridge", transformation="Log", new.data=df)
#' df$predicted <-
#'   exp(predict(model, newx=makeMatrix(df, includeInteraction=FALSE)$mat)[,"s0"])
#'
#' rocTrain <- roc(df$ObesityClass, df$predicted,
#'                 levels=c("Normal weight", "Obese"))
#' cutoff <- coords(rocTrain, x="best")[1,"threshold"]
#' df <- df %>%
#'   mutate(predictionGroup = ifelse(ObesityClass=="Normal weight"&predicted<cutoff, "NN",
#'                              ifelse(ObesityClass=="Normal weight"&predicted>cutoff, "NO",
#'                                ifelse(ObesityClass=="Obese"&predicted>cutoff, "OO",
#'                                  ifelse(ObesityClass=="Obese"&predicted<cutoff, "ON", NA)))))
#'
#' ggplot(df, aes(x=predicted, y=exp(BMI), col=predictionGroup)) +
#'   geom_point() +
#'   labs(x="predicted BMI", y="observed BMI")
#'
#' plotANOVA(df, "White")
#'
plotANOVA <- function(data, ethnicity) {
  for (met in c("Age", metabolites)) {
    data[[met]] <- scale(data[[met]])
  }
  groupLevels <- subset(data, subset=!is.na(predictionGroup),
                        select=c("Age", metabolites, "predictionGroup"))

  levelDiffsNormal <- matrix(nrow=length(metabolites)+1, ncol=6)
  rownames(levelDiffsNormal) <- c("Age", metabolites)
  colnames(levelDiffsNormal) <- c("NO-NN", "NO-NN (lCI)", "NO-NN (uCI)",
                                  "OO-NO", "OO-NO (lCI)", "OO-NO (uCI)")

  levelDiffsObese <- matrix(nrow=length(metabolites)+1, ncol=6)
  rownames(levelDiffsObese) <- c("Age", metabolites)
  colnames(levelDiffsObese) <- c("OO-ON", "OO-ON (lCI)", "OO-ON (uCI)",
                                 "NN-ON", "NN-ON (lCI)", "NN-ON (uCI)")

  legend <- c()
  for (met in c("Age", metabolites)) {
    diffs <- TukeyHSD(aov(groupLevels[[met]]~groupLevels$predictionGroup))$`groupLevels$predictionGroup`
    levelDiffsNormal[met,] <- c(diffs["NO-NN", "diff"], diffs["NO-NN", "lwr"],
                                diffs["NO-NN", "upr"], diffs["OO-NO", "diff"],
                                diffs["OO-NO", "lwr"], diffs["OO-NO", "upr"])
    levelDiffsObese[met,] <- c(diffs["OO-ON", "diff"], diffs["OO-ON", "lwr"],
                               diffs["OO-ON", "upr"], -diffs["ON-NN", "diff"],
                               -diffs["ON-NN", "upr"], -diffs["ON-NN", "lwr"])
    if (abs(diffs["OO-NN", "diff"])>0.5) {
      legend <- c(legend, met)
    } else {
      legend <- c(legend, "other")
    }
  }

  levelDiffsNormal <- data.frame(met=rownames(levelDiffsNormal),
                                 legend = legend,
                                 as.data.frame(levelDiffsNormal),
                                 check.names=FALSE)
  levelDiffsObese <- data.frame(met=rownames(levelDiffsObese),
                                legend=legend,
                                as.data.frame(levelDiffsObese),
                                check.names=FALSE)

  pN <- ggplot(data=levelDiffsNormal,
               aes(x=`OO-NO`, y=`NO-NN`, label=met, col=met)) +
    theme(legend.position = "none") +
    coord_fixed(ratio=1) +
    geom_abline(slope=0, intercept=0, lty="dashed") +
    geom_abline(slope=-1, intercept=0, lty="dashed") +
    geom_vline(xintercept=0, lty="dashed") +
    geom_point() +
    #geom_pointrange(aes(ymin=`NO-NN (lCI)`, ymax=`NO-NN (uCI)`), fatten=1, lwd=0.1) +
    #geom_pointrange(aes(xmin=`OO-NO (lCI)`, xmax=`OO-NO (uCI)`), fatten=1, lwd=0.1) +
    geom_text(vjust=0, nudge_y=0.05, size=2) +
    labs(title=paste0("NO (", ethnicity, " ethnicity)"),
         x="Scaled log conc. diff. OO-NO",
         y="Scaled log conc. diff. NO-NN")

  pO <- ggplot(data=levelDiffsObese,
               aes(x=`NN-ON`, y=`OO-ON`, label=met, col=met)) +
    theme(legend.position = "none") +
    coord_fixed(ratio=1) +
    geom_abline(slope=0, intercept=0, lty="dashed") +
    geom_abline(slope=1, intercept=0, lty="dashed") +
    geom_vline(xintercept=0, lty="dashed") +
    geom_point() +
    #geom_pointrange(aes(ymin=`OO-ON (lCI)`, ymax=`OO-ON (uCI)`), fatten=1, lwd=0.1) +
    #geom_pointrange(aes(xmin=`NN-ON (lCI)`, xmax=`NN-ON (uCI)`), fatten=1, lwd=0.1) +
    geom_text(vjust=0, nudge_y=0.05, size=2) +
    labs(title=paste0("ON (", ethnicity, " ethnicity)"),
         x="Scaled log conc. diff. NN-ON",
         y="Scaled log conc. diff. OO-ON")

  list(normal=pN, obese=pO)
}

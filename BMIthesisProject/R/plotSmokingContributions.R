#' Visualize the main metabolite contributions in the discrepancies in prediction
#' between smokers and non-smokers.
#'
#' @param data_smoking A data frame with metabolite measurements and predicted BMI of smoking patients
#' @param data_test A data frame with metabolite measurements and predicted BMI of non-smoking test set patients
#' @return A ggplot object visualizing the main metabolic contributions of smoking to the discrepancy in predicted BMI
#'
plotSmokingContributions <- function(data_smoking, data_test, model) {
  metDiff <- vector(mode="numeric", length=length(metabolites)+1)
  names(metDiff) <- c("Age", metabolites)
  metDiffSd <- vector(mode="numeric", length=length(metabolites)+1)
  names(metDiffSd) <- c("Age", metabolites)
  for (met in c("Age", metabolites)) {
    nonSmokers <- subset(data_test, subset = ObesityClass=="Normal weight")[[met]]
    smokers <- subset(data_smoking, subset = ObesityClass=="Normal weight")[[met]]
    mean1 <- mean(nonSmokers)
    mean2 <- mean(smokers)
    met_beta <- model$beta[met, "s0"]
    metDiff[met] <- (mean(smokers)-mean(nonSmokers))*met_beta
    metDiffSd[met] <- sqrt(var(smokers)/length(smokers)+var(nonSmokers)/length(nonSmokers))*met_beta

  }

  fracs <- data.frame(met=names(metDiff),
                      diffs=metDiff,
                      lCI=metDiff-1.96*metDiffSd,
                      uCI=metDiff+1.96*metDiffSd)
  w <- order(abs(fracs$diffs), decreasing=TRUE)   # from large effects to small effects
  fracs <- fracs[w,]
  fracs <- fracs[1:20,]

  ggplot(data=fracs,
         aes(x=reorder(met,abs(diffs)), y=diffs, ymin=lCI, ymax=uCI)) +
    geom_pointrange(stat="identity") +
    geom_hline(yintercept=0, lty="dashed") +
    #theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          legend.title = element_blank()) +
    xlab("") + ylab("contribution") +
    coord_flip()
}

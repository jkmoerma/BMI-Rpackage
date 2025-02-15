#' Illustrates the interpretation of the figures on metabolic differences of the ON outliers
#'
#' @return A ggplot object
illustrON <- function() {
  ggplot(data=data.frame("OO-NO"=c(0,-1, 1, -2), "NO-NN"=c(1,0, 2, -1), check.names=FALSE),
         mapping=aes(x=`OO-NO`, y=`NO-NN`)) +
    geom_point() +
    geom_abline(intercept=0, slope=1, lty="dashed") +
    geom_abline(intercept=1, slope=1, lty="dashed") +
    geom_vline(xintercept=0) +
    geom_hline(yintercept=0) +
    geom_segment(data=data.frame("OO-NO"=0, "NO-NN"=0,
                                 "OO-NO end"=0, "NO-NN end"=1,
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0.1, "inches"),
                             ends="both", type="open"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="purple")) +
    geom_text(data=data.frame("OO-NO"=0, "NO-NN"=0.5,
                              check.names=FALSE),
              aes(col="purple"),
              label="OO-NN",
              nudge_x= -0.05, angle=90) +
    geom_segment(data=data.frame("OO-NO"=c(0,-1), "NO-NN"=c(1,0),
                                 "OO-NO end"=c(0, -1)-0.4, "NO-NN end"=c(1,0)+0.4,
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0.1, "inches"),
                             ends="last", type="closed"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="blue")) +
    geom_text(data=data.frame("OO-NO"=c(0, -1)-0.4, "NO-NN"=c(1, 0)+0.4,
                              check.names=FALSE),
              aes(col="blue"),
              label=c("ideal model", "random guessing BMI"),
              nudge_x= -0.15, nudge_y= 0.15) +
    geom_segment(data=data.frame("OO-NO"=-1, "NO-NN"=0,
                                 "OO-NO end"=0, "NO-NN end"=1,
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0, "inches"),
                             ends="last", type="closed"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="green")) +
    geom_text(data=data.frame("OO-NO"=-0.5, "NO-NN"=0.5,
                              check.names=FALSE),
              aes(col="green"),
              label="suboptimal model",
              nudge_x= -0.1, nudge_y= 0.1, angle=45) +
    geom_segment(data=data.frame("OO-NO"=c(-1,0), "NO-NN"=c(0,1),
                                 "OO-NO end"=c(-2, 1), "NO-NN end"=c(-1,2),
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0, "inches"),
                             ends="last", type="closed"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="red")) +
    geom_text(data=data.frame("OO-NO"=c(0.5, -1.5), "NO-NN"=c(1.5, -0.5),
                              check.names=FALSE),
              aes(col="red"),
              label=c("abberant in ON", "abberant in ON"),
              nudge_x= -0.1, nudge_y= 0.1, angle=45) +
    scale_color_manual(values=c("blue", "green", "purple", "red")) +
    theme(legend.position = "none") +
    coord_fixed(ratio=1) +
    labs(title="ON outliers: observations and expectations", x="NN-ON", y="OO-ON")
}

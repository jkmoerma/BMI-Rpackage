#' Train a stepwise OLS model from the data with a specification about main/interaction effects and a BMI transformation.
#'
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "OLS" for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "balanced" if oversampling is needed.
#' @param new.data Data to build a regression model from
#' @return An OLS regression model
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' metabolites <- c("met1", "met2", "met_110", "met_coll")
#'
#' trainOLS(effect="main", type="OLS", transformation="Log", balancing="balanced", new.data=df)
#'
trainOLS <- function(effect, type, transformation, balancing="", new.data) {

  n0 <- nrow(new.data)
  if (balancing=="Balanced") {new.data <- oversample(new.data)}
  amount <- nrow(new.data) - n0

  new.data$transBMI <- new.data$BMI
  if (transformation=="Inv") {new.data$transBMI <- exp(-new.data$BMI)}

  interactionEffect <- FALSE
  if (effect=="interaction") {
    interactionEffect <- TRUE
  }

  # main effects, log(BMI), ethnicity="White"
  intercept_only <- lm(transBMI ~ 1, data=new.data)
  all <- lm(data=new.data,
            formula = eval(parse(text=paste0("transBMI~`", paste(c("Age", metabolites),
                                                                 collapse="`+`"), "`"))))
  model <- step(intercept_only, direction='both', scope=formula(all), trace=0,
                k=2*(n0+amount)/n0)

  if (interactionEffect) {
    # interaction effects, log(BMI), ethnicity="White"
    variables <-
      colnames(model$model)[-which(colnames(model$model)=="transBMI")]
    interactionSet <- c()
    for (i in 1:(length(variables)-1)) {
      for (j in (i+1):length(variables)) {
        term <- paste0(variables[i], "`*`", variables[j])
        addedModel <- lm(data=new.data,
                         formula=eval(parse(text=paste0("transBMI~`",
                                                        paste0(c(variables, term),
                                                               collapse="`+`"), "`"))))
        AIC0 <- aic(model, oversampled=balancing=="Balanced", amount=amount)
        AIC1 <- aic(addedModel, oversampled=balancing=="Balanced", amount=amount)
        if (AIC1 < AIC0) {interactionSet <- c(interactionSet, term)}
      }
    }
    all <- lm(data=new.data,
              formula = eval(parse(text=paste0("transBMI~`", paste0(c(variables, interactionSet),
                                                                    collapse="`+`"), "`"))))
    model <- step(model, direction='both', scope=formula(all), trace=0,
                  k=2*(n0+amount)/n0)
  }
  model
}

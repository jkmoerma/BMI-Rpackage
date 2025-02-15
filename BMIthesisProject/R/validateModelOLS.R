#' Use 4-fold cross validation inside new.data to score a regression formulation with type OLS.
#' Calculate the AUC of every fold building a model from the other 3 folds.
#' Return the average AUC over the 4 folds.
#'
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "OLS" for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "Balanced" if oversampling is needed.
#' @param new.data Data to build a regression model from
#' @return An averaged AUC value representing a score of the model formulation
#'
validateModelOLS <- function(effect, type, transformation, balancing, new.data) {

  stopifnot("Function is for Ridge and LASSO regression only" = type=="OLS")

  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}

  outcome <- "BMI"
  if (transformation=="Inv") {outcome <- "exp(-BMI)"}

  # initiate data frame collecting data + predictions
  AUC <- 0

  # split data in 4 folds
  set.seed(4)
  foldid <- sample(1:4, size=nrow(new.data), replace=TRUE)

  # for every fold:
  for (i in 1:4) {
    w <- which(foldid==i)
    # train on other folds
    data_train <- new.data[-w, ]
    n0 <- nrow(data_train)
    if (balancing=="Balanced") {data_train <- oversample(data_train)}
    amount <- nrow(data_train) - n0

    intercept_only <- lm(eval(parse(text=paste0(outcome, " ~ 1"))), data=data_train)
    all <- lm(data=data_train,
              formula = eval(parse(text=paste0(outcome, "~`", paste(c("Age",  metabolites),
                                                                    collapse="`+`"), "`"))))
    mainModel <- step(intercept_only, direction='both', scope=formula(all), trace=0,
                      k=2*(n0+amount)/n0)

    if (effect=="main") {model <- mainModel}
    if (effect=="interaction") {
      variables <-
        colnames(mainModel$model)[-which(colnames(mainModel$model)==outcome)]
      interactionSet <- c()
      for (i in 1:(length(variables)-1)) {
        for (j in (i+1):length(variables)) {
          term <- paste0(variables[i], "`*`", variables[j])
          addedModel <- lm(data=data_train,
                           formula=eval(parse(text=paste0(outcome, "~`",
                                                          paste0(c(variables, term),
                                                                 collapse="`+`"), "`"))))
          AIC0 <- aic(mainModel, oversampled=balancing=="Balanced", amount=amount)
          AIC1 <- aic(addedModel, oversampled=balancing=="Balanced", amount=amount)
          if (AIC1 < AIC0) {interactionSet <- c(interactionSet, term)}
        }
      }
      all <- lm(data=data_train,
                formula = eval(parse(text=paste0(outcome, "~`", paste0(c(variables, interactionSet),
                                                                       collapse="`+`"), "`"))))
      model <- step(mainModel, direction='both', scope=formula(all), trace=0,
                    k=2*(n0+amount)/n0)
    }

    # predict left out fold
    foldPreds <- predict(object=model, newdata=new.data[w, ])

    # store error classes
    if (transformation=="Log") {
      obs <- exp(new.data[w, "BMI"])
      preds <- exp(foldPreds)
    }
    if (transformation=="Inv") {
      obs <- exp(new.data[w, "BMI"])
      preds <- 1/foldPreds
    }
    w_normal <- which(obs < 25)
    w_obese <- which(obs > 30)
    roc_j <- roc(controls=preds[w_normal], cases=preds[w_obese])

    AUC <- AUC + auc(roc_j)
  }

  return(round(AUC/4, digits=3))

}

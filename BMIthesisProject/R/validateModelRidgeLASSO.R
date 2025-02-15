#' Use 4-fold cross validation inside new.data to score a regression formulation with type Ridge or LASSO.
#' Calculate the AUC of every fold building a model from the other 3 folds.
#' Return the average AUC over the 4 folds.
#'
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "Ridge" or "LASSO for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "Balanced" if oversampling is needed.
#' @param new.data Data to build a regression model from
#' @return An averaged AUC value representing a score of the model formulation
#'
validateModelRidgeLASSO <- function(effect, type, transformation, balancing,
                                    new.data) {

  stopifnot("Function is for Ridge and LASSO regression only" = type=="Ridge"|type=="LASSO")

  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}

  includeInteraction <- FALSE
  if (effect=="interaction") {includeInteraction <- TRUE}

  outcomes <- new.data$BMI
  if (transformation=="Inv") {outcomes <- exp(-new.data$BMI)}

  # initiate pooled AUC
  AUC <- 0

  # split data in 4 folds
  set.seed(4)
  foldid <- sample(1:4, size=nrow(new.data), replace=TRUE)

  # for every fold:
  for (i in 1:4) {
    w <- which(foldid==i)

    # train on other folds (no lambda specified, therefore tuned)
    pruneModel <- trainRidgeLASSO(effect, type, transformation,
                                  new.data=new.data[-w,])

    # predict left out fold
    RidgeLASSOtest <- makeMatrix(new.data[w,], includeInteraction)
    foldPreds <- predict(object=pruneModel,
                         newx=RidgeLASSOtest$mat,
                         type="response")[,"s0"]

    # store error classes
    if (transformation=="Log") {
      obs <- exp(outcomes[w])
      preds <- exp(foldPreds)
    }
    if (transformation=="Inv") {
      obs <- 1/outcomes[w]
      preds <- 1/foldPreds
    }
    w_normal <- which(obs < 25)
    w_obese <- which(obs > 30)
    roc_j <- roc(controls=preds[w_normal], cases=preds[w_obese])

    AUC <- AUC + auc(roc_j)

  }

  return(round(AUC/4, digits=3))

}

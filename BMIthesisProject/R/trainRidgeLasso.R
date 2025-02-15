#' Train a Ridge or LASSO model from the data with a specification about main/interaction effects and a BMI transformation.
#'
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "Ridge" or "LASSO" for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "Balanced" if oversampling is needed.
#' @param new.data Data to build a regression model from
#' @return A glmnet regression model with 1 lambda parameter
#' @examples
#' metabolites <- c("met1", "met2", "met_coll", "met_110")
#'
#' data_train <- oversample(df)
#' trainRidgeLASSO(effect="main", type="Ridge", transformation="Log", new.data=data_train)
#'
trainRidgeLASSO <- function(effect, type, transformation,
                            new.data=NULL, new.x=NULL, new.y=NULL, lambda=NULL) {

  if (!is.null(new.x)) {
    stopifnot("Outcomes must be provided separately if data matrix is used" =
                !is.null(new.y))
  }

  interactionEffect <- FALSE
  if (effect=="interaction") {interactionEffect <- TRUE}

  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}

  if (is.null(new.x)) {
    generateX <- makeMatrix(new.data, includeInteraction = effect=="interaction")
    regression_data <- generateX$mat
    interactions <- generateX$interactions
    new.y <- new.data$BMI
    outcomes <- new.data$BMI
    if (transformation=="Inv") {outcomes <- exp(-new.data$BMI)}
  } else {
    regression_data <- new.x
    outcomes <- new.y
    if (transformation=="Inv") {outcomes <- exp(-new.y)}
  }

  # if not provided, tune a lambda parameter for the final model
  if (is.null(lambda)) {
    lambda <- tuneLambda(regression_data, alpha, new.y, effect,
                         transformation, returnAll=FALSE)
  }

  # return model trained on all the received data with tuned lambda parameter
  glmnet(x=regression_data,
         y=outcomes,
         alpha=alpha,
         lambda=lambda,
         family="gaussian")
}

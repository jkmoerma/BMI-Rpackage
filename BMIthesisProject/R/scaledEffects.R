#' Calculates regression coefficients of a model formulation and bootstrap replicates of the data set.
#' The predictors are scaled to unit variance in order to compare the regression coefficients of different predictors equally
#'
#' @param data Data to draw bootstrap replicates from
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "Ridge" or "LASSO for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "Balanced" if oversampling is needed.
#' @param boot.n the amount of bootstrap replicates to be drawn
#' @return A data tibble with regression coefficients for every replicate.
#'
scaledEffects <- function(data, effect, type, transformation, balancing="",
                          boot.n=100) {

  # if balancing the data set was part of the model formulation, scaling data on an oversampled replicate was necessary for ...
  if (balancing=="") {data_ref <- data}
  if (balancing=="Balanced") {data_ref <- oversample(data)}

  # specify alpha for the type of regression
  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}

  # generate a scaled version of the data set
  # do not standardize BMI on type of transformation
  data_scaled <- data
  for (met in c("Age", metabolites)) {
    data_scaled[[met]] <- data_scaled[[met]]/sd(data_ref[[met]])
  }

  if (balancing=="") {
    tune_outcome <- data_scaled$BMI
    dataMatrix <- makeMatrix(data_scaled, includeInteraction = effect=="interaction")$mat
  }
  if (balancing=="Balanced") {
    tune_data <- oversample(data_scaled)
    tune_outcome <- tune_data$BMI
    dataMatrix <- makeMatrix(tune_data, includeInteraction = effect=="interaction")$mat
  }

  lambda <- tuneLambda(dataMatrix, alpha = alpha, new.y=tune_outcome,
                       effect, transformation, returnAll=FALSE)

  # for boot.n bootstrap replicates, recalculate the regression coefficients
  # Register parallel backend
  numCores <- detectCores() - 1  # Use one less than the total number of cores
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  # Parallelized loop
  coeffs <- foreach(i=1:boot.n, .combine=bind_rows, .packages=c("dplyr", "glmnet"),
                    .export=c("trainOLS", "trainRidgeLASSO",
                              "metabolites", "aic", "oversample",
                              "SMOTE", "makeMatrix", "roc", "auc")) %dopar% {
                                set.seed(i)
                                # retrieve stratified bootstrap replicate
                                w <- sample(1:nrow(data_scaled), size=nrow(data_scaled), replace=TRUE)
                                data_rep <- data_scaled[w,]

                                if (type=="OLS") {
                                  if (balancing=="") {
                                    model <- trainOLS(effect, type, transformation, data_rep)
                                  }
                                  if (balancing=="Balanced") {
                                    model <- trainOLS(effect, type, transformation, oversample(data_rep))
                                  }
                                  return(model$coefficients)
                                }
                                if (type=="LASSO"|type=="Ridge") {
                                  if (balancing=="") {
                                    new.x <- dataMatrix[w,]
                                    new.y <- data_scaled$BMI[w]
                                    model <- trainRidgeLASSO(effect, type, transformation, new.x=new.x,
                                                             new.y=new.y, lambda=lambda)
                                  }
                                  if (balancing=="Balanced") {
                                    new.x <- makeMatrix(data_rep, includeInteraction = effect=="interaction")$mat
                                    new.y <- data_rep$BMI
                                    model <- trainRidgeLASSO(effect, type, transformation, new.x=new.x,
                                                             new.y=new.y, lambda=lambda)
                                  }
                                  return(model$beta[,"s0"])
                                }
                              }

  # Stop the cluster
  stopCluster(cl)

  # return coefficients
  return(coeffs)

}

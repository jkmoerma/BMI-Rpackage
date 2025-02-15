#' Hyper parameter tuning of a penalty lambda parameter for LASSO and ridge regression. A 4-fold cross validation is performed with the received data.
#' For every fold, a model is trained on the 3 other folds with 100 probed lambda parameters to predict the values of the left-out fold.
#' Per fold, the AUC is calculated for normal weight versus obese patients for every probed lambda.
#' The AUC values are averaged over the 4 folds for each probed lambda.
#' The lambda with the highest AUC fold average is chosen
#'
#' @param data_matrix Measured metabolite values and possibly in matrix format
#' @param alpha A parameter passed to glmnet. Set to 0 for ridge regression and to 1 for LASSO regression
#' @param new.y A vector having the same length as the number of rows of data_matrix. It represents the log(BMI) values according the data in data_matrix
#' @param effect "main" or "interaction" dependent on which effects are wished in the model
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param returnAll Logical. If FALSE, return the tuned lambda parameter, if TRUE return all probed lambdas and probed AUCs as a data frame
#' @return If returnAll=TRUE, returns a data frame with all probed lanmba values and the according AUC to distinguish normal weight and obese patients.
#'         If returnAll=FALSE, returns a single numeric value for the best lambda
#' @examples
#'
#' metabolites <- c("met1", "met2", "met_coll", "met_110")
#'
#' df_oversample <- oversample(df)
#' AUCs <- tuneLambda(data_matrix=makeMatrix(df_oversample)$mat,
#'                    new.y = df_oversample$BMI, alpha=0,
#'                    effect="main", transformation="Log",
#'                    returnAll=TRUE)
#' bestLambda <- AUCs$lambda[which.max(AUCs$AUC)]
#' plot(AUC~log(lambda), AUCs)
#' abline(v=log(bestLambda), col="blue")
#'
tuneLambda <- function(data_matrix, alpha, new.y, effect,
                       transformation, returnAll=FALSE) {

  # divide data_matrix in 10 folds for cross-validation
  set.seed(15)
  tuneSelect <- sample(x=1:nrow(data_matrix), size=round(0.25*nrow(data_matrix)))

  # define outcomes for regression
  if (transformation=="Log") {outcomes <- new.y}
  if (transformation=="Inv") {outcomes <- exp(-new.y)}

  ## calculate a start set of lambda values

  # the highest lambda parameter meaningful for the model
  model0 <- glmnet(x=data_matrix[-tuneSelect, ],
                   y=outcomes[-tuneSelect],
                   alpha=alpha,
                   nlambda=3,
                   lambda.min.ratio=0.1,
                   family="gaussian")
  max.lambda <- max(model0$lambda)

  # a lower limit lambda parameter defined by where the percentage deviance explained is almost 1
  min.lambda.series <- max.lambda * 10**(-(0:10))
  model1 <- glmnet(x=data_matrix[-tuneSelect, ],
                   y=outcomes[-tuneSelect],
                   alpha=alpha,
                   lambda=min.lambda.series,
                   family="gaussian")
  minIndex <- min(which(1-model1$dev.ratio/max(model1$dev.ratio)<1e-2)) + 2
  min.lambda <- model1$lambda[minIndex]

  # a sequence of desired lambda values from the model
  lambda_sequence <- exp(seq(from=log(max.lambda), to=log(min.lambda), length=100))

  # define folds for lambda cross validation
  set.seed(7)
  foldid <- sample(1:4, size=nrow(data_matrix), replace=TRUE)

  AUCs <- 0
  for (i in 1:4) {
    w <- which(foldid==i)

    # train models with probed lambda values
    model <- glmnet(x=data_matrix[-w, ],
                    y=outcomes[-w],
                    alpha=alpha,
                    lambda=lambda_sequence,
                    family="gaussian")
    preds <- predict(object=model,
                     newx=data_matrix[w, ],
                     type="response")
    BMIs <- exp(new.y)[w]

    # calculate AUC from predictions of the left-out data at the different lambda values probed
    AUCs_i <- sapply(X=1:model$dim[2],
                     FUN=function(j) {
                       w_normal <- which(BMIs < 25)
                       w_obese <- which(BMIs > 30)
                       roc_j <- roc(controls=preds[w_normal,j],
                                    cases=preds[w_obese,j])
                       auc(roc_j)
                     }
    )
    AUCs <- AUCs + AUCs_i
  }
  AUCs <- AUCs/4


  # return lambda (returnAll=FALSE) or whole evaluation (returnAll=TRUE)
  if (returnAll) {
    return(data.frame(lambda=model$lambda, AUC=AUCs))
  } else {
    # choose lambda parameter tuned to the highest cross-validated AUC
    tuneIndex <- which.max(AUCs)
    lambda <- model$lambda[tuneIndex]
    return(lambda)
  }

}

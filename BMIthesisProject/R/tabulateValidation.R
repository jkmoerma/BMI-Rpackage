#' Use 4-fold cross validation inside new.data to score the provided regression formulations
#' For every regression formulation, calculate the AUC of every fold building a model from the other 3 folds.
#' For every regression formulation, return the average AUC over the 4 folds.
#'
#' @param effects a vector of "main" and/or "interaction" effects of every formulation
#' @param types a vector of regression types "OLS", "Ridge" and/or "LASSO" of every formulation
#' @param transformation a vector of BMI transformations "Log" and/or "Inv" of every formulation
#' @param balancing defaults to NULL, or a vector of "" and/or "Balanced" reflecting balancing requests of every formulation
#' @param new.data Data to build a regression model from
#' @return An averaged AUC value representing a score of the model formulation
#' @examples
#' metabolites <- c("met1", "met2", "met_coll", "met_110")
#'
#' formulations <- expand.grid(effects=c("main", "interaction"),
#'                             types=c("OLS", "Ridge", "LASSO"),
#'                             transformations=c("Log", "Inv"),
#'                             balancing=c("", "Balanced"),
#'                             stringsAsFactors=FALSE)
#'
#' formulationScores <- tabulateValidation(formulations$effects, formulations$types,
#'                                         formulations$transformations,
#'                                         balancing=formulations$balancing, df)
#'
#' View(formulationScores)
#'
tabulateValidation <- function(effects, types, transformations,
                               balancing=NULL, new.data) {

  stopifnot("length of arguments 'effects', 'types' and 'transformations' must be equal" = length(effects)==length(types))
  stopifnot("length of arguments 'effects', 'types' and 'transformations' must be equal" = length(effects)==length(transformations))
  if (!is.null(balancing)) {
    stopifnot("length of arguments 'effects' and 'balancing' must be equal" = length(effects)==length(balancing))
  }
  if (is.null(balancing)) {
    balancing <- rep("", times=length(effects))
  }

  # Register parallel backend
  numCores <- detectCores() - 1  # Use one less than the total number of cores
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  # Parallelized loop
  validations <- foreach(i=1:length(effects), .combine=c,
                         .packages=c("dplyr", "glmnet"),
                         .export=c("validateModelOLS", "validateModelRidgeLASSO",
                                   "metabolites", "aic", "oversample",
                                   "SMOTE", "makeMatrix", "roc", "auc",
                                   "trainRidgeLASSO", "tuneLambda")) %dopar% {
                                     effect <- effects[i]
                                     type <- types[i]
                                     transformation <- transformations[i]
                                     balance <- balancing[i]

                                     # Predict values of the given data with the requested model
                                     if (type == "OLS") {
                                       validation <- validateModelOLS(effect, type, transformation, balance, new.data)
                                     }
                                     if (type == "Ridge" || type == "LASSO") {
                                       validation <- validateModelRidgeLASSO(effect, type, transformation, balance,
                                                                             new.data)
                                     }

                                     return(validation)
                                   }

  # Stop the cluster
  stopCluster(cl)

  formatTypes <- types
  formatTypes[which(formatTypes=="Ridge")] <- "ridge"
  formatTransformations <- tolower(transformations)
  formatTransformations[which(formatTransformations=="inv")] <- "1/x"

  return(bind_cols("effect"=effects, "type"=formatTypes, "transf."=formatTransformations,
                   "balancing"=balancing, "AUC.cv"=validations))
}

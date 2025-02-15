#' This function makes a matrix of a data frame with metabolite measurements and Age. It is usually used for training and predicting with the glmnet-function.
#' Columns for interaction terms can be included if requested
#'
#' @param data A data frame of metabolite values and "Age".
#' @param includeInteraction Logical value for optionally including interaction terms if requested
#' @return A list containing a matrix (mat) and the interaction terms included in the matrix (interactions)
#' @examples
#' metabolites <- c("met1", "met2", "met_110", "met_coll")
#' makeMatrix(df, includeInteraction = TRUE)
#'
makeMatrix <- function(data, includeInteraction=FALSE) {
  relevant <- which(colnames(data) %in% c("Age", metabolites))
  RidgeLASSO_data <- as.matrix(data[, relevant])
  vars0 <- colnames(data)[relevant]
  interactions <- c()

  if (includeInteraction) {
    n <- ncol(RidgeLASSO_data)
    num_interactions <- n * (n - 1) / 2

    # Pre-allocate a matrix for interactions
    interaction_matrix <- matrix(NA, nrow = nrow(RidgeLASSO_data), ncol = num_interactions)

    interaction_index <- 1
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        interaction_matrix[, interaction_index] <- RidgeLASSO_data[,i] * RidgeLASSO_data[,j]
        interactions <- c(interactions, paste(vars0[i], vars0[j], sep="*"))
        interaction_index <- interaction_index + 1
      }
    }

    # Combine original data with interactions
    RidgeLASSO_data <- cbind(RidgeLASSO_data, interaction_matrix)
    colnames(RidgeLASSO_data) <- c(vars0, interactions)
  }

  return(list(mat = RidgeLASSO_data, interactions = interactions))
}

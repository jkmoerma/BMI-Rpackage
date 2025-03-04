% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeMatrix.R
\name{makeMatrix}
\alias{makeMatrix}
\title{This function makes a matrix of a data frame with metabolite measurements and Age. It is usually used for training and predicting with the glmnet-function.
Columns for interaction terms can be included if requested}
\usage{
makeMatrix(data, includeInteraction = FALSE)
}
\arguments{
\item{data}{A data frame of metabolite values and "Age".}

\item{includeInteraction}{Logical value for optionally including interaction terms if requested}
}
\value{
A list containing a matrix (mat) and the interaction terms included in the matrix (interactions)
}
\description{
This function makes a matrix of a data frame with metabolite measurements and Age. It is usually used for training and predicting with the glmnet-function.
Columns for interaction terms can be included if requested
}
\examples{
metabolites <- c("met1", "met2", "met_110", "met_coll")
makeMatrix(df, includeInteraction = TRUE)

}

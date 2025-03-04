% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trainOLS.R
\name{trainOLS}
\alias{trainOLS}
\title{Train a stepwise OLS model from the data with a specification about main/interaction effects and a BMI transformation.}
\usage{
trainOLS(effect, type, transformation, balancing = "", new.data)
}
\arguments{
\item{effect}{"main" or "interaction" dependent on which effects are wished}

\item{type}{Regression type, must be "OLS" for this function}

\item{transformation}{"Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required}

\item{balancing}{Defaults to "", not oversampling the data. Specify "balanced" if oversampling is needed.}

\item{new.data}{Data to build a regression model from}
}
\value{
An OLS regression model
}
\description{
Train a stepwise OLS model from the data with a specification about main/interaction effects and a BMI transformation.
}
\examples{
library(dplyr)
library(tidyr)
library(tibble)

metabolites <- c("met1", "met2", "met_110", "met_coll")

trainOLS(effect="main", type="OLS", transformation="Log", balancing="balanced", new.data=df)

}

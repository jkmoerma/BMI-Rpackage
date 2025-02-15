#' Oversamples a given data set to balance out the "Overweight" and "Obese" obesity classes with respect to the "Normal weight" using SMOTE.
#' The data set must contain the two columns called "ObesityClass" and "met_110",  and some columns of the globally specified "metabolites", "BMI" and/or "Age".
#'
#' @param data_train A data frame. One of the colums must be called "ObesityClass", consisting of values "Normal weight", "Overweight" and "Obese". Some columns must be
#' @return A SMOTE oversampled version of the data set as a data frame. Patient ID and other irrelevant columns are not included
#'
#' @examples
#' metabolites <- c("met_001", "met_002", "met_110", "met_coll")
#'
#' oversampled_data <- oversample(df)
#'
#' table(oversampled_data$ObesityClass)
#' boxplot(BMI~ObesityClass, oversampled_data)
#' boxplot(BMI~ObesityClass, df)
#'
oversample <- function(data_train) {

  stopifnot("met_110 is used for SMOTE interpolation and is not included" = "met_110"%in%names(data_train))

  irrelevant <- which(! names(data_train) %in% c(metabolites, "BMI", "Age", "ObesityClass"))

  for (name in names(data_train)[irrelevant]) {
    data_train[[name]] <- NULL
  }

  overweight_data <- subset(data_train, subset = !ObesityClass=="Obese")
  overweight_data$ObesityClass <- factor(ifelse(overweight_data$ObesityClass == "Overweight",
                                                "Overweight","Normal weight"))
  counts <- table(overweight_data$ObesityClass)
  minority <- min(counts)
  majority <- max(counts)
  overweight_resampled <- SMOTE(ObesityClass ~ met_110, overweight_data,
                                perc.over=100*(majority/minority-1),
                                perc.under=100/(1-minority/majority))

  obese_data <- subset(data_train, subset = !ObesityClass=="Overweight")
  obese_data$ObesityClass <- factor(ifelse(obese_data$ObesityClass == "Obese",
                                           "Obese","Normal weight"))
  counts <- table(obese_data$ObesityClass)
  minority <- min(counts)
  majority <- max(counts)
  obese_resampled <- SMOTE(ObesityClass ~ met_110, obese_data,
                           perc.over=100*(majority/minority-1),
                           perc.under=100/(1-minority/majority))

  data_balanced <- bind_rows(subset(data_train, subset = ObesityClass=="Normal weight"),
                             subset(overweight_resampled, subset = ObesityClass=="Overweight"),
                             subset(obese_resampled, subset = ObesityClass=="Obese"))

  # recalculate metabolite ratios ON LOG SCALE !!
  ratios <- metabolites[which(grepl(pattern="/", metabolites))]
  for (ratio in ratios) {
    mets <- strsplit(ratio, split="/")[[1]]
    data_balanced[[ratio]] <- with(data_balanced, eval(parse(text=paste(mets[1], "-", mets[2]))))
  }

  data_balanced

}

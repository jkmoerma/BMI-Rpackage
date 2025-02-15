#' This function was written to impute and explore missing data. It uses Multivariate Imputation Chained Equations to perform this.
#'
#' @param df A data frame as given by the owner of the data
#' @param makePlot logical, if TRUE, a plot will be generated and stored in file "boxplotMissing.pdf"
#' @return A similar data frame with MICE imputed missing values for met_002, met_010 and met_068.
#'
boxplotMissing <- function(df, makePlot=TRUE) {

  # one missing value of met_010 will be imputed with the ethnicity stratified mean
  # metabolite ratios with met_010 will be recalculated
  met10colnames <- names(df)[which(grepl(names(df), pattern="met_010"))]
  stratified_means_met010 <- df %>% group_by(Race) %>% summarise(exp(mean(log(met_010), na.rm=TRUE)))
  df_complete <- merge(df, stratified_means_met010, by="Race")
  w <- which(is.na(df_complete$met_010))
  df_complete$met_010[w] <- df_complete$`exp(mean(log(met_010), na.rm = TRUE))`[w]
  df_complete$`exp(mean(log(met_010), na.rm = TRUE))` <- NULL
  for (name in met10colnames) {
    df_complete[[name]][w] <- with(df_complete[w,], eval(parse(text=name)))
  }

  data_model <- df_complete
  data_model$Smoking <- df_complete$`Smoking status`
  data_model$`Smoking status` <- NULL
  data_model$`Maternal Age` <- NULL
  data_model$AgeGroup <- NULL
  data_model$ObesityClass <- NULL
  for (met in c("BMI", metabolites)) {
    data_model[[met]] <- log(df_complete[[met]])
  }

  # impute met_002 and met_068 stratified on ethnicity
  log_all <- c()
  for (race in unique(data_model$Race)) {
    data_race <- subset(data_model, subset=Race==race)
    if (nrow(data_race) > length(metabolites)+10) {
      imputed_data <- mice(data_race, method="norm.predict")
    } else {
      imputed_data <- mice(data_race, method="pmm")
    }

    log_race <- complete(imputed_data) %>%
      mutate(met_002C=met_002) %>%
      mutate(met_068C=met_068)
    log_race <- subset(log_race, select=c("ID", "met_002C", "met_068C"))
    log_all <- bind_rows(log_all, log_race)

  }

  df_complete <- merge(df_complete, log_all, by="ID")
  df_complete$met_002 <- exp(df_complete$met_002C)
  df_complete$met_068 <- exp(df_complete$met_068C)

  if (makePlot) {

    missing02 <- wilcox.test(formula = value ~ missing,
                             data=data.frame(value = df_complete$met_002,
                                             missing = df_complete$ID %in% df$ID[which(is.na(df$met_002))]))
    missing68 <- wilcox.test(formula = value ~ missing,
                             data=data.frame(value=df_complete$met_068,
                                             missing = df_complete$ID %in% df$ID[which(is.na(df$met_068))]))

    pdf(file="boxplotMissing.pdf", width=8.27, height=5.845)

    layout(matrix(1:2, ncol=2))

    title <- sprintf("met_002 \n Wilcox. p-val. %.3f \n P(met[miss] > met) = %.3f",
                     missing02$p.value, 1-missing02$statistic/(38*1597))
    boxplot(formula = value ~ missing,
            data=data.frame(value=df_complete$met_002,
                            missing = df_complete$ID %in% df$ID[which(is.na(df$met_002))]),
            xlab="", ylab="log(met_002)", main=title, xaxt="n")
    axis(side=1, at=1:2, labels=c("measured", "missing"))

    title <- sprintf("met_068 \n Wilcox. p-val. %.3f \n P(met[miss] > met) = %.3f",
                     missing68$p.value, 1-missing68$statistic/(49*1586))
    boxplot(formula = value ~ missing,
            data=data.frame(value=df_complete$met_068,
                            missing = df_complete$ID %in% df$ID[which(is.na(df$met_068))]),
            xlab="", ylab="log(met_068)", main=title, xaxt="n")
    axis(side=1, at=1:2, labels=c("measured", "missing"))

    dev.off()

  }

  w <- match(df$ID, df_complete$ID)
  df_complete[w,]

}

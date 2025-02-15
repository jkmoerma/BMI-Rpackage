#' Make a table of patient counts of smokers and non-smokers for every ethnicity
#'
#' @param df A data frame containing columns "Smoking status" and "Race".
#' @return A matrix with all the observed ethnicities as columns, counted "Smoking" and "Non-smoking" patients as entries.
#' @export
#'
tableSmokingStatus <- function(df) {
  count_table <- df %>% group_by(`Smoking status`, Race) %>% summarise(n())
  SmokingStatus <- spread(count_table, key = Race, value = `n()`, fill=0)
  answer <- as.matrix(SmokingStatus[,-1])
  answer <- cbind(answer, "all"=margin.table(answer, margin=1))
  rownames(answer) <- c("Non-smoking", "Smoking")[1+(SmokingStatus$`Smoking status`)]
  answer
}

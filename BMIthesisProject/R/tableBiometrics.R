#' Make a table of patient counts for every ethnicity inside every obesity class
#'
#' @param df A data frame containing columns "ObesityClass" and "Race".
#' @return A matrix with all the observed ethnicities as columns, the different obesity classes and the ethnicity total as rows and the counts as entries.
#' @examples
#'
#' df_normal <- data.frame(ObesityClass=rep("Class1", 40), Race=c(rep("Race1", 15), rep("Race2", 25))
#' df_obese <- data.frame(ObesityClass=rep("Class2", 50), Race=c(rep("Race1", 20), rep("Race2", 30)))
#' df_all <- bind_rows(df_normal, df_obese)
#' tableBiometrics(df_all)
#'
tableBiometrics <- function(df) {
  count_table <- df %>% group_by(ObesityClass, Race) %>% summarise(n())
  biometrics <- spread(count_table, key = Race, value = `n()`)

  count_all <- df %>% group_by(Race) %>% summarise(n())
  biometrics1 <- spread(count_all, key = Race, value = `n()`)

  answer_tibble <- bind_rows(biometrics, biometrics1)
  answer <- as.matrix(answer_tibble[,-1])
  rownames(answer) <- c(as.character(answer_tibble$ObesityClass[-nrow(answer_tibble)]), "total")
  answer
}

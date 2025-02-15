#' Summarizes a patients and their reclassification status to counts and fractions reclassified as syndromatic
#'
#' @param data A data frame containing columns "ObesityClass" and "Reclassified", the latter consisting of values "Healthy" and "Syndromatic"
#' @param race Ethnicity of the patients in the data set
#' @param strat Stratification of the model used for reclassification
#' @return A tibble with columns "Race", "model", "ObesityClass", "pred.: Healthy", "pred. Syndr." and "frac. syndromatic"
#' @examples
#' normals <- data.frame(ObesityClass=rep("Normal weight", times=50),
#'                       Reclassified=c("Healthy", "Syndromatic")[1+rbinom(n=50, size=1, prob=0.2)])
#' overweights <- data.frame(ObesityClass=rep("Overweight", times=50),
#'                           Reclassified=c("Healthy", "Syndromatic")[1+rbinom(n=50, size=1, prob=0.5)])
#' obeses <- data.frame(ObesityClass=rep("Obese", times=50),
#'                           Reclassified=c("Healthy", "Syndromatic")[1+rbinom(n=50, size=1, prob=0.8)])
#' data <- bind_rows(normals, overweights, obeses)
#' countClasses(data, race="Not specified", strat="Unknown")
#'
countClasses <- function(data, race, strat) {

  stopifnot('countClasses uses colums "ObesityClass" and "Reclassified" to summarize' = all(c("ObesityClass", "Reclassified") %in% names(data)))
  stopifnot('Patients must be reclaassified as either "Healthy" or "Syndromatic"' = all(unique(data$Reclassified) %in% c("Healthy", "Syndromatic")))

  classCounts <- data %>%
    group_by(ObesityClass, Reclassified) %>% summarise(n=n())
  classCounts <- pivot_wider(classCounts, names_from="Reclassified",
                             values_from="n")
  classCounts <- classCounts %>%
    mutate(Race=race) %>%
    mutate(model=strat) %>%
    mutate(ObesityClass=as.character(ObesityClass)) %>%
    mutate(`pred.: Healthy`=Healthy) %>%
    mutate(`pred.: Syndr.`=Syndromatic) %>%
    mutate("frac. syndromatic"=calculateFraction(Healthy, Syndromatic))
  classCounts$Healthy <- NULL
  classCounts$Syndromatic <- NULL
  classCounts
}

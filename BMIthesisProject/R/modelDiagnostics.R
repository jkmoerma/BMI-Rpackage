#' Stores a pdf file called paste0("modelDiagnostics", ethnicity, ".pdf") with plots of the residuals and squared residuals as a function of the predicted BMI.
#' The functional form of the predictors significantly deviating from linearity are included.
#' All p-values are calculated from F-tests between an loess-estimator and an intercept estimator of the residuals
#'
#' @param effects a vector of "main" and/or "interaction" effects of formulations the model diagnostics need to be examined of
#' @param types a vector of regression types "OLS", "Ridge" and/or "LASSO" of formulations the model diagnostics need to be examined of
#' @param transformation a vector of BMI transformations "Log" and/or "Inv" of formulations the model diagnostics need to be examined of
#' @param balancing defaults to NULL, or a vector of "" and/or "Balanced" reflecting balancing requests of formulations the model diagnostics need to be examined of
#' @param ethnicity the ethnicity of patients in new.data
#' @param new.data Data to build a regression model from
#' @return NULL
#'
modelDiagnostics <- function(effects, types, transformations, balancing=NULL,
                             ethnicity, model=NULL, new.data) {

  pdf(file=paste0("modelDiagnostics", ethnicity, ".pdf"), width=8.27, height=11.69)

  ethnicities <- rep(ethnicity, times=length(effects))
  if (is.null(balancing)) {balancing <- rep("", times=length(effects))}

  irrelevant <- which(colnames(new.data)%in%c("Race", "ID", "ObesityClass", "BMI", "transBMI", "Smoking"))
  ridge_data <- as.matrix(new.data[,-irrelevant])
  vars0 <- colnames(new.data)[-irrelevant]

  init_ridge_data <- makeMatrix(new.data, includeInteraction=TRUE)
  ridge_data <- init_ridge_data$mat
  interactions <- init_ridge_data$interactions


  for (i in 1:length(effects)) {

    layout(matrix(1:2, ncol=1))

    effect <- effects[i]
    type <- types[i]
    transformation <- transformations[i]
    balance <- balancing[i]
    modeltitle <- paste0(effect, type, transformation, "Model", ethnicity, balance)

    if (is.null(model)) {
      model <- eval(parse(text=modeltitle))
    }

    if (is.null(new.data$predicted)) {
      # predict values of the given data with the requested model
      if (type=="OLS") {
        df_model <- model$df
        valuePreds <- predict(model, newdata=new.data)
      }
      if (type=="Ridge"|type=="LASSO") {
        if (effect=="main") {
          valuePreds <-
            predict(newx=ridge_data[, c(vars0)],
                    object=model)[, "s0"]
        } else {
          valuePreds <-
            predict(newx=ridge_data[, c(vars0, interactions)],
                    object=model)[, "s0"]
        }

      }
    } else {
      if (transformation=="Log") {valuePreds <- log(new.data$predicted)}
      if (transformation=="Inv") {valuePreds <- 1/new.data$predicted}
    }

    if (transformation=="Log") {valueOuts <- new.data$BMI
    SSTOT <- sum((new.data$BMI-mean(new.data$BMI))**2)
    xlabel <- "predicted log(BMI)"}
    if (transformation=="Inv") {valueOuts <- exp(-new.data$BMI)
    SSTOT <- sum((exp(-new.data$BMI)-mean(exp(-new.data$BMI)))**2)
    xlabel <- "predicted 1/BMI"}

    valueResiduals <- valueOuts - valuePreds
    residuals <- data.frame(prediction=valuePreds, residual=valueResiduals,
                            residualSquare=(valueResiduals)**2)


    # plot residuals as function of predicted transformed BMI to assess linearity

    plot(x=valuePreds, y=valueResiduals, main=paste0(modeltitle, ": linearity assessment"),
         xlab=xlabel, ylab="residual")
    abline(h=0, lty="dashed")
    linearity <- loess(formula = residual~prediction, data=residuals)
    linearitySmooth <- predict(linearity, se=TRUE,
                               newdata=data.frame(prediction=sort(valuePreds)))
    lines(sort(valuePreds), linearitySmooth$fit, col="red")
    lines(sort(valuePreds), linearitySmooth$fit - 1.96*linearitySmooth$se.fit,
          col="red", lty="dashed")
    lines(sort(valuePreds), linearitySmooth$fit + 1.96*linearitySmooth$se.fit,
          col="red", lty="dashed")

    linearityFit <- lm(formula=residual~1, data=residuals)
    linearitySSE <- sum(linearityFit$residuals**2)
    linearityDOF <- linearityFit$df.residual

    smoothSSE <- sum(linearity$residuals**2)
    smoothDOF <- linearity$n - linearity$enp

    Fstat <- (linearitySSE-smoothSSE)*smoothDOF/(smoothSSE*(linearityDOF-smoothDOF))
    Rsquare <- 1 - linearitySSE/SSTOT

    pval <- 1-pf(q=Fstat, df1=linearityDOF-smoothDOF, df2=smoothDOF)
    text(x=sum(range(valuePreds))/2, y=max(valueResiduals), pos=1, col="red",
         labels=sprintf("linearity: p = %.4f \n R-square = %.2f", pval, Rsquare))


    # plot squared residuals as function of predicted transformed BMI to assess homoscedasticity

    plot(x=range(valuePreds), y=c(1e-5, max((valueResiduals)**2)),
         col="white", main=paste0(modeltitle, ": homoscedasticity assessment"),
         xlab=xlabel, ylab="squared residual", log="y")
    points(x=valuePreds, y=(valueResiduals)**2)
    abline(h=mean((valueResiduals)**2), lty="dashed", col="blue")
    homoscedacity <- loess(formula = residualSquare~prediction, data=residuals)
    homoscedacitySmooth <- predict(homoscedacity, se=TRUE,
                                   newdata=data.frame(prediction=sort(valuePreds)))
    lines(sort(valuePreds), homoscedacitySmooth$fit, col="red")
    lines(sort(valuePreds), homoscedacitySmooth$fit - 1.96*homoscedacitySmooth$se.fit,
          col="red", lty="dashed")
    lines(sort(valuePreds), homoscedacitySmooth$fit + 1.96*homoscedacitySmooth$se.fit,
          col="red", lty="dashed")

    assess_homoscedasticity <- lm(residualSquare~prediction, data=residuals)

    pval <- summary(assess_homoscedasticity)$coefficients["prediction", "Pr(>|t|)"]
    Qpreds <- quantile(valuePreds, probs=c(0.25, 0.75))
    sigma2 <- predict(homoscedacity, Qpreds)
    factorIncrease <- sigma2[2]/sigma2[1]

    text(x=sum(range(valuePreds))/2, y=max(valueResiduals**2), pos=1, col="red",
         labels=sprintf("homoscedacity: p = %.4f\nQ1(pred) vs Q3(pred) rel. var. incr. = %.2f",
                        pval, factorIncrease))


    # make second page plotting the covariates with a non-linear functional form

    # checking the functional form of the variables entering the model
    nonlinears <- c()
    for (met in c("Age", metabolites)) {
      residuals <- data.frame(feature=new.data[[met]], residual=valueResiduals,
                              residualSquare=(valueResiduals)**2)

      linearity <- loess(formula = residual~feature, data=residuals)
      linearitySmooth <- predict(linearity, se=TRUE,
                                 newdata=data.frame(feature=sort(new.data[[met]])))

      linearityFit <- lm(formula=residual~feature, data=residuals)
      linearitySSE <- sum(linearityFit$residuals**2)
      linearityDOF <- linearityFit$df.residual

      smoothSSE <- sum(linearity$residuals**2)
      smoothDOF <- linearity$n - linearity$enp

      Fstat <- (linearitySSE-smoothSSE)*smoothDOF/(smoothSSE*(linearityDOF-smoothDOF))

      pval <- 1-pf(q=Fstat, df1=linearityDOF-smoothDOF, df2=smoothDOF)

      if (pval<0.05) {nonlinears <- c(nonlinears, met)}

    }

    if (length(nonlinears) > 0) {

      layout(matrix(1:6, ncol=2, byrow=TRUE))

      variables <- c("Age", metabolites)

      reduced_vars <- variables[-which(variables %in% nonlinears)]
      reduced_model <- lm(formula = eval(parse(text=paste0("BMI~`", paste(reduced_vars, collapse="`+`"), "`"))),
                          data=new.data)
      for (met in nonlinears) {

        residuals <- data.frame(feature=new.data[[met]], residual=valueResiduals,
                                residualSquare=(valueResiduals)**2)

        plot(x=new.data[[met]], y=valueResiduals, main=paste0("model residuals vs. ", met),
             xlab=met, ylab="residual")
        abline(h=0, lty="dashed")
        linearity <- loess(formula = residual~feature, data=residuals)
        linearitySmooth <- predict(linearity, se=TRUE,
                                   newdata=data.frame(feature=sort(new.data[[met]])))
        lines(sort(new.data[[met]]), linearitySmooth$fit, col="red")
        lines(sort(new.data[[met]]), linearitySmooth$fit - 1.96*linearitySmooth$se.fit,
              col="red", lty="dashed")
        lines(sort(new.data[[met]]), linearitySmooth$fit + 1.96*linearitySmooth$se.fit,
              col="red", lty="dashed")

        linearityFit <- lm(formula=residual~feature, data=residuals)
        linearitySSE <- sum(linearityFit$residuals**2)
        linearityDOF <- linearityFit$df.residual

        smoothSSE <- sum(linearity$residuals**2)
        smoothDOF <- linearity$n - linearity$enp

        Fstat <- (linearitySSE-smoothSSE)*smoothDOF/(smoothSSE*(linearityDOF-smoothDOF))

        pval <- 1-pf(q=Fstat, df1=linearityDOF-smoothDOF, df2=smoothDOF)
        text(x=sum(range(new.data[[met]]))/2, y=max(valueResiduals), pos=1,
             labels=sprintf("linearity: p = %.2e", pval), col="red")

        modelVar <- lm(formula = eval(parse(text=paste(met,"~",paste(reduced_vars, collapse="+")))),
                       data=new.data)

        plot(x=modelVar$residuals, y=reduced_model$residuals,
             xlab=paste0(met, "-pred(", met, ")"), ylab="residual",
             main=paste0("Desired functional form ", met))
        abline(h=0, lty="dashed")
        residuals <- data.frame(feature=modelVar$residuals, residual=reduced_model$residuals)
        linearity <- loess(formula = residual~feature, data=residuals)
        linearitySmooth <- predict(linearity, se=TRUE,
                                   newdata=data.frame(feature=sort(modelVar$residuals)))
        lines(sort(modelVar$residuals), linearitySmooth$fit, col="red")
        lines(sort(modelVar$residuals), linearitySmooth$fit - 1.96*linearitySmooth$se.fit,
              col="red", lty="dashed")
        lines(sort(modelVar$residuals), linearitySmooth$fit + 1.96*linearitySmooth$se.fit,
              col="red", lty="dashed")

      }

    }

  }

  dev.off()

}

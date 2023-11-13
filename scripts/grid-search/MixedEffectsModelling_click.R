## Author: Surabhi S Nath
## Description: This script implements survival analysis (Cox regression) models at click level.
## on grid-search exploration data.
## Helper functions and code for plotting in utils/Utils.R
## model analysis table written to model_fits/

# Imports
{
    library(lme4)
    library(nlme)
    library("afex")
    library(lmerTest)
    library(sjPlot)
    library(sjmisc)
    library(ggplot2)
    library(insight)
    library(lattice)
    library(fitdistrplus)
    library(corrplot)
    library(MASS)
    library(merTools)
    library(car)
    library("survival")
    library("survminer")
    library(jtools)
    library(interactions)

    source("utils/Utils.R")
}

# Setup
{
    # Read data
    data = read.csv("utils/click_data.csv")

    # Scale all variables in the data
    my_scale <- function(x) {
        as.numeric(scale(x))
    }

    data$Subject <- factor(data$Subject)
    data$click_id <- as.numeric(data$click_id)
    data$current_LSC_sq <- my_scale(data$current_LSC ^ 2)
    data$current_LSC <- my_scale(data$current_LSC)
    data$current_intricacy_sq <- my_scale(data$current_intricacy ^ 2)
    data$current_intricacy <- my_scale(data$current_intricacy)
    data$underlying_LSC_sq <- my_scale(data$underlying_LSC ^ 2)
    data$underlying_LSC <- my_scale(data$underlying_LSC)
    data$underlying_intricacy_sq <- my_scale(data$underlying_intricacy ^ 2)
    data$underlying_intricacy <- my_scale(data$underlying_intricacy)
    data$current_change_in_LSC <- my_scale(data$current_change_in_LSC)
    data$current_change_in_intricacy <- my_scale(data$current_change_in_intricacy)
    data$avg_change_in_LSC <- my_scale(data$avg_change_in_LSC)
    data$avg_change_in_intricacy <- my_scale(data$avg_change_in_intricacy)
}


{
  num_folds <- 3

  folds <- data %>%
  group_by(Subject) %>%
  modelr::crossv_kfold(k = num_folds)

  temp_train_1 <- folds$train$`1`
  data_train_fold1 <- data[temp_train_1$idx, ]

  temp_test_1 <- folds$test$`1`
  data_test_fold1 <- data[temp_test_1$idx, ]

  temp_train_2 <- folds$train$`2`
  data_train_fold2 <- data[temp_train_2$idx, ]

  temp_test_2 <- folds$test$`2`
  data_test_fold2 <- data[temp_test_2$idx, ]

  temp_train_3 <- folds$train$`3`
  data_train_fold3 <- data[temp_train_3$idx, ]

  temp_test_3 <- folds$test$`3`
  data_test_fold3 <- data[temp_test_3$idx, ]

  data_train_folds = list(data_train_fold1, data_train_fold2, data_train_fold3)
  data_test_folds = list(data_test_fold1, data_test_fold2, data_test_fold3)
}

# Cox regression models
models <- list(
    "Subject",
    "underlying_LSC + underlying_intricacy",
    "current_LSC + current_intricacy",
    "current_change_in_LSC + current_change_in_intricacy",
    "underlying_LSC + underlying_intricacy + current_LSC + current_intricacy",
    "underlying_LSC + underlying_intricacy + current_LSC + current_intricacy + current_change_in_LSC + current_change_in_intricacy",
    "underlying_LSC_sq + underlying_intricacy_sq + current_LSC_sq + current_intricacy_sq",
    "underlying_LSC * current_LSC",
    "underlying_intricacy * current_intricacy",
    "underlying_LSC * current_LSC + underlying_intricacy * current_intricacy",
    "underlying_LSC * current_LSC + underlying_intricacy * current_intricacy + Subject"
)

# Save all results to model_fits/Table_2_click_level.csv
{
    df <- data.frame(matrix(ncol = 8, nrow = 0, dimnames =
    list(NULL, c("Id", "model", "AIC", "BIC",
        "Concordance train", "Concordance test", "log rank p", "likelihood ratio p"))))

    id <- 0
    for(formula in models) 
    {
        id <- id + 1
        fullformula <- paste("Surv(click_id, move_on) ~ ", formula)
        formula <- as.formula(fullformula)
        
        # calcultate the following metrics
        AICs <- numeric(num_folds)
        BICs <- numeric(num_folds)
        concordances_train <-numeric(num_folds)
        concordances_test <-numeric(num_folds)
        logrankp <- numeric(num_folds)
        likratiop <- numeric(num_folds)

        for(x in 1:num_folds)
        {
            model <- coxph(formula, data = data_train_folds[[x]])
            
            # calculate metrics and store in arrays
            modsum = summary(model)
            logrankp[[x]] <- modsum$logtest["pvalue"]
            likratiop[[x]] <- modsum$sctest["pvalue"]

            pred_train <- predict(model, newdata = data_train_folds[[x]], type = 'risk')
            concordances_train[[x]] <- survConcordance(Surv(click_id, move_on) ~ pred_train, data = data_train_folds[[x]])$concordance
            pred_test <- predict(model, newdata = data_test_folds[[x]], type = 'risk')
            concordances_test[[x]] <- survConcordance(Surv(click_id, move_on) ~ pred_test, data = data_test_folds[[x]])$concordance

            AICs[x] <- AIC(model)
            BICs[x] <- BIC(model)
        }

        df[nrow(df) + 1, ] <- c(id, noquote(fullformula), c(mean(AICs), mean(BICs), mean(concordances_train), mean(concordances_test), max(logrankp), max(likratiop)))
    }

    write.csv(df, "model_fits/Table_2_click_level.csv", row.names = FALSE)
}


# Save fixed effects to model_fits/ - Table 3
{
  bestformula <- as.formula("Surv(click_id, move_on) ~  underlying_LSC * current_LSC + underlying_intricacy * current_intricacy")
  model <- coxph(bestformula, data = data_train_folds[[1]])
  write.csv(coef(model), file = "model_fits/Table_3_coeff_move_on.csv")
}

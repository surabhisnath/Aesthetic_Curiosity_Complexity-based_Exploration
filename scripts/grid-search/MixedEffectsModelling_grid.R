## Author: Surabhi S Nath
## Description: This script implements mixed effects models at grid level.
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
  library(jtools)
  library(interactions)
  library(modelr)

  source("utils/Utils.R")

  # Set seed 
  set.seed(20) # Seed set randomly for reproducibility
}

# Setup
{
  # Read data
  data <- read.csv("utils/grid_data.csv")

  # Scale all variables in the data
  my_scale <- function(x) {
    as.numeric(scale(x))
  }

  data$Subject <- factor(data$Subject)
  data$grid_id <- factor(data$grid_id)
  data$pattern <- factor(data$pattern)
  data$underlying_LSC_sq <- my_scale(data$underlying_LSC ^ 2)
  data$underlying_LSC <- my_scale(data$underlying_LSC)
  data$underlying_intricacy_sq <- my_scale(data$underlying_intricacy ^ 2)
  data$underlying_intricacy <- my_scale(data$underlying_intricacy)
  data$final_LSC_sq <- my_scale(data$final_LSC ^ 2)
  data$final_LSC <- my_scale(data$final_LSC)
  data$final_intricacy_sq <- my_scale(data$final_intricacy ^ 2)
  data$final_intricacy <- my_scale(data$final_intricacy)
  data$avg_LSC <- my_scale(data$avg_LSC)
  data$avg_intricacy <- my_scale(data$avg_intricacy)
  data$final_change_in_LSC <- my_scale(data$final_change_in_LSC)
  data$final_change_in_intricacy <- my_scale(data$final_change_in_intricacy)
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
}

# Mixed Effects Models
models <- list(
  "1 + (1 | Subject)",
  "underlying_LSC + underlying_intricacy + (1 | Subject)",
  "final_LSC + final_intricacy + (1 | Subject)",
  "final_change_in_LSC + final_change_in_intricacy + (1 | Subject)",
  "underlying_LSC + underlying_intricacy + final_LSC + final_intricacy + (1 | Subject)",
  "underlying_LSC_sq + underlying_intricacy_sq + final_LSC_sq + final_intricacy_sq + (1 | Subject)",
  "underlying_LSC + underlying_intricacy + final_LSC + final_intricacy + (underlying_LSC + underlying_intricacy| Subject)",
  "underlying_LSC + underlying_intricacy + final_LSC + final_intricacy + (final_LSC + final_intricacy| Subject)",
  "avg_LSC + avg_intricacy + (1 | Subject)",
  "underlying_LSC + underlying_intricacy + final_LSC + final_intricacy + final_change_in_LSC + final_change_in_intricacy + (1 | Subject)",
  "underlying_LSC * underlying_intricacy + (1 | Subject)",
  "final_LSC * final_intricacy + (1 | Subject)",
  "underlying_LSC * final_LSC + (1 | Subject)",
  "underlying_intricacy * final_intricacy + (1 | Subject)",
  "underlying_LSC * final_LSC + underlying_intricacy * final_intricacy + (1 | Subject)",
  "underlying_LSC_sq * underlying_intricacy_sq + final_LSC_sq * final_intricacy_sq + (1 | Subject)",
  "underlying_LSC * final_LSC * underlying_intricacy * final_intricacy + (1 | Subject)",
  "underlying_LSC * final_change_in_LSC + underlying_intricacy * final_change_in_intricacy + (1 | Subject)"
)

# Save all results to model_fits/Table_4_grid_level.csv
{
  df <- data.frame(matrix(ncol = 13, nrow = 0, dimnames =
  list(NULL, c("Id", "model", "AIC", "BIC", "AIC/BIC Var",
  "Rsq train mean", "Rsq train var", "Rsq test mean", "Rsq test var",
  "RMSE train mean", "RMSE train var", "RMSE test mean", "RMSE test var"))))

  id <- 0
  for (formula in models) {
    id <- id + 1
    fullformula <- paste("num_clicks ~", formula)
    f1 <- lmer(fullformula, data = data_train_fold1,
    control = lmerControl(optimizer = "bobyqa"))
    f2 <- lmer(fullformula, data = data_train_fold2,
    control = lmerControl(optimizer = "bobyqa"))
    f3 <- lmer(fullformula, data = data_train_fold3,
    control = lmerControl(optimizer = "bobyqa"))

    # Analyse model fit
    metrics <- modelanalysis(num_folds,
    list(f1, f2, f3), list(data_train_fold1, data_train_fold2,
    data_train_fold3), list(data_test_fold1, data_test_fold2, data_test_fold3),
    FALSE, FALSE, fullformula) # set second last param to TRUE for printing
    df[nrow(df) + 1, ] <- c(id, noquote(fullformula), metrics)
  }

  write.csv(df, "model_fits/Table_4_grid_level.csv", row.names = FALSE)
}

# Evaluate performance of best model and make plots
# for train, test and random effects
# Plots saved to ./plots/
{
  bestformula <- "num_clicks ~ underlying_LSC * final_LSC + underlying_intricacy * final_intricacy + (1 | Subject)"
  f1 <- lmer(bestformula, data = data_train_fold1,
  control = lmerControl(optimizer = "bobyqa"))
  f2 <- lmer(bestformula, data = data_train_fold2,
  control = lmerControl(optimizer = "bobyqa"))
  f3 <- lmer(bestformula, data = data_train_fold3,
  control = lmerControl(optimizer = "bobyqa"))
  f <- lmer(bestformula, data = data_train_fold1,
  control = lmerControl(optimizer = "bobyqa"))

  # Plots data vs predictions on train and test data
  metrics <- modelanalysis(num_folds,
  list(f1, f2, f3), list(data_train_fold1, data_train_fold2,
  data_train_fold3), list(data_test_fold1, data_test_fold2, data_test_fold3),
  FALSE, TRUE)

  # Plots random effects
  ggCaterpillar(ranef(f, condVar = TRUE))
  ggsave("plots/random_effects_num_clicks.pdf")
}

# Save fixed effects to model_fits/ - Table 5
{
  write.csv(fixef(f), file = "model_fits/Table_5_coeff_num_clicks.csv")
}
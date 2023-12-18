## Author: Surabhi S Nath
## Description: This script implements mixed effects models of complexity
## ratings on 27x27, semi-revealed patterns.
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

  source("utils/Utils.R")

  # Set seed
  set.seed(20) # Seed set randomly for reproducibility
}

# Setup
{
  # Read data
  data <- read.csv("utils/data.csv")
  # remove subjects who failed all attention checks
  all_attention_failed <- c("fafl6odp", "s7hnc00q", "p60mbj11")
  data <- data[! data$subject %in% all_attention_failed]

  # Scale all variables in the data
  my_scale <- function(x) {
    as.numeric(scale(x))
  }

  data$rtime <- data$rtime / 1000
  data$subject <- factor(data$subject)
  data$set <- factor(data$set, levels = c(1, 2, 4))
  data$pattern <- factor(data$pattern)
  data$num_colours <- as.numeric(data$num_colours)
  data$num_revealed <- as.numeric(data$num_revealed)
  data$visible_complexity_rating <- my_scale(data$visible_complexity_rating)
  data$imagined_complexity_rating <- my_scale(data$imagined_complexity_rating)
  data$LSC_sq <- my_scale(data$LSC ^ 2)
  data$intricacy_sq <- my_scale(data$intricacy ^ 2)
  data$underlying_LSC_sq <- my_scale(data$underlying_LSC ^ 2)
  data$underlying_intricacy_sq <- my_scale(data$underlying_intricacy ^ 2)
  data$LSC <- my_scale(data$LSC)
  data$intricacy <- my_scale(data$intricacy)
  data$underlying_LSC <- my_scale(data$underlying_LSC)
  data$underlying_intricacy <- my_scale(data$underlying_intricacy)
}

# K-fold stratified sampling
{
  num_folds <- 3

  folds <- data %>%
    group_by(subject) %>%
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

# Modelling Mixed Models
models <- list(
  "1 + (1 | subject)",
  "LSC + (1 | subject)",
  "intricacy + (1 | subject)",
  "imagined_complexity_rating + (1 | subject)",
  "LSC + intricacy + (1 | subject)",
  "LSC + intricacy + imagined_complexity_rating + (1 | subject)",
  "underlying_LSC + underlying_intricacy + (1 | subject)",
  "LSC + intricacy + num_revealed + num_colours + (1 | subject)",
  "LSC + intricacy + underlying_LSC + underlying_intricacy + (1 | subject)",
  "LSC * intricacy + underlying_LSC * underlying_intricacy + (1 | subject)",
  "LSC * underlying_LSC + intricacy * underlying_intricacy + (1 | subject)",
  "LSC_sq + intricacy_sq + underlying_LSC_sq + underlying_intricacy_sq + (1 | subject)",
  "LSC + intricacy + underlying_LSC + underlying_intricacy + ((LSC + underlying_LSC) | subject)",
  "LSC * underlying_LSC + intricacy * underlying_intricacy + ((LSC + intricacy) | subject)",
  "LSC * underlying_LSC + intricacy * underlying_intricacy + ((LSC + underlying_LSC) | subject)",
  "LSC * underlying_LSC + intricacy * underlying_intricacy + ((intricacy + underlying_intricacy) | subject)"
)

# Save all results to model_fits/
{
  df <- data.frame(matrix(ncol = 13, nrow = 0,
                          dimnames =
                            list(NULL,
                                 c("Id", "model", "AIC", "BIC", "AIC/BIC Var",
                                   "Rsq train mean", "Rsq train var",
                                   "Rsq test mean", "Rsq test var",
                                   "RMSE train mean", "RMSE train var",
                                   "RMSE test mean", "RMSE test var"))))

  id <- 0
  for (formula in models) {
    id <- id + 1
    fullformula <- paste("visible_complexity_rating ~", formula)
    f1 <- lmer(fullformula,
               data = data_train_fold1,
               control = lmerControl(optimizer = "bobyqa"))
    f2 <- lmer(fullformula,
               data = data_train_fold2,
               control = lmerControl(optimizer = "bobyqa"))
    f3 <- lmer(fullformula,
               data = data_train_fold3,
               control = lmerControl(optimizer = "bobyqa"))

    # Analyse model fit
    metrics <-
      modelanalysis(num_folds,
                    list(f1, f2, f3),
                    list(data_train_fold1, data_train_fold2, data_train_fold3),
                    list(data_test_fold1, data_test_fold2, data_test_fold3),
                    FALSE, FALSE, fullformula) # set second last param to TRUE for printing
    df[nrow(df) + 1, ] <- c(id, noquote(fullformula), metrics)
  }

  write.csv(df, "model_fits/models_complexity-extended.csv", row.names = FALSE)
}

# Evaluate performance of best model and make plots
# for train, test and random effects
# Plots saved to ./plots/
{
  bestformula <- "visible_complexity_rating ~ LSC + intricacy + underlying_LSC + underlying_intricacy + ((LSC + underlying_LSC) | subject)"
  f1 <- 
    lmer(bestformula, data = data_train_fold1,
         control = lmerControl(optimizer = "bobyqa"))
  f2 <-
    lmer(bestformula, data = data_train_fold2,
         control = lmerControl(optimizer = "bobyqa"))
  f3 <-
    lmer(bestformula, data = data_train_fold3,
         control = lmerControl(optimizer = "bobyqa"))
  f <- 
    lmer(bestformula, data = data_train_fold1,
         control = lmerControl(optimizer = "bobyqa"))

  # Plots data vs predictions on train and test data
  metrics <-
    modelanalysis(num_folds,
                  list(f1, f2, f3),
                  list(data_train_fold1, data_train_fold2, data_train_fold3),
                  list(data_test_fold1, data_test_fold2, data_test_fold3),
                  FALSE, TRUE)

  # Plots random effects
  ggCaterpillar(ranef(f, condVar = TRUE))
  ggsave("plots/random_effects_complexity-extended.pdf")

  # Save fixed effects to model_fits/
  write.csv(fixef(f), file = "model_fits/complexity-extended_coeff.csv")
}

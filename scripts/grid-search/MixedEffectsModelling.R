## Author: Surabhi S Nath
## Description: This script implements mixed effects models.
## on grid-search exploration data.
## Helper functions and code for plotting in utils/Utils.R
## model analysis table written to model_fits/
## model plots written to plots/

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

  data$uLSC <- my_scale(data$underlying_LSC)
  data$uLSCsq <- data$uLSC ^ 2

  data$uInt <- my_scale(data$underlying_intricacy)
  data$uIntsq <- data$uInt ^ 2

  data$vLSC <- my_scale(data$final_LSC)
  data$vLSCsq <- data$vLSC ^ 2

  data$vInt <- my_scale(data$final_intricacy)
  data$vIntsq <- data$vInt ^ 2

  data$avg_LSC <- my_scale(data$avg_LSC)
  data$avg_intricacy <- my_scale(data$avg_intricacy)
  data$final_change_in_LSC <- my_scale(data$final_change_in_LSC)
  data$final_change_in_intricacy <- my_scale(data$final_change_in_intricacy)
  data$avg_change_in_LSC <- my_scale(data$avg_change_in_LSC)
  data$avg_change_in_intricacy <- my_scale(data$avg_change_in_intricacy)
}

# Split the data into 3 stratified folds
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

  "uLSC + uInt + (1 | Subject)",
  "vLSC + vInt + (1 | Subject)",
  "uLSC + vLSC + (1 | Subject)",
  "uInt + vInt + (1 | Subject)",

  "uLSC + uInt + vLSC + vInt + (1 | Subject)",

  "uLSC * vLSC + (1 | Subject)",
  "uInt * vInt + (1 | Subject)",
  "uLSC * uInt + (1 | Subject)",
  "vLSC * vInt + (1 | Subject)",

  "uLSC * vLSC + uInt * vInt + (1 | Subject)",  # best
  "uLSC * uInt + vLSC * vInt + (1 | Subject)",
  "uLSC + uInt + vLSC + vInt + uLSC:vInt + vLSC:uInt + (1 | Subject)",
  "uLSC + vLSC + uInt + vInt + uLSC:vLSC + uInt:vInt + uLSC:uInt + vLSC:vInt + uLSC:vInt + vLSC:uInt + (1 | Subject)",


  # Supplementary analysis

  # change in complexity
  "final_change_in_LSC + final_change_in_intricacy + (1 | Subject)",
  "uLSC + uInt + vLSC + vInt + final_change_in_LSC + final_change_in_intricacy + (1 | Subject)",

  # quadratic effects
  "uLSCsq + uIntsq + vLSCsq + vIntsq + (1 | Subject)",
  "uLSCsq + vLSCsq + uIntsq + vIntsq + uLSC:vLSC + uInt:vInt + (1 | Subject)",
  "uLSCsq + uIntsq + vLSCsq + vIntsq + uLSC:vLSC + uInt:vInt + (1 | Subject)",
  "(uLSCsq + vLSCsq + uLSC:vLSC) + (uIntsq + vIntsq + uInt:vInt) + (uLSC + vLSC):(uInt + vInt) + (1 | Subject)",
  "uLSC + vLSC + uInt + vInt + uLSCsq + uIntsq + vLSCsq + vIntsq + (1 | Subject)",
  "uLSC + vLSC + uInt + vInt + uLSCsq + uIntsq + vLSCsq + vIntsq + uLSC:vLSC + uInt:vInt + (1 | Subject)",
  "uLSC + vLSC + uInt + vInt + uLSCsq + uIntsq + vLSCsq + vIntsq + uLSC:vLSC + uInt:vInt + uLSC:uInt + vLSC:vInt + (1 | Subject)",
  "uLSC + vLSC + uInt + vInt + uLSCsq + uIntsq + vLSCsq + vIntsq + uLSC:vLSC + uInt:vInt + uLSC:uInt + vLSC:vInt + uLSC:vInt + uInt:vLSC + (1 | Subject)"
)

# Save all results to model_fits/Table_3_mixedeffects.csv
{
  df <-
    data.frame(matrix(ncol = 13, nrow = 0,
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
    fullformula <- paste("num_clicks ~", formula)
    f1 <-
      lmer(fullformula,
           data = data_train_fold1,
           control = lmerControl(optimizer = "bobyqa"))
    f2 <-
      lmer(fullformula, 
           data = data_train_fold2,
           control = lmerControl(optimizer = "bobyqa"))
    f3 <- 
      lmer(fullformula,
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

  write.csv(df, "model_fits/Table_3_mixedeffects.csv", row.names = FALSE)
}

# Evaluate performance of best model and make plots
# for train, test and random effects
# Plots saved to ./plots/
{
  bestformula <- "num_clicks ~ uLSC + vLSC + uInt + vInt + uLSC:vLSC + uInt:vInt + (1 | Subject)"

  f <- lmer(bestformula,
            data = data,
            control = lmerControl(optimizer = "nloptwrap"))

  # make and save model vs predictions of train and test - Figure 3a,b
  make_plots(f, data_train_fold1, data_test_fold1)

  # Plot and save random effects
  ggCaterpillar(ranef(f, condVar = TRUE))
  ggsave("plots/random_effects_num_clicks.pdf")
}

# Save fixed effects to model_fits/ - Table 4
{
  model_summary <- summary(f)
  coefficients <- fixef(f)
  standard_errors <- sqrt(diag(vcov(f)))
  variable_names <- rownames(summary(f)$coefficients)
  p_values <- coef(summary(f))[, "Pr(>|t|)"]
  signif_levels <-
    ifelse(p_values < 0.001, "***",
           ifelse(p_values < 0.01, "**",
                  ifelse(p_values < 0.05, "*",
                         ifelse(p_values < 0.1, ".", "NS"))))
  results_df <-
    data.frame(Variable = variable_names,
               Coefficients = coefficients,
               StdError = standard_errors,
               PValue = p_values,
               Significance = signif_levels)
  write.csv(results_df,
            file = "model_fits/Table_4_coeff_num_clicks.csv",
            row.names = FALSE)
}

# plot interactions - Figure 3c,d
{
  p <- interact_plot(f,
    pred = vLSC, modx = uLSC, interval = TRUE,
    x.label = "vLSC", y.label = "Number of Clicks",
    legend.main = "uLSC", colors = "seagreen"
  ) + theme(
    axis.title = element_text(family = "serif", size = 44),
    axis.text = element_text(family = "serif", size = 26),
    legend.text = element_text(family = "serif", size = 30),
    legend.title = element_text(family = "serif", size = 40),
    strip.text = element_text(family = "serif")
  )

  # Figure 3c
  ggsave(filename = "plots/uLSC_vLSC_interaction.pdf",
         plot = p, width = 10, height = 10, units = "in")

  p <- interact_plot(f,
    pred = vInt, modx = uInt, interval = TRUE,
    x.label = "vInt", y.label = "Number of Clicks",
    legend.main = "uInt", colors = "seagreen"
  ) + theme(
    axis.title = element_text(family = "serif", size = 44),
    axis.text = element_text(family = "serif", size = 26),
    legend.text = element_text(family = "serif", size = 30),
    legend.title = element_text(family = "serif", size = 44),
    strip.text = element_text(family = "serif")
  )

  # Figure 3d
  ggsave(filename = "plots/uInt_vInt_interaction.pdf",
         plot = p, width = 10, height = 10, units = "in")
}

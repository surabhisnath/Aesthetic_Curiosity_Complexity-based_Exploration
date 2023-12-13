library(lme4)
library(ggplot2)
library(interactions)
library(dplyr)

# Catterpillar plot for random effects
# Code adapted from: https://stackoverflow.com/questions/13847936/plot-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot-how-to-mak, answer by caracal
ggCaterpillar <- function(re, QQ = FALSE, likeDotplot = TRUE, detailedFacetLabs = TRUE) {
    f <- function(x, nm = "ranef plot") {
    pv <- attr(x, "postVar")
    cols <- 1:(dim(pv)[1])
    se <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
    ord <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each = nrow(x))
    pDf <- data.frame(
        y = unlist(x)[ord],
        ci = 1.96 * se[ord],
        nQQ = rep(stats::qnorm(stats::ppoints(nrow(x))), ncol(x)),
        ID = factor(rep(rownames(x), ncol(x))[ord], levels = rownames(x)[ord]),
        ind = gl(ncol(x), nrow(x), labels = names(x))
    )

    if (detailedFacetLabs) {
        pDf$ind <- ifelse(grepl("(Intercept)", pDf$ind), "intercept adjustment", paste0("slope adj: ", pDf$ind))
    }

    if (QQ) { ## normal QQ-plot
        p <- ggplot(pDf, aes_string(x = "nQQ", y = "y"))
        p <- p + facet_wrap(~ind, scales = "free")
        p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
    } else { ## caterpillar dotplot
        p <- ggplot(pDf, aes_string(x = "ID", y = "y")) +
        coord_flip()
        if (likeDotplot) { ## imitate dotplot() -> same scales for random effects
        p <- p + facet_wrap(~ind)
        } else { ## different scales for random effects
        p <- p + facet_grid(ind ~ ., scales = "free_y")
        }
        p <- p + xlab(nm) + ylab("Random effects")
        scale <- 12 - log(length(levels(pDf$ID)), 2)
        p <- p + theme(axis.text.y = element_text(size = scale))
    }

    p <- p + theme(plot.title = element_text(family = "serif", size = 16),
        axis.title = element_text(family = "serif", size = 26),
        axis.text.x = element_text(family = "serif", size = 18),
        axis.text.y = element_text(family = "serif", size = 18)
        )
    p <- p + geom_hline(yintercept = 0, lwd = I(7 / 12), colour = I(grDevices::hsv(0 / 12, 7 / 12, 7 / 12)), alpha = I(5 / 12))
    p <- p + geom_errorbar(aes_string(ymin = "y - ci", ymax = "y + ci"), width = 0, colour = "black")
    p <- p + geom_point(aes())
    return(p)
  }

    lapply(seq_along(re), function(y, n, i) {
    f(y[[i]], n[[i]])
    }, y = re, n = names(re)) # adds plot names
}

# Analysis function
modelanalysis <- function(num_folds, models, data_train_folds, data_test_folds, to_print, to_plot, fullformula) {
    # This function performs all analyses on the model and returns all performance metrics
    # Arguements:
    # num_folds: number of CV folds
    # models - array of lmer models on each fold
    # data_train_folds - ground truth train data for each fold
    # data_test_folds - ground truth test data for each fold
    # to_print - logical - indicates if metrics need to be printed (TRUE means yes)
    # to_plot - logical - indicates if plots should be made (TRUE means yes)

    VIFs <- numeric(num_folds)
    AICs <- numeric(num_folds)
    BICs <- numeric(num_folds)
    rsqtrains <- numeric(num_folds)
    rsqtests <- numeric(num_folds)
    RMSEtrains <- numeric(num_folds)
    RMSEtests <- numeric(num_folds)

    for (x in 1:num_folds) {

    AICs[x] <- AIC(models[[x]])
    BICs[x] <- BIC(models[[x]])

    rsqtrains[x] <- cor(predict(models[[x]], data_train_folds[[x]]), data_train_folds[[x]]$num_clicks)^2
    rsqtests[x] <- cor(predict(models[[x]], data_test_folds[[x]]), data_test_folds[[x]]$num_clicks)^2
    RMSEtrains[x] <- sqrt(mean(residuals(models[[x]])^2))
    RMSEtests[x] <- sqrt(mean((predict(models[[x]], data_test_folds[[x]]) - data_test_folds[[x]]$num_clicks)^2))
    }

    mean_aic <- mean(AICs)
    mean_bic <- mean(BICs)
    var_aic_bic <- var(AICs)
    mean_rsq_train <- mean(rsqtrains)
    var_rsq_train <- var(rsqtrains)
    mean_rsq_test <- mean(rsqtests)
    var_rsq_test <- var(rsqtests)
    mean_rmse_train <- mean(RMSEtrains)
    var_rmse_train <- var(RMSEtrains)
    mean_rmse_test <- mean(RMSEtests)
    var_rmse_test <- var(RMSEtests)

    if (to_print) {

        print(noquote(fullformula))

        print(noquote(paste("Mean AIC =", mean_aic)))
        print(noquote(paste("Mean BIC =", mean_bic)))
        print(noquote(paste("Var AIC, BIC =", var_aic_bic)))

        print(noquote(paste("Mean R^2 train =", mean_rsq_train)))
        print(noquote(paste("Var R^2 train =", var_rsq_train)))
        print(noquote(paste("Mean R^2 test =", mean_rsq_test)))
        print(noquote(paste("Var R^2 test =", var_rsq_test)))

        print(noquote(paste("Mean train RMSE =", mean_rmse_train)))
        print(noquote(paste("Var train RMSE =", var_rmse_train)))
        print(noquote(paste("Mean test RMSE =", mean_rmse_test)))
        print(noquote(paste("Var test RMSE =", var_rmse_test)))
    }

    if (to_plot) # Make plot if plot is TRUE
    {
        make_plots(models[[1]], data_train_folds[[1]], data_test_folds[[1]])
    }

    metrics1 <- c(mean_aic, mean_bic, var_aic_bic)
    metrics1 <- round(metrics1, digits = 1)
    metrics2 <- c(mean_rsq_train, var_rsq_train, mean_rsq_test, var_rsq_test, mean_rmse_train, var_rmse_train, mean_rmse_test, var_rmse_test)
    metrics2 <- signif(metrics2, digits = 4)

    # returns all metrics
    return(c(metrics1, metrics2))
}

make_plots <- function(model, data_train_fold, data_test_fold) {
    # Make and save plots
    # Save plot of num_clicks vs predictions (on both train and test data)

    pdf(file = "plots/num_clicks_train.pdf", width = 10, height = 10, family = "Times")
    par(mar = c(5, 6, 4, 1) + .1)
    p1 <- plot(predict(model, data_train_fold), data_train_fold$num_clicks, xlab = "Predictions on training data", ylab = "Number of Clicks", cex.axis = 3, cex.lab = 3, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3), pch = 16, cex = 1.5, axes = FALSE)
    x_range <- range(predict(model, data_train_fold))
    axis(1, at = seq(from = floor(x_range[1]/20)*20, to = ceiling(x_range[2]/20)*20, by = 20), cex.axis = 2)
    # axis(1, cex.axis = 2)
    axis(2, cex.axis = 2)
    p1 + theme(
    plot.title = element_text(family = "serif", size = 18),
    axis.title = element_text(family = "serif", size = 50),
    axis.text = element_text(family = "serif", size = 30)
    )
    abline(a = 0, b = 1, col = "blue", lwd = 5, lty = 2)
    dev.off()

    pdf(file = "plots/num_clicks_test.pdf", width = 10, height = 10, family = "Times")
    par(mar = c(5, 6, 4, 1) + .1)
    p2 <- plot(predict(model, data_test_fold), data_test_fold$num_clicks, xlab = "Predictions on test data", ylab = "Number of Clicks", cex.axis = 3, cex.lab = 3, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3), pch = 16, cex = 1.5, axes = FALSE)
    axis(1, cex.axis = 2)
    axis(2, cex.axis = 2)
    p2 + theme(
        plot.title = element_text(family = "serif", size = 18),
        axis.title = element_text(family = "serif", size = 50),
        axis.text = element_text(family = "serif", size = 30)
    )
    abline(a = 0, b = 1, col = "blue", lwd = 5, lty = 2)
    dev.off()
}

create_label_expression <- function(var) 
{
    # function for formatting purposes
    if (grepl("\\^2$", var)) {
        base_var <- sub("\\^2$", "", var)
        return(bquote(italic(.(base_var))^2))
    } else {
        return(bquote(italic(.(var))))
    }
}


survival_prob <- function(model, newdata, filename, var, lvl)
# Make and save survival plots - Figure 2
# Save plots to plots/
{
    label_expression <- create_label_expression(var)
    surv_fit <- survfit(model, newdata = newdata)
    levels <- lvl

    surv_plot <- ggsurvplot(
        surv_fit, 
        data = newdata, 
        xlab = "Clicks", 
        ylab = "Probability of Staying on the Grid"
    )

    surv_plot$plot <- surv_plot$plot +
    theme(
        axis.text.x = element_text(size = 38, family = "serif"),    # Customize axes text sizes
        axis.text.y = element_text(size = 38, family = "serif"),    
        axis.title.x = element_text(size = 44, family = "serif"),   
        axis.title.y = element_text(size = 44, family = "serif"),   
        legend.title = element_text(size = 44, family = "serif"),
        legend.text = element_text(size = 44, family = "serif"),
        legend.position = "top",                                    # Position the legend at the top
    ) + labs(
        fill = label_expression,
        color = label_expression
    ) + scale_color_manual(
        values = c("red", "green", "blue", "orange"),
        labels = levels
    ) + scale_fill_manual(
        values = c("red", "green", "blue", "orange"),
        labels = levels
    )

    full_file_path <- paste0("plots/", filename, ".pdf")
    ggsave(full_file_path, plot = surv_plot$plot, width = 10, height = 10, device = "pdf")
}
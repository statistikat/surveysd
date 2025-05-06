#library(surveysd)
library(data.table)

surveysd:::demo.eusilc
surveysd:::recalib

source("R/rescaled.bootstrap.kommis.R")
source("R/draw.bootstrap.R")


setDTthreads(1)
set.seed(1234)
eusilc <- demo.eusilc(n = 3, prettyNames = TRUE)

dat_boot <- draw.bootstrap(eusilc,  REP = 1000 , hid = "hid",
                           weights = "pWeight",
                           strata = "region", period = "year")

dat_boot_rw <- draw.bootstrap(eusilc, method = "Rao-Wu", REP = 1000 , hid = "hid",
                           weights = "pWeight",
                           strata = "region", period = "year")

# calibrate weight for bootstrap replicates
dat_boot_calib <- recalib(dat_boot, conP.var = "gender", conH.var = "region",
                          verbose = TRUE)

dat_boot_calib_rw <- recalib(dat_boot_rw,  conP.var = "gender", conH.var = "region",
                          verbose = TRUE)








# # calibrate on other variables
# dat_boot_calib <- recalib(dat_boot, conP.var = c("gender", "age"),
#                           conH.var = c("region", "hsize"), verbose = TRUE)
# 
# dat_boot_calib_rw <- recalib(dat_boot_rw, conP.var = c("gender", "age"),
#                           conH.var = c("region", "hsize"), verbose = TRUE)
# 
# # supply contingency tables directly
# conP1 <- xtabs(pWeight ~ age + year, data = eusilc)
# conP2 <- xtabs(pWeight ~ gender + year, data = eusilc)
# conH1 <- xtabs(pWeight ~ region + year,
#                data = eusilc[!duplicated(paste(hid,year))])
# conH2 <- xtabs(pWeight ~ hsize + year,
#                data = eusilc[!duplicated(paste(hid,year))])
# 
# conP <- list(conP1,conP2)
# conH <- list(conH1,conH2)
# dat_boot_calib <- recalib(dat_boot, conP.var = NULL,
#                           conH.var = NULL, conP = conP,
#                           conH = conH, verbose = TRUE)
# 
# dat_boot_calib_rw <- recalib(dat_boot_rw, conP.var = NULL,
#                           conH.var = NULL, conP = conP,
#                           conH = conH, verbose = TRUE)
# 
#
# # calibrate on gender x age
# dat_boot_calib <- recalib(dat_boot, conP.var = list(c("gender", "age")),
#                           conH.var = NULL, verbose = TRUE)
# 
# # identical
# conP1 <- xtabs(pWeight ~ age + gender + year, data = eusilc)
# conP <- list(conP1)
# dat_boot_calib <- recalib(dat_boot, conP.var = NULL,
#                           conH.var = NULL, conP = conP,
#                           conH = NULL, verbose = TRUE)
















######################################################
########## check #####################################

# weight 1 - summary
summary(eusilc$pWeight)
summary(dat_boot_calib$w1)
summary(dat_boot_calib_rw$w1)

#
# all weights - summary
library(data.table)

original_stats <- data.table(
  Methode = "Original",
  Min = min(eusilc$pWeight),
  Median = median(eusilc$pWeight),
  Mean = mean(eusilc$pWeight),
  Max = max(eusilc$pWeight)
)

weight_cols <- grep("^w[0-9]+$", names(dat_boot_calib), value = TRUE)
preston_stats <- t(sapply(weight_cols, function(w) {
  summary(dat_boot_calib[[w]])[c("Min.", "Median", "Mean", "Max.")]
}))
preston_final <- data.table(
  Methode = "Preston (alle Gewichte)",
  Min = mean(preston_stats[, "Min."]),
  Median = mean(preston_stats[, "Median"]),
  Mean = mean(preston_stats[, "Mean"]),
  Max = mean(preston_stats[, "Max."])
)

weight_cols_rw <- grep("^w[0-9]+$", names(dat_boot_calib_rw), value = TRUE)
raowu_stats <- t(sapply(weight_cols_rw, function(w) {
  summary(dat_boot_calib_rw[[w]])[c("Min.", "Median", "Mean", "Max.")]
}))
raowu_final <- data.table(
  Methode = "Rao-Wu (alle Gewichte)",
  Min = mean(raowu_stats[, "Min."]),
  Median = mean(raowu_stats[, "Median"]),
  Mean = mean(raowu_stats[, "Mean"]),
  Max = mean(raowu_stats[, "Max."])
)

results <- rbind(original_stats, preston_final, raowu_final)
print(results)

# all weighs - all variables - mean + variance
vars_to_compare <- c("eqIncome", "age",  "gender", "hy130n", "py010n", "py050n", "py090n", "py100n")
w_cols_standard <- grep("^w[0-9]+$", names(dat_boot), value = TRUE)
w_cols_rw <- grep("^w[0-9]+$", names(dat_boot_rw), value = TRUE)

compare_variable_means_and_vars <- function(varname) {
  x <- eusilc[[varname]]
  weights_orig <- eusilc$pWeight
  
  wm <- function(x, w) {
    if (is.numeric(x) || is.integer(x)) {
      return(weighted.mean(x, w, na.rm = TRUE))
    } else if (is.factor(x) || is.logical(x)) {
      x_bin <- as.numeric(x == levels(factor(x))[1])
      return(weighted.mean(x_bin, w, na.rm = TRUE))
    } else if (is.character(x)) {
      x_fact <- factor(x)
      x_bin <- as.numeric(x_fact == levels(x_fact)[1])
      return(weighted.mean(x_bin, w, na.rm = TRUE))
    } else {
      return(NA_real_)
    }
  }
  
  original <- wm(x, weights_orig)
  standard_means <- sapply(w_cols_standard, function(w) wm(x, dat_boot[[w]]))
  rao_wu_means <- sapply(w_cols_rw, function(w) wm(x, dat_boot_rw[[w]]))
  
  dt <- data.table(
    Variable = varname,
    Method = c("Original", "Standard", "Rao-Wu"),
    Mean = c(
      original,
      mean(standard_means, na.rm = TRUE),
      mean(rao_wu_means, na.rm = TRUE)
    ),
    Variance = c(
      NA_real_,  # Original ist Punktwert
      var(standard_means, na.rm = TRUE),
      var(rao_wu_means, na.rm = TRUE)
    )
  )

  dt[, `:=`(
    Mean = formatC(Mean, format = "f", digits = 4),
    Variance = formatC(Variance, format = "f", digits = 4)
  )]
  
  return(dt)
}
compare_all_stats <- rbindlist(lapply(vars_to_compare, compare_variable_means_and_vars))
print(compare_all_stats)


# all weights - all variables - verteilung mittelwerte
library(ggplot2)
get_weighted_means <- function(varname) {
  x <- eusilc[[varname]]
  orig <- weighted.mean(x, eusilc$pWeight, na.rm = TRUE)
  
  means_std <- sapply(w_cols_standard, function(w) weighted.mean(x, dat_boot[[w]], na.rm = TRUE))
  means_rw  <- sapply(w_cols_rw, function(w) weighted.mean(x, dat_boot_rw[[w]], na.rm = TRUE))
  
  data.table(
    Variable = varname,
    Mean = c(means_std, means_rw),
    Methode = rep(c("Standard", "Rao-Wu"), c(length(means_std), length(means_rw))),
    Original = orig
  )
}
vars_to_plot <- c("eqIncome", "hy130n", "py010n", "py050n", "py090n", "py100n")
dt_means_all <- rbindlist(lapply(vars_to_plot, get_weighted_means))

ggplot(dt_means_all, aes(x = Mean, fill = Methode, color = Methode)) +
  geom_density(alpha = 0.3) +
  geom_vline(aes(xintercept = Original), linetype = "dashed", color = "black") +
  facet_wrap(~Variable, scales = "free") +
  labs(
    title = "Verteilung der gewichteten Mittelwerte (Bootstrap)",
    x = "Gewichteter Mittelwert",
    y = "Dichte"
  ) +
  theme_minimal()

















#######################












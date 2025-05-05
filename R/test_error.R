library(surveysd)
library(data.table)
setDTthreads(1)
set.seed(1234)
eusilc <- demo.eusilc(n = 3, prettyNames = TRUE)

dat_boot <- draw.bootstrap(eusilc,  REP = 2, hid = "hid",
                           weights = "pWeight",
                           strata = "region", period = "year")

dat_boot_rw <- draw.bootstrap(eusilc, method = "Rao-Wu", REP = 2, hid = "hid",
                           weights = "pWeight",
                           strata = "region", period = "year")

# calibrate weight for bootstrap replicates
dat_boot_calib <- recalib(dat_boot, conP.var = "gender", conH.var = "region",
                          verbose = TRUE)

dat_boot_calib_rw <- recalib(dat_boot_rw,  conP.var = "gender", conH.var = "region",
                          verbose = TRUE)

# calibrate on other variables
dat_boot_calib <- recalib(dat_boot, conP.var = c("gender", "age"),
                          conH.var = c("region", "hsize"), verbose = TRUE)

dat_boot_calib_rw <- recalib(dat_boot_rw, conP.var = c("gender", "age"),
                          conH.var = c("region", "hsize"), verbose = TRUE)

# supply contingency tables directly
conP1 <- xtabs(pWeight ~ age + year, data = eusilc)
conP2 <- xtabs(pWeight ~ gender + year, data = eusilc)
conH1 <- xtabs(pWeight ~ region + year,
               data = eusilc[!duplicated(paste(hid,year))])
conH2 <- xtabs(pWeight ~ hsize + year,
               data = eusilc[!duplicated(paste(hid,year))])

conP <- list(conP1,conP2)
conH <- list(conH1,conH2)
dat_boot_calib <- recalib(dat_boot, conP.var = NULL,
                          conH.var = NULL, conP = conP,
                          conH = conH, verbose = TRUE)

dat_boot_calib_rw <- recalib(dat_boot_rw, conP.var = NULL,
                          conH.var = NULL, conP = conP,
                          conH = conH, verbose = TRUE)


# calibrate on gender x age
dat_boot_calib <- recalib(dat_boot, conP.var = list(c("gender", "age")),
                          conH.var = NULL, verbose = TRUE)

# identical
conP1 <- xtabs(pWeight ~ age + gender + year, data = eusilc)
conP <- list(conP1)
dat_boot_calib <- recalib(dat_boot, conP.var = NULL,
                          conH.var = NULL, conP = conP,
                          conH = NULL, verbose = TRUE)
















######################################################
########## check #####################################

# income
compare_means <- data.table(
  Method = c("Original", "Standard", "Rao-Wu"),
  Mean = c(
    weighted.mean(eusilc$eqIncome, eusilc$pWeight),
    mean(c(
      weighted.mean(eusilc$eqIncome, dat_boot$w1),
      weighted.mean(eusilc$eqIncome, dat_boot$w2)
    )),
    mean(c(
      weighted.mean(eusilc$eqIncome, dat_boot_rw$w1),
      weighted.mean(eusilc$eqIncome, dat_boot_rw$w2)
    ))
  )
)
print(compare_means)

# gender
eusilc[, female := gender == "female"]
compare_gender <- data.table(
  Method = c("Original", "Standard", "Rao-Wu"),
  FemaleShare = c(
    weighted.mean(eusilc$female, eusilc$pWeight),
    mean(c(
      weighted.mean(eusilc$female, dat_boot$w1),
      weighted.mean(eusilc$female, dat_boot$w2)
    )),
    mean(c(
      weighted.mean(eusilc$female, dat_boot_rw$w1),
      weighted.mean(eusilc$female, dat_boot_rw$w2)
    ))
  )
)

# region
prop.table(xtabs(pWeight ~ region, data = eusilc))
prop.table(xtabs(pWeight * w1 ~ region, data = dat_boot))



#######################

boot_means_std <- c(
  weighted.mean(eusilc$eqIncome, dat_boot$w1),
  weighted.mean(eusilc$eqIncome, dat_boot$w2)
)

# Rao-Wu-Bootstrap
boot_means_rw <- c(
  weighted.mean(eusilc$eqIncome, dat_boot_rw$w1),
  weighted.mean(eusilc$eqIncome, dat_boot_rw$w2)
)

var_comparison <- data.table(
  Method = c("Standard", "Rao-Wu"),
  Variance = c(var(boot_means_std), var(boot_means_rw))
)
print(var_comparison)











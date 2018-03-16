############################################################
#
# BATCH-Script for sensitivity analysis
#
#

# draw from udb_short_calib (contains 1000 bootstrap replicates)
# then calc estimates
# divide result into groups number of observations in each group
# <5 - 5-15 - 15-30 - 31-50 - 51-80 - 80-120 - 121+
# look at standard error as well as tails of distribution (.025,.975)
# plot results

library(data.table)
library(surveysd)
library(ggplot2)
library(mountSTAT)
pfad_meth <- mountWinShare("DatenREG","REG_METHODIK","meth")[1]

dat_boot_calib <- fread(paste0(mountO(),"/B/Datenaustausch/NETSILC3/udb_ES_calib.csv"))

B <- 1000
res <- list()
for(i in 1:B){

  res_i <- calc.stError(dat=dat_boot_calib,weights="RB050",year="RB010",b.weights=paste0("w",1:i),
                        var="arose",cross_var=list(c("DB040","DB100")),year.diff=c("2014-2008"),
                        p=c(.01,0.025,.05,.1,.9,.95,0.975,.99))
  res_i <- res_i$Estimates
  res_i[,NumberWeights:=i]

  res <- c(res,list(res_i))
}

save(res,file=paste0(pfad_meth,"/Gussenbauer/surveysd/Results_i1000_arose.RData"))


B <- 100
nB <- seq(50,800,by=50)

for(b in nB){
  res_boot <- list()
  for(i in 1:B){

    random_w <- sample(1:1000,b)

    res_i <- calc.stError(dat=dat_boot_calib,weights="RB050",year="RB010",b.weights=paste0("w",random_w),
                          var="HX080",cross_var=list(c("DB040","DB100")),year.diff=c("2014-2008"),
                          p=c(.01,0.025,.05,.1,.9,.95,0.975,.99))
    res_i <- res_i$Estimates
    res_i[,RUN:=i]

    res_boot <- c(res_boot,list(res_i))

  }
  save(res_boot,file=paste0(pfad_meth,"/Gussenbauer/surveysd/Results_boot",b,"_arose.RData"))
}


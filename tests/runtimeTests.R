 library(surveysd)
 library(laeken)
 library(data.table)
 library(microbenchmark)

 eusilc <- surveysd:::demo.eusilc()

 dat_boot <- draw.bootstrap(eusilc,REP=1000,hid="db030",weights="rb050",strata=c("db040"),
                            period="year")

 # calibrate weight for bootstrap replicates
 dat_boot_calib <- recalib(dat=copy(dat_boot),hid="db030",weights="rb050",
                           period="year",b.rep=paste0("w",1:1000),
                           conP.var=c("rb090"),conH.var = c("db040"))

 # estimate weightedRatio for povmd60 per period
 
 b.weights <- paste0("w",1:50)
 var <- "povmd60"
 weights <- "rb050"
 groups <- c("year","db040")
 by.eval <- paste(c(groups),collapse=",")
 
 res.names <- c(t(outer(var, 1:length(c(weights,b.weights)), surveysd:::paste_)))
 eval.fun <- paste0(res.names,"=fun(",c(t(outer(var,c(weights,b.weights), surveysd:::paste_c))),")")
 eval.fun <- paste0(".(",paste(eval.fun,collapse=","),")")

 fun <- weightedRatio
 a <- surveysd:::dt.eval("dat_boot_calib[,",eval.fun,",by=list(",by.eval,")]")
 b <- dat_boot_calib[, lapply(c(weights,b.weights), function(w) { fun(.SD[[var]],.SD[[w]]) }), by = groups]
 
 WeightedRatio_return <- function(condition, na.rm = FALSE) { 
   function(subset, w) { 
     subset[, sum(w * eval(parse(text = condition)), na.rm = na.rm)/sum(w)*100] 
   } 
 } 
 FUN_return <- WeightedRatio_return("povmd60 == 1")
 c <- dat_boot_calib[, lapply(c(weights,b.weights), function(w) { FUN_return(.SD,.SD[[w]]) }), by = groups]
 
 all.equal(a,b,check.attributes=FALSE)
 all.equal(a,c,check.attributes=FALSE)

 
 # Ergebins für 50 Bootstrap replicates + orig weight
 # Übergebe ganzen Datensatz:
 dat_test <- dat_boot_calib[,mget(c(groups,var,weights,paste0("w",1:1000)))]

 microbenchmark(parse_text={surveysd:::dt.eval("dat_boot_calib[,",eval.fun,",by=list(",by.eval,")]")},
                lapply={ dat_boot_calib[, lapply(c(weights,b.weights), function(w) { fun(.SD[[var]],.SD[[w]]) }), by = groups]})
 
 
 # Übergebe teil von Datensatz:
 dat_test <- dat_boot_calib[,mget(c(groups,var,weights,b.weights))]
 
 microbenchmark(parse_text={surveysd:::dt.eval("dat_boot_calib[,",eval.fun,",by=list(",by.eval,")]")},
                lapply={ dat_boot_calib[, lapply(c(weights,b.weights), function(w) { fun(.SD[[var]],.SD[[w]]) }), by = groups]},
                FUN_return={dat_boot_calib[, lapply(c(weights,b.weights), function(w) { FUN_return(.SD,.SD[[w]]) }), by = groups]})
 
 
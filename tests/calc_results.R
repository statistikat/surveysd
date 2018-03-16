#############################################################################
# Ergebnistabellen Erstellen

library(data.table)
library(surveysd)

if(Sys.info()['sysname']=="Windows"){
  pfadO <- "O:"
}else{
  library(mountSTAT)
  pfadO <- mountO()[1]
}

############################################################
# Ergebnisse für Spanien
#

dat_boot_calib <- fread(file=paste0(pfadO,"/B/Datenaustausch/NETSILC3/udb_ES_calib.csv"))

es_results <- calc.stError(dat_boot_calib,weights="RB050",b.weights=paste0("w",1:250),
                           year="RB010",var="arose",cross_var=list(c("DB040","DB100")),year.diff = c("2014-2008","2015-2009"))


results_part1 <- es_results$Estimates[!is.na(DB040)&!is.na(DB100)][RB010%in%c("2009","2015","2015-2009",
                                                                              "2008_2009_2010","2014_2015_2016",
                                                                              "2015-2009-mean")]

# Differenz zu nationalem 3JDS
results_nat <- es_results$Estimates[is.na(DB040)&is.na(DB100)&RB010%in%c("2008_2009_2010","2014_2015_2016")]
results_part2 <- results_part1[RB010%in%c("2008_2009_2010","2014_2015_2016")]
results_part2[,val_arose:=val_arose-results_nat[RB010==.BY,val_arose],by=RB010]
results_part2[,RB010:=paste0(RB010,"_Nat")]

# Differenz zu EU 3JDS
results_part3 <- results_part1[RB010%in%c("2008_2009_2010","2014_2015_2016")]
results_part3[RB010%in%c("2008_2009_2010"),val_arose:=val_arose-21.8]
results_part3[RB010%in%c("2014_2015_2016"),val_arose:=val_arose-23.2]
results_part3[,RB010:=paste0(RB010,"_EU")]

# Werte verwerfen von (geschätztes) KI größer als 33% vom Punktschätzer ist
results_ES <- rbindlist(list(results_part1,results_part2,results_part3),use.names=TRUE)
# results_ES[abs(stE_arose*1.96/val_arose)>=.33,val_arose:=NA]
results_ES[abs(stE_arose*1.65/val_arose)>1,val_arose:=NA]
results_ES[RB010%in%c("2009","2015","2008_2009_2010","2014_2015_2016")&abs(stE_arose*1.96/val_arose)>=.33,val_arose:=NA]

results_part1[RB010=="2015-2009-mean"][sign(p0.025_arose)==sign(p0.975_arose)]

############################################################
# Ergebnisse für AT
#

dat_boot_calib <- fread(file=paste0(pfadO,"/B/Datenaustausch/NETSILC3/udb_AT_calib.csv"))

at_results <- calc.stError(dat_boot_calib,weights="RB050",b.weights=paste0("w",1:250),
                           year="RB010",var="arose",cross_var=list(c("DB040","DB100")),year.diff = c("2014-2008","2015-2009"))


results_part1 <- at_results$Estimates[!is.na(DB040)&!is.na(DB100)][RB010%in%c("2009","2015","2015-2009",
                                                                              "2008_2009_2010","2014_2015_2016",
                                                                              "2015-2009-mean")]

# Differenz zu nationalem 3JDS
results_nat <- at_results$Estimates[is.na(DB040)&is.na(DB100)&RB010%in%c("2008_2009_2010","2014_2015_2016")]
results_part2 <- results_part1[RB010%in%c("2008_2009_2010","2014_2015_2016")]
results_part2[,val_arose:=val_arose-results_nat[RB010==.BY,val_arose],by=RB010]
results_part2[,RB010:=paste0(RB010,"_Nat")]

# Differenz zu EU 3JDS
results_part3 <- results_part1[RB010%in%c("2008_2009_2010","2014_2015_2016")]
results_part3[RB010%in%c("2008_2009_2010"),val_arose:=val_arose-21.8]
results_part3[RB010%in%c("2014_2015_2016"),val_arose:=val_arose-23.2]
results_part3[,RB010:=paste0(RB010,"_EU")]

# Werte verwerfen von (geschätztes) KI größer als 33% vom Punktschätzer ist
results_AT <- rbindlist(list(results_part1,results_part2,results_part3),use.names=TRUE)
# results_AT[abs(stE_arose*1.96/val_arose)>=.33,val_arose:=NA]
results_AT[abs(stE_arose*1.65/val_arose)>1,val_arose:=NA]
results_AT[RB010%in%c("2009","2015","2008_2009_2010","2014_2015_2016")&abs(stE_arose*1.96/val_arose)>=.33,val_arose:=NA]

###############################################################
# Ergebnisse kombinieren und in xlsx herausschreiben
results_all <- rbind(results_AT,results_ES)
results_all[,val_arose:=round(val_arose,digits=1)]
results_all <- dcast(results_all,DB040+DB100~RB010,value.var="val_arose")

oldnames <- colnames(results_all)[-(1:2)]
newnames <- c("2009(AAA)","DEA2009(AAA)","DNA2009(AAA)","2009","2015(AAA)","DEA2015(AAA)","DNA2015(AAA)","2015","diff","diff(AAA)")

setnames(results_all,oldnames,newnames)

setcolorder(results_all,c("DB040","DB100","2015","2009","diff","2015(AAA)","2009(AAA)","diff(AAA)","DNA2015(AAA)","DNA2009(AAA)","DEA2015(AAA)","DEA2009(AAA)"))

fwrite(results_all,file=paste0(pfadO,"/B/Datenaustausch/NETSILC3/udb_ES_AT_Ergebnisse.csv"),sep=";",dec=",")


fwrite(es_results$Estimates,file=paste0(pfadO,"/B/Datenaustausch/NETSILC3/udb_ES_Resultate.csv"),sep=";",dec=",")
fwrite(at_results$Estimates,file=paste0(pfadO,"/B/Datenaustausch/NETSILC3/udb_AT_Resultate.csv"),sep=";",dec=",")

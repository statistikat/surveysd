###########################################################
# SENSITIVITY ANALYSIS FOR NUMBER OF BOOTSTRAP REPLICATES
#
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

load("udb_short_calib.RData")

B <- 1000
res <- list()
for(i in 282:B){

  res_i <- calc.stError(dat=dat_boot_calib,weights="rb050",year="rb010",b.weights=paste0("w",1:i),
               var="hx080",cross_var=list(c("db040","db100")),year.diff=c("2014-2008"),
               p=c(.01,0.025,.05,.1,.9,.95,0.975,.99))
  res_i <- res_i$Estimates
  res_i[,NumberWeights:=i]

  res <- c(res,list(res_i))
}

save(res,file="Results_i281.RData")

res <- rbindlist(res)
res <- res[rb010%in%dat_boot_calib[,unique(rb010)]]
res[,rb010:=as.numeric(rb010)]

num_obs <- rbindlist(list(dat_boot_calib[,.N,by=rb010],dat_boot_calib[,.N,by=list(rb010,db040,db100)]),use.names=TRUE,fill=TRUE)
num_obs[,GROUP:=cut(N,c(-Inf,10,25,50,100,200,500,1000,Inf))]

res <- merge(res,num_obs,by=c("rb010","db040","db100"))
res[,GROUP_VAR:=.GRP,by=c("rb010","db040","db100")]

estim <- c("stE_hx080","p0.01_hx080","p0.025_hx080","p0.05_hx080","p0.1_hx080","p0.9_hx080","p0.95_hx080","p0.975_hx080","p0.99_hx080")
sel_var <- c("rb010","db040","db100","NumberWeights","GROUP_VAR","GROUP")

for( i in 1:length(estim)){
  dat_plot <- subset(res,select=c(sel_var,estim[i]))
  # berechne Veränderungsrate vom estimate über NumberWeights
  setkeyv(dat_plot,c("GROUP_VAR","NumberWeights"))
  dat_plot <- dat_plot[,.(EST_CHANGE=diff(get(estim[i]))/get(estim[i])[-1],NumberWeights=NumberWeights[-1]),by=list(GROUP_VAR,GROUP)]

  p1 <- ggplot(dat_plot,aes(NumberWeights,EST_CHANGE,group=GROUP_VAR))+
    geom_line()+coord_cartesian(xlim = c(-1, 1))+ggtitle(estim[i])+
    facet_grid(GROUP~.)
  plot(p1)

  readline(prompt="Press [enter] to continue")
}




for(i in 1:B){

  random_w <- sample(1:1000,50)

  calc.stError(dat=dat_boot_calib,weights="rb050",year="rb010",b.weights=paste0("w",random_w),
               var="hx080",cross_var=list(c("db040","db100")),year.diff=c("2014-2008"),
               p=c(.01,0.025,.05,.1,.9,.95,0.975,.99))



}






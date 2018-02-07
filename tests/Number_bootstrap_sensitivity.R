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
for(i in 456:B){

  res_i <- calc.stError(dat=dat_boot_calib,weights="rb050",year="rb010",b.weights=paste0("w",1:i),
               var="hx080",cross_var=list(c("db040","db100")),year.diff=c("2014-2008"),
               p=c(.01,0.025,.05,.1,.9,.95,0.975,.99))
  res_i <- res_i$Estimates
  res_i[,NumberWeights:=i]

  res <- c(res,list(res_i))
}

save(res,file="Results_i1000.RData")
save(res,file="Results_i455.RData")

load("Results_i455.RData")

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
  dat_plot <- na.omit(dat_plot[,.(EST_CHANGE=diff(get(estim[i])),NumberWeights=NumberWeights[-1]),by=list(GROUP_VAR,GROUP)])
  dat_plot <- dat_plot[,.(Q=c("min","med","max"),quantile(EST_CHANGE,c(0,.5,1))),by=list(GROUP,NumberWeights)]
  dat_plot <- dcast(dat_plot,GROUP+NumberWeights~Q,value.var="V2")
  dat_plot[,min:=med-min]
  dat_plot[,max:=max-med]
  setnames(dat_plot,c("GROUP"),c("NumberObs"))

  p1 <- ggplot(dat_plot,aes(NumberWeights,med,group=NumberObs))+
    geom_line(aes(color=NumberObs),size=.25)+xlab("Number Bootstrap weights")+ylab("Min - Median - Max")+
    geom_ribbon(aes(ymin = med-min, ymax = med+max,fill=NumberObs,colour=NumberObs),linetype = 2, alpha= 0.1)+
    coord_cartesian(ylim=c(1,-1))+
    facet_grid(NumberObs~.)+ggtitle(estim[i])

  plot(p1)

  readline(prompt="Press [enter] to continue")
}


# draw boostrap weights randomly and look at results

B <- 100
nB <- seq(50,800,by=50)

for(b in nB){
  res_boot <- list()
  for(i in 1:B){

    random_w <- sample(1:1000,b)

    res_i <- calc.stError(dat=dat_boot_calib,weights="rb050",year="rb010",b.weights=paste0("w",random_w),
                          var="hx080",cross_var=list(c("db040","db100")),year.diff=c("2014-2008"),
                          p=c(.01,0.025,.05,.1,.9,.95,0.975,.99))
    res_i <- res_i$Estimates
    res_i[,RUN:=i]

    res_boot <- c(res_boot,list(res_i))

  }
  save(res_boot,file=paste0("Results_boot",b,".RData"))
}


res_all <- list()
for(b in nB){
  load(paste0("Results_boot",b,".RData"))
  res_boot <- rbindlist(res_boot)
  res_boot[,B:=b]
  res_all <- c(res_all,list(res_boot))
}
res_all <- rbindlist(res_all)


res_all  <- res_all[rb010%in%dat_boot_calib[,unique(rb010)]]
res_all[,rb010:=as.numeric(rb010)]

num_obs <- rbindlist(list(dat_boot_calib[,.N,by=rb010],dat_boot_calib[,.N,by=list(rb010,db040,db100)]),use.names=TRUE,fill=TRUE)
num_obs[,GROUP:=cut(N,c(-Inf,10,25,50,100,200,500,1000,Inf))]

res_all  <- merge(res_all ,num_obs,by=c("rb010","db040","db100"))
res_all[,GROUP_VAR:=.GRP,by=c("rb010","db040","db100")]

estim <- c("stE_hx080","p0.01_hx080","p0.025_hx080","p0.05_hx080","p0.1_hx080","p0.9_hx080","p0.95_hx080","p0.975_hx080","p0.99_hx080")

for( i in 1:length(estim)){
  plot_boot <- res_all[,mget(c(estim[i],"GROUP","GROUP_VAR","B"))]
  plot_boot <- plot_boot[,diff(range(get(estim[i]))),by=list(GROUP_VAR,GROUP,B)]
  plot_boot <- plot_boot[,.(Q=c("p.01","p.5","p.99"),quantile(V1,c(.01,.5,.99))),by=list(GROUP,B)]
  plot_boot <- dcast(plot_boot,GROUP+B~Q,value.var="V2")
  plot_boot[,p.01:=p.5-p.01]
  plot_boot[,p.99:=p.99-p.5]
  setnames(plot_boot,c("GROUP"),c("NumberObs"))

  p1 <- ggplot(plot_boot,aes(factor(B),p.5,group=NumberObs))+
    geom_line(aes(color=NumberObs),size=.25)+xlab("Number Bootstrap weights")+ylab("1% - Median - 99%")+
    geom_ribbon(aes(ymin = p.5-p.01, ymax = p.5+p.99,fill=NumberObs,colour=NumberObs),linetype = 2, alpha= 0.1)+
    facet_grid(NumberObs~.)+ggtitle(estim[i])

  plot(p1)

  readline(prompt="Press [enter] to continue")
}










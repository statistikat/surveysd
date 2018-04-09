#' Plot surveysd-Objects
#'
#' Plot results of \code{calc.stError()}
#'
#' @param dat object of class 'surveysd' output of function \code{\link{calc.stError}}
#' @param variable Name of the variable for which standard errors have been calcualated in \code{dat}
#' @param type can bei either 'summary' or 'grouping', default value is 'summary'.
#' For 'summary' a barplot is created giving an overview of the number of estimates having the flag \code{smallGroup}, \code{cvHigh}, both or none of them.
#' For 'grouping' results for point estimate and standard error are plotted for pre defined groups.
#' @param groups If \code{type='grouping'} variables must be defined by which the data is grouped. Only 2 levels are supported as of right now.
#' If only one group is defined the higher group will be the estimate over the whole year.
#' Results are plotted for the first argument in \code{groups} as well as for the combination of \code{groups[1]} and \code{groups[2]}.
#' @param sd.type can bei either \code{'ribbon'} or \code{'dot'} and is only used if \code{type='grouping'}. Default is \code{"dot"}
#' For \code{sd.type='dot'} point estimates are plotted and flagged if the corresponding standard error and/or the standard error using the mean over k-years exceeded the value \code{cv.limit} (see \code{\link{calc.stError}}).
#' For \code{sd.type='ribbon'} the point estimates including ribbons, defined by point estimate +- estimated standard error are plotted.
#' The calculated standard errors using the mean over k years are plotted using less transparency. Results for the higher level (~\code{groups[1]}) are coloured grey.
#'
#' @export
#'

plot.surveysd <- function(dat,variable=dat$param$var[1],type=c("summary","grouping"),
                          groups=NULL,sd.type=c("dot","ribbon")){

  #################
  # Input checking
  type <- type[1]
  sd.type <- sd.type[1]
  if(class(dat)[1]!="surveysd"){
    stop("The data needs to be an object of class 'surveysd'!")
  }
  if(!type%in%c("summary","grouping")){
    stop("Parameter type can only take values 'summary' or 'grouping'!")
  }
  if(!variable%in%dat$param$var){
    stop("No results for ",variable," present in the data!")
  }
  if(type=="grouping"&is.null(groups)){
    other_var <- unique(unlist(dat$param$cross_var))
    if(is.null(groups)){
      stop("Paramter 'groups' cannot be NULL if type='grouping'!")
    }
    if(any(!groups%in%other_var)){
      stop("Variables in 'groups' must contain variables from dat$params$cross_var!")
    }
  }
  if(!sd.type%in%c("ribbon","dot")){
    stop("Parameter 'sd.type' can only take values 'ribbon' or 'dot'!")
  }

  #################
  # get variables from dat and merge tables together
  year <- dat$param$year

  poss.values <- c("Missing","cvHigh","SmallGroup+cvHigh","SmallGroup","OK")
  poss.color <- c("grey","orangered3","yellow2","dodgerblue1","green4")

  plot.dat <- dat$Estimates
  plot.dat <- define_type(plot.dat,dat,variable=variable)
  values.plot <- plot.dat[,unique(res_type)]
  color.plot <- poss.color[poss.values%in%values.plot]

  # define groups
  dt.eval("plot.dat[grepl('^[[:digit:]]+$',",year,"),GROUPING:='Years']")
  dt.eval("plot.dat[grepl('-',",year,"),GROUPING:='Other']")
  dt.eval("plot.dat[grepl('_',",year,"),GROUPING:='",dat$param$year.mean,"-Year-Mean']")
  plot.dat[,GROUPING:=factor(GROUPING,levels=c("Years",paste0(dat$param$year.mean,"-Year-Mean"),"Other"))]

  #################
  # plot type 'summary'
  if(type=="summary"){
    title <- paste("Results of Groups per year for variable",variable)
    p1 <- ggplot(plot.dat,aes(get(year),fill=res_type))+
      geom_bar()+xlab("")+ylab("Count")+coord_flip()+
      facet_grid(GROUPING~.,scales="free_y",space="free_y")+
      theme(legend.title = element_blank())+
      scale_fill_manual(breaks = values.plot,
                        values=color.plot)+
      ggtitle(title)
    plot(p1)
  }else if(type=="grouping"){
    #################
    # plot type 'grouping'
    if(length(groups)==1){
      groups <- c(year,groups)
    }else{
      groups <- groups[1:2]
    }

    exp1_2 <- paste(paste0("is.na(",other_var[!other_var%in%groups],")"),collapse="&")

    exp1 <- paste0("c(!is.na(",groups[1],")&!is.na(",groups[2],")&",exp1_2,")")

    exp2_2 <- paste(paste0("is.na(",other_var[!other_var%in%groups[1]],")"),collapse="&")
    exp2 <- paste0("c(!is.na(",groups[1],")&",exp2_2,")")
    plot.group1 <- dt.eval("plot.dat[",exp1,"]")
    plot.group2 <- dt.eval("plot.dat[",exp2,"]")

    if(nrow(plot.group1)==0|nrow(plot.group2)==0){
      stop("Cannot create plot since results have not been calculated by ",groups[1]," and/or ",groups[2],"!")
    }

    val_var <- paste0("val_",variable)
    ste_var <- paste0("stE_",variable)
    ste_var_mean <- paste0(ste_var,"_mean")
    ste_bool <- variable
    ste_bool_mean <- paste0(ste_bool,"_mean")

    # prepare plot.group1 for plotting
    plot.group1.year <- plot.group1[GROUPING=="Years"]
    plot.group1.yearmean <- plot.group1[GROUPING==paste0(dat$param$year.mean,"-Year-Mean")]
    dt.eval("plot.group1.yearmean[,",year,":=tstrsplit(",year,",split='_',keep=",round((dat$param$year.mean+1)/2),")]")

    setnames(plot.group1.yearmean,c(ste_var,ste_bool),c(ste_var_mean,ste_bool_mean))
    plot.group1.year <- merge(plot.group1.year,plot.group1.yearmean[,mget(c(year,other_var,ste_var_mean,ste_bool_mean))],all.x=TRUE)

    # prepare plot.group2 for plotting
    plot.group2.year <- plot.group2[GROUPING=="Years"]
    plot.group2.yearmean <- plot.group2[GROUPING==paste0(dat$param$year.mean,"-Year-Mean")]
    dt.eval("plot.group2.yearmean[,",year,":=tstrsplit(",year,",split='_',keep=",round((dat$param$year.mean+1)/2),")]")

    setnames(plot.group2.yearmean,c(ste_var,ste_bool),c(ste_var_mean,ste_bool_mean))
    plot.group2.year <- merge(plot.group2.year,plot.group2.yearmean[,mget(c(year,other_var,ste_var_mean,ste_bool_mean))],all.x=TRUE)

    ste_var2 <- paste0(ste_var,"2")
    ste_var_mean2 <- paste0(ste_var_mean,"2")
    val_var2 <- paste0(val_var,"2")
    ste_bool2 <- paste0(variable,"2")
    ste_bool_mean2 <- paste0(ste_bool_mean,"2")
    setnames(plot.group2.year,c(ste_var,ste_var_mean,val_var,ste_bool,ste_bool_mean),c(ste_var2,ste_var_mean2,val_var2,ste_bool2,ste_bool_mean2))
    # merge with plot.group

    if(year==groups[1]){
      plot.group1.year <- merge(plot.group1.year,plot.group2.year[,mget(c(groups[1],ste_var2,ste_var_mean2,val_var2,ste_bool2,ste_bool_mean2))],
                                all.x=TRUE,by=groups[1])
    }else{
      plot.group1.year <- merge(plot.group1.year,plot.group2.year[,mget(c(groups[1],year,ste_var2,ste_var_mean2,val_var2,ste_bool2,ste_bool_mean2))],
                                all.x=TRUE,by=c(groups[1],year))
    }

    dt.eval("plot.group1.year[,",groups[1],":=factor(",groups[1],")]")
    dt.eval("plot.group1.year[,",groups[2],":=factor(",groups[2],")]")

    # plot results of subgroups and add results for mean over k years
    p1 <- ggplot(plot.group1.year,aes(get(year),get(val_var),group=get(groups[2])))+
      geom_line(aes(colour=get(groups[2])))

    if(sd.type=="ribbon"){
      p1 <- p1 + geom_ribbon(aes(ymin = get(val_var)-get(ste_var),
                                 ymax = get(val_var)+get(ste_var),
                                 fill=get(groups[2]),colour=get(groups[2])),linetype = 2, alpha= 0.1)

      # add results for k-mean-years
      p1 <- p1 + geom_ribbon(aes(ymin = get(val_var)-get(ste_var_mean),
                                 ymax = get(val_var)+get(ste_var_mean),
                                 fill=get(groups[2]),colour=get(groups[2])),linetype = 2, alpha= 0.5)

      # add results for group[1] -> higher level grouping
      p1 <- p1 + geom_line(aes(get(year),get(val_var2)),color="black")+
        geom_ribbon(aes(ymin = get(val_var2)-get(ste_var2),
                        ymax = get(val_var2)+get(ste_var2)),fill="grey",linetype = 2, alpha= 0.1)

      # add results for group[1] and k-year-mean
      p1 <- p1 + geom_ribbon(aes(ymin = get(val_var2)-get(ste_var_mean2),
                                 ymax = get(val_var2)+get(ste_var_mean2)),fill="grey",linetype = 2, alpha= 0.5)
    }else{
      plot.group1.year[,shape_bool:="Low st.Error"]
      plot.group1.year[get(ste_bool)==TRUE,shape_bool:="high st.Error"]
      plot.group1.year[get(ste_bool_mean)==TRUE,shape_bool:="high st.Error for mean"]

      plot.group1.year[,shape_bool2:="Low st.Error"]
      plot.group1.year[get(ste_bool2)==TRUE,shape_bool2:="high st.Error"]
      plot.group1.year[get(ste_bool_mean2)==TRUE,shape_bool2:="high st.Error for mean"]

      p1 <- p1 + geom_point(data=plot.group1.year[get(ste_bool)==TRUE],
                            aes(get(year),get(val_var),shape=shape_bool))

      # add results for k-mean-years
      p1 <- p1 + geom_point(data=plot.group1.year[get(ste_bool_mean)==TRUE],
                            aes(get(year),get(val_var),shape=shape_bool))

      # add results for group[1] -> higher level grouping
      p1 <- p1 + geom_line(aes(get(year),get(val_var2),linetype="main"))
      p1 <- p1 + geom_point(data=plot.group1.year[get(ste_bool2)==TRUE],
                            aes(get(year),get(val_var),shape=shape_bool2))

      # add results for group[1] and k-year-mean
      p1 <- p1 + geom_point(data=plot.group1.year[get(ste_bool_mean2)==TRUE],
                            aes(get(year),get(val_var),shape=shape_bool2))

    }


    # define paramters for plot
    ylabel <- paste0(dat$param$fun," of ",variable)
    p1 <- p1 + xlab("")+ylab(ylabel)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))


    # split the plots into facets
    if(groups[1]==year){
      p1 <- p1 + facet_wrap(~get(groups[2]))
      if(sd.type=="ribbon"){
        p1 <- p1 +  theme(legend.position = "none")
      }else{
        shapes <- c(1,16)
        names(shapes) <- c("high st.Error","high st.Error for mean")
        cols <- c("main"="grey")
        p1 <- p1 + scale_shape_manual(values=shapes)
        p1 <- p1 + scale_colour_discrete(guide = FALSE)+
        # p1 <- p1 + scale_colour_discrete(value="grey")
           theme(legend.title=element_blank())
        p1 <- p1 + scale_linetype_manual(values=c("dotted"),labels=paste0("Result for ",groups[1]))
      }
    }else{
      p1 <- p1 + facet_wrap(~get(groups[1]))
      if(sd.type=="ribbon"){
        p1 <- p1 + guides(fill=guide_legend(title=groups[2]))+
          guides(colour=guide_legend(title=groups[2]))
      }else{
        shapes <- c(1,16)
        names(shapes) <- c("high st.Error","high st.Error for mean")
        cols <- c("main"="grey")
        p1 <- p1 + scale_shape_manual(values=shapes)+
          guides(colour=guide_legend(title=groups[2]),
                 linetype=guide_legend(title=""),
                 shape=guide_legend(title=""))

        p1 <- p1 + scale_linetype_manual(values=c("dotted"),labels=paste0("Result for ",groups[1]))
      }
    }

    plot(p1)

  }
}


define_type <- function(plot.dat,dat,variable="HX080"){

  val_variable <- paste("val",variable,sep="_")
  poss.values <- c("Missing","cvHigh","SmallGroup+cvHigh","SmallGroup","OK")

  # merge data
  plot.dat <- merge(plot.dat,dat$smallGroups,all.x=TRUE)
  setkeyv(dat$cvHigh,unique(c(dat$param$year,unlist(dat$param$cross_var))))
  plot.dat <- dat$cvHigh[plot.dat]

  # define type
  plot.dat[!is.na(N)&get(variable)==FALSE,res_type:="SmallGroup"]
  plot.dat[!is.na(N)&get(variable)==TRUE,res_type:="SmallGroup+cvHigh"]
  plot.dat[is.na(get(val_variable)),res_type:="Missing"]
  plot.dat[get(variable)==TRUE&is.na(res_type),res_type:="cvHigh"]
  plot.dat[is.na(res_type),res_type:="OK"]
  plot.dat[,res_type:=factor(res_type,levels=poss.values)]

  return(plot.dat)
}


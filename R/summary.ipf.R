#' Generate summary output for a ipf calibration
#'
#' @param x object of class ipf
#' @param ... additional arguments
#'
#' @author Laura Gruber
#' @return a list of the following outputs
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # load data
#' eusilc <- demo.eusilc(n = 1, prettyNames = TRUE)
#'
#' # personal constraints
#' conP1 <- xtabs(pWeight ~ age, data = eusilc)
#' conP2 <- xtabs(pWeight ~ gender + region, data = eusilc)
#' conP3 <- xtabs(pWeight*eqIncome ~ gender, data = eusilc)
#'
#' # household constraints
#' conH1 <- xtabs(pWeight ~ hsize + region, data = eusilc)
#'
#' # simple usage ------------------------------------------
#'
#' calibweights1 <- ipf(
#'   eusilc,
#'   conP = list(conP1, conP2, eqIncome = conP3),
#'   bound = NULL,
#'   verbose = TRUE
#' )
#' output <- summary(calibweights1)
#' # the output can easily be exported to an Excel file, e.g. with
#' # library(openxlsx)
#' # write.xlsx(output, "SummaryIPF.xlsx")
#' }
summary.ipf <- function(object, ...){
  terms <- variable <- NULL
  dots <- list(...)
  extraExport <- dots$extraExport
  av <- attributes(x)
  w <- av$baseweight
  hid <- av$hid
  output <- vector(mode = "list", length = 3)
  av <- attributes(x)

  ## ---- FormP/FormH - Variablen etc. ---- ##
  vars <- c()
  formulas <- c()
  if(any(names(av) == "conP")){
    calibWeightName <- as.character(av$formP[1][[1]])[2]
    varsP_help <- c()
    for(i in seq_along(av$formP)){
      varsP_help <- c(varsP_help,labels(terms(av$formP[[i]])))
    }
    vars <- c(vars,unique(varsP_help))

    #
    formPBase <- lapply(av$formP,function(x){
      x <- as.character(x)
      as.formula(paste(w,"~",x[3]))
    })
    #
    formulas <- c(formulas,as.character(av$formP))

  }

  if(any(names(av) == "conH")){
    calibWeightName <- as.character(av$formH[1][[1]])[2]
    varsH_help <- c()
    for(i in seq_along(av$formH)){
      varsH_help <- c(varsH_help,labels(terms(av$formH[[i]])))
    }
    vars <- c(vars,unique(varsH_help))

    #
    formHBase <- lapply(av$formH,function(x){
      x <- as.character(x)
      as.formula(paste(w,"~",x[3]))
    })
    #
    formulas <- c(formulas,as.character(av$formH))

  }
  vars <- as.list(unique(vars))

  dist_gew <- list()
  if(!is.null(w)){
    for(i in vars){
      part_i <- x[,list(cv=sd(get(w))/mean(get(w)),
                        min=min(get(w)),
                        max=max(get(w)),
                        quotient=min(get(w))/max(get(w))),keyby=c(i)]
      setnames(part_i,i,paste0("value",1:length(i)))
      part_i[,variable:=paste(i,collapse=" x ")]
      setcolorder(part_i,c("variable",paste0("value",1:length(i)),"cv","min","max","quotient"))
      dist_gew <- c(dist_gew,list(part_i))
    }  
    dist_gew <- rbindlist(dist_gew,use.names=TRUE,fill=TRUE)
    setcolorder(dist_gew,c("variable",paste0("value",1:length(i)),"cv","min","max","quotient"))
  }
  
  

  cols <- c(which(colnames(x) %in% c(hid,vars,w,calibWeightName,extraExport)))
  names(output)[1] <- "weighted data"
  output[[1]] <- subset(x,select=cols)


  #beinhaltet ebenfalls formP und formH
  formulaOutput <- data.table(Verteilungen=rep(NA,length(formulas)))
  formulaOutput$Verteilungen[seq_along(formulas)] <- formulas

  ##Speichern der Tabellen
  names(output)[2] <- "margin tables"
  output[[2]] <- formulaOutput

  names(output)[3] <- "distribution of the weights"
  output[[3]] <- dist_gew

  #### ---- ####
  if (any(names(av) == "conP")) {
    for (i in seq_along(av$conP)) {
      output4 <- list(
        "conP_%i" = as.data.table(av$conP[[i]]),
        "conP_%i_adjusted" = as.data.table(av$conP_adj[[i]]),
        "conP_%i_original" = as.data.table(xtabs(formPBase[[i]], data = x)),
        "conP_%i_rel_diff_original" = as.data.table(round(
          100*(av$conP[[i]] - xtabs(formPBase[[i]], data = x)) / av$conP[[i]], 2)),
        "conP_%i_rel_diff_calib" = as.data.table(round(
          100*(av$conP[[i]] - av$conP_adj[[i]]) / av$conP[[i]], 2))
      )
      names(output4) <- sprintf(names(output4), i)
      output <- append(output, output4)
    }
  }


  if(any(names(av) == "conH")){
    output5 <- vector(mode = "list", length = length(av$conH)*5)
    for(i in seq_along(av$conH)){
      j <- ((i-1)*5)+1
      names(output5)[(j):(j+4)] <- c(paste0("conH_",i),paste0("conH_",i,"_adjusted"),paste0("conH_",i,"_original"),
                                     paste0("conH_",i,"_rel_diff_original"),paste0("conH_",i,"_rel_diff_calib"))
      output5[[j]] <- as.data.table(av$conH[[i]]) #conH_i
      output5[[j+1]] <- as.data.table(av$conH_adj[[i]]) #conH_i_adjusted
      output5[[j+2]] <- as.data.table(xtabs(formHBase[[i]],data=x)) #conH_i_original
      output5[[j+3]] <- as.data.table(round(100*(av$conH[[i]]-xtabs(formHBase[[i]],data=x))/av$conH[[i]],2)) #conH_i_rel_diff_original
      output5[[j+4]] <- as.data.table(round(100*(av$conH[[i]]-av$conH_adj[[i]])/av$conH[[i]],2)) #conH_i__rel_diff_calib

    }
    output <- append(output,output5)
  }
  
  return(output)
}

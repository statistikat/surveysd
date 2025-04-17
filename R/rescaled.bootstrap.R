#' @title Draw bootstrap replicates
#'
#' @description Draw bootstrap replicates from survey data using the rescaled
#'   bootstrap for stratified multistage sampling, presented by Preston, J.
#'   (2009).
#'
#' @param dat either data frame or data table containing the survey sample
#' @param REP integer indicating the number of bootstraps to be drawn
#' @param strata string specifying the column name in `dat` that is used for
#'   stratification. For multistage sampling multiple column names can be
#'   specified by `strata=c("strata1","strata2","strata3")` or 
#'   `strata=c("strata1>strata2>strata3")`. See Details for more
#'   information.
#' @param cluster string specifying the column name in `dat` that is used for
#'   clustering. For instance given a household sample the column containing
#'   the household ID should be supplied.
#'   For multistage sampling multiple column names can be specified
#'   by `cluster=c("cluster1","cluster2","cluster3")` or 
#'   `cluster=c("cluster1>cluster2>cluster3")`.
#'   See Details for more information.
#' @param fpc string specifying the column name in `dat` that contains the
#'   number of PSUs at the first stage. For multistage sampling the number of
#'   PSUs at each stage must be specified by `strata=c("fpc1","fpc2","fpc3")` 
#'   or `strata=c("fpc1>fpc2>fpc3")`.
#' @param single.PSU either "merge" or "mean" defining how single PSUs need to
#'   be dealt with. For `single.PSU="merge"` single PSUs at each stage are
#'   merged with the strata or cluster with the next least number of PSUs. If
#'   multiple of those exist one will be select via random draw. For
#'   `single.PSU="mean"` single PSUs will get the mean over all bootstrap
#'   replicates at the stage which did not contain single PSUs.
#' @param return.value either "data", "replicates" and/or "selection" 
#'   specifying the return value of the function. For "data" the survey data is 
#'   returned as class `data.table`, for "replicates" only the bootstrap replicates 
#'   are returned as `data.table`. For "selection" list of data.tables with 
#'   length of `length(strata)` is returned containing 1:`REP` 0-1 columns 
#'   indicating if a PSU was selected for each sampling stage.
#' @param run.input.checks logical, if TRUE the input will be checked before applying
#'   the bootstrap procedure
#' @param already.selected list of data.tables or `NULL` where each data.table contains 
#'   columns in `cluster`, `strata` and additionally 1:`REP` columns containing
#'   0-1 values which indicate if a PSU was selected for each bootstrap replicate.
#'   Each of the data.tables corresponds to one of the sampling stages. First entry
#'   in the list corresponds to the first sampling stage and so on.
#' @param seed integer specifying the seed for the random number generator.
#' 
#' @details For specifying multistage sampling designs the column names in
#' `strata`,`cluster` and `fpc` need to be seperated by ">".\cr
#' For multistage sampling the strings are read from left to right meaning that
#' the first vector entry or column name before the first ">" is taken as the column for
#' stratification/clustering/number of PSUs at the first and the last vector entry
#' or column after
#' the last ">" is taken as the column for stratification/clustering/number of
#' PSUs at the last stage.
#' If for some stages the sample was not stratified or clustered one must
#' specify this by "1" or "I", e.g. `strata=c("strata1","I","strata3")` or 
#' `strata=c("strata1>I>strata3")` if there was
#' no stratification at the second stage or 
#' `cluster=c("cluster1","cluster2","I")` respectively 
#' `cluster=c("cluster1>cluster2>I")`
#' if there were no clusters at the last stage.\cr
#' The number of PSUs at each stage is not calculated internally and must be
#' specified for any sampling design.
#' For single stage sampling using stratification this can usually be done by
#' adding over all sample weights of each PSU by each strata-code.\cr
#' Spaces in each of the strings will be removed, so if column names contain
#' spaces they should be renamed before calling this procedure!\cr 
#' If `already.selected` is supplied the sampling of bootstrap replicates 
#' considers if speficif PSUs have already been selected by a previous survey wave.
#' For a specific `strata` and `cluster` this could lead to more than `floor(n/2)` 
#' records selected. In that case records will be de-selected such that `floor(n/2)` records,
#' with `n` as the total number of records, are selected for each 
#' `strata` and `cluster`. This parameter ist mostly used by [draw.bootstrap] in
#' order to consider the rotation of the sampling units over time.
#'
#' @return returns the complete data set including the bootstrap replicates or
#'   just the bootstrap replicates, depending on `return.value="data"` or
#'   `return.value="replicates"` respectively.
#' @export rescaled.bootstrap
#'
#' @references Preston, J. (2009). Rescaled bootstrap for stratified multistage
#'   sampling. Survey Methodology. 35. 227-234.
#'
#' @author Johannes Gussenbauer, Statistics Austria
#'
#' @examples
#' 
#' library(surveysd)
#' library(data.table)
#' setDTthreads(1)
#' set.seed(1234)
#' eusilc <- demo.eusilc(n = 1,prettyNames = TRUE)
#' 
#' eusilc[,N.households:=uniqueN(hid),by=region]
#' eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=10,strata="region",
#'                                        cluster="hid",fpc="N.households")
#' 
#' eusilc[,new_strata:=paste(region,hsize,sep="_")]
#' eusilc[,N.housholds:=uniqueN(hid),by=new_strata]
#' eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=10,strata=c("new_strata"),
#'                                        cluster="hid",fpc="N.households")
#'
#'



rescaled.bootstrap <- function(
    dat, 
    method = c("Preston", "Rao-Wu"), # preston oder Rao-Wu
    REP = 1000, 
    strata = "DB050>1",   # strata = spalten die Gruppen definieren die sicher abgedeckt werden sollen -> resampling innerhalb der Schicht (Subgruppen repräsentieren)
    cluster = "DB060>DB030",   # cluster = spalten mit natürlichen Gruppen -> resampling ganze cluster (Gruppeneffekte identifizieren)
    fpc = "N.cluster>N.households", # finite population correction: Korrigiert Varianz wenn Stichprobe groß im Verhältnis zur GG ist (Spalten mit der Gesamtzahl der möglichen Gruppen) -> Ansonsten wird Varianz zu groß geschätzt. Ohne FCP wird Varianz so geschätzt als gäbe es unendlich viele Schulen
    single.PSU = c("merge", "mean"), # was tun wenn eine Scicht nur eine Gruppe hat ? mit anderer Schicht zusammenführen vs. mean aus andere bootstrap Stichproben nehmen
    return.value = c("data", "replicates"), # was am ende rauskommt
    run.input.checks = TRUE, # prüft vorab ob alle Eingaben korrekt sind
    already.selected = NULL, # falls man schon eine bootstrap Auswahl hat
    seed = NULL) {
  
  
  InitialOrder <- N <- SINGLE_BOOT_FLAG <- SINGLE_BOOT_FLAG_FINAL <- f <-
    n_prev <- n_draw_prev <- sum_prev <- n_draw <- NULL  # diese ganzen variablen als Null deklarieren
  
  dat <- copy(dat)
  
  # check run.input.checks
  if (!is.logical(run.input.checks)) {
    stop("run.input.checks can only be logical")
  }
  if (run.input.checks) {
    # check input data
    if (!is.data.table(dat) & !is.data.frame(dat)) {
      stop("dat needs to be a data frame or data table")
    } else {
      dat <- data.table(dat)
    }
    # check seed
    if(!is.null(seed)){
      check.input(seed, input.name = "seed", input.length=1, input.type="numeric")
    }
  }
  set.seed        ################################################################################################# @Johannes: Zahl ? #####################################
  
  # prepare input
  removeCols <- c()  #temporär
  if (is.null(cluster)) {
    cluster <- generateRandomName(20, colnames(dat))
    removeCols <- c(removeCols, cluster)
    dat[, c(cluster) := 1:.N]
  }
  # wenn keine cluster/strata explizit angegeben werden dann funktioniert der algo trotzdem auch wenn der Benutzer keine cluster/strata definiert hat
  if (is.null(strata)) {
    strata <- generateRandomName(20, colnames(dat))
    removeCols <- c(removeCols, strata)
    dat[, c(strata) := 1]
  }
  
  
  # check inputs strata, cluster, fpc
  check.input(strata, input.name = "strata", input.type="character")
  check.input(fpc, input.name = "fpc", input.type="character")
  check.input(cluster, input.name = "cluster", input.type="character")
  
  # -> überprüft ob Parameter den richtigen Typ haben
  if(length(strata)==1 && strata %like% ">"){
    strata <- unlist(strsplit(strata, ">"))
    strata <- gsub("\\s", "", strata)
  } # wenn strata ein element hat UND es ein > enthält -> trenne zeichenkette bei >, entferne Leerzeichen
  
  if(length(cluster)==1 && cluster %like% ">"){
    cluster <- unlist(strsplit(cluster, ">"))
    cluster <- gsub("\\s", "", cluster)
  }
  if(length(fpc)==1 && fpc %like% ">"){
    fpc <- unlist(strsplit(fpc, ">"))
    fpc <- gsub("\\s", "", fpc)
  }
  
  single.PSU <- single.PSU[1]  #falls ausversehen mehrere werte übergeben wurden nimm nur den ersten
  # return.value <- return.value[1]
  
  # continue input checks
  if (run.input.checks) {
    
    # check REP
    check.input(REP, input.name = "REP", input.length=1, input.type="numeric",
                decimal.possible = FALSE)
    
    # check method
    method <- match.arg(method)
    
    # Warning for Rao-Wu method
    if (method == "Rao-Wu") {
      warning("The 'Rao-Wu' method should only be used when the first stage sampling is conducted without replacement. Ensure that your sampling design follows this requirement.")
    }
    
    # check design variables
    # check length of strata, cluster, fpc
    if(length(strata) != length(fpc) | length(strata) != length(cluster)){
      stop("strata, cluster, and fpc need to have the same number of arguments",
           " either separated with '>' or be character vectors of the same length")
    }
    
    # cluster can only be "1" or "I" in the final stage
    if(length(cluster)>1 && any(cluster[1:(length(cluster)-1)] %in% c("I","1"))){
      stop("Specifying 'I' or '1' for parameter 'cluster' 
           can only be set in the final sampling stage!")
    }
    
    check.values <- c(strata, cluster, fpc)
    check.values <- check.values[!check.values %in% c("1", "I")]
    check.values <- check.values[!check.values %in% colnames(dat)]
    if (length(check.values) > 0) {
      stop("dat does not contain the column(s): ", check.values)
    }
    
    # check if there are missings in fpc
    if (any(is.na(dat[, .SD, .SDcols = c(fpc)]))) {
      stop("missing values not allowed in fpc")
    }
    
    # check return.value
    if (any(!return.value %in% c("data", "replicates", "selection"))) {
      stop("return.value can only take the values 'data', 'replicates' or 'selection'")
    }
    
    # check single.PSU
    if (is.null(single.PSU) || !single.PSU %in% c("merge", "mean")) {
      warning("single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap",
              " replicates for single PSUs cases will be missing!")
      single.PSU <- FALSE
    }
    
    # check for each stage that PSUs are not in mutiple strata
    for(i in seq_along(strata)){
      if(!strata[i]%in%c("1","I") & !cluster[i]%in%c("1","I")){
        countMultiple <- dat[,uniqueN(strata_i), by=c(cluster[i]),
                             env = list(strata_i = strata[i])]
        countMultiple <- countMultiple[V1>1]
        if(nrow(countMultiple)>0){
          stop("Some sampling units in ",cluster[i]," occur in multiple strata of ",strata[i])
        }
      }
    }
    
  }
  
  # check if variable f, N, n are in data.table
  overwrite.names <- c("f", "N", "n", "n_prev", "n_draw", "n_draw_prev") #namenskonflikte vermeiden, ursprüngliche Daten werden umbenannt in mit Präfix "ORIGINAL_"
  overwrite.names <- overwrite.names[overwrite.names %in% colnames(dat)]
  if (length(overwrite.names) > 0) {
    overwrite.names.new <- paste0("ORIGINAL_", overwrite.names)
    setnames(dat, overwrite.names, overwrite.names.new)
    
    strata[strata %in% overwrite.names] <-
      overwrite.names.new[strata %in% overwrite.names]
    cluster[cluster %in% overwrite.names] <-
      overwrite.names.new[cluster %in% overwrite.names]
    fpc[fpc %in% overwrite.names] <-
      overwrite.names.new[fpc %in% overwrite.names]
  }
  
  # set index for data to return dat in correct order
  dat[, InitialOrder := .I]  # neue spalte in dat mit Namen InitialOrder
  
  # calculate bootstrap replicates
  stages <- length(strata) # Anzahl an Stufen basierend auf Strata (verschiedene Schichten)
  print("stages: ")
  print(stages)
  
  n <- nrow(dat)
  # define values for final calculation
  n.calc <- matrix(0, nrow = n, ncol = stages) # n = Stichprobengröße
  N.calc <- matrix(0, nrow = n, ncol = stages) #  N = Gesamtpopulation
  n_draw.calc <- matrix(0, nrow = n, ncol = stages)  # n_draw = Ziehungszahlen in den verschiedenen Stufen
  delta.calc <- array(0, dim = c(n, stages, REP)) # um Änderungen von delta über alle stufen und reps hinweg zu speichern
  
  print("n.calc:")
  print(head(n.calc))
  
  print("N.calc:")
  print(head(N.calc)) 
  
  print("n_draw.calc:")
  print(head(n_draw.calc)) 
  
  print("delta.calc:")
  print(head(delta.calc))
  
  
  
  delta_selection <- list()
  
  for (i in 1:stages) {  # iteriere über jede stufe im mehrstufigen Prozess
    # iterate over sampling stages
    # select units in each sampling stage and
    # calucalte factors to finally generate
    # boostrap replicates
    
    # define by.val
    by.val <- strata[i]  # by.val bestimmt die Gruppierung für die aktuelle Stufe
    if (i > 1) {
      by.val <- c(strata[1:(i - 1)], cluster[1:(i - 1)], strata[i]) # kombiniert strata der aktuellen und vorherigen stufe
    }
    by.val <- by.val[!by.val %in% c("1", "I")] # platzhater herausfiltern
    
    # define cluster value
    clust.val <- cluster[i]
    if (clust.val %in% c("1", "I")) {
      clust.val <- paste0("ID_help_", i)  # erstellt hilfsids wenn kein echtes Cluster definiert ist, Sampling erfolgt auf Individualebene (weil immer eine Cluster-Definition benötigt wird)
      dat[, c(clust.val) := .I, by = c(by.val)]
    }
    
    # check of fpc[i] is well defined
    # fpc[i] can not have different values per each by.val
    check.fpc <- dat[,uniqueN(fpc_i), by=c(by.val), 
                     env = list(fpc_i = fpc[i])]  # prüfen ob fpc werte innerhalb jeder stratum-Cluster-kombination konsistent sind
    check.fpc <- check.fpc[V1>1]
    if (nrow(check.fpc) > 0) {
      stop("values in ", fpc[i], " do vary in some strata-cluster ",
           "combinations at sampling stage ", i)
    }
    
    singles <- dat[fpc_i>1, sum(!duplicated(clust.val)), by=c(by.val),   # findet gruppen in denen nur eine einzige psu vorhanden ist -> das ist problematisch für bootstrap (mit nur einem element kann nicht geresamplet werden)
                   env = list(fpc_i = fpc[i],
                              clust.val = clust.val)]
    singles <- singles[V1==1]
    if (nrow(singles) > 0) {
      # if singel.PSU=="merge" change the coding of the current stage of the
      #   single PSU
      # to the coding in the same subgroup, according the next higher stage,
      #   with the smallest number of PSUs
      # if multiple smallest exist choose one per random draw
      higher.stages <- c(strata[1:(i - 1)], cluster[1:(i - 1)])
      by.val.tail <- tail(by.val,1)
      firstStage <- length(higher.stages) == 0
      if (firstStage) {
        higher.stages <- by.val
      }
      singles <- unique(subset(singles, select = higher.stages))
      
      ################ SINGLE ANFANG ################ ##############################################################################################
      
      if (method == "Preston") {
        # wenn merge angegeben wird: Diese Gruppen werden mit einer ähnlichen (übergeordneter Struktur, z.B. Region) zusammengelegt, sodass mehr als 1 PSU entsteht.
        if (single.PSU == "merge") {
          # save original PSU coding and fpc values to replace changed values
          #   bevore returning the data.table
          if (return.value == "data") {
            by.val.orig <- paste0(by.val.tail, "_ORIGINALSINGLES")
            fpc_orig <- paste0(fpc[i], "_ORIGINALSINGLES")
            
            set(dat, i = by.val.orig, value = dat[[by.val.tail]])
            set(dat, i = fpc_orig, value = dat[[fpc[i]]])
          }
          
          setkeyv(dat, higher.stages)
          next.PSU <- dat[singles, .(N = sum(!duplicated(clust.val))), by = c(by.val),
                          on = c(higher.stages),
                          env = list(clust.val = clust.val)]
          
          if(nrow(next.PSU)==1){
            stop("Only 1 single PSU in first sampling stage! Cannot combine single PSU.\nPlease manually combine sampling STRATA to eliminate single PSU!")
          }
          
          new.var <- paste0(tail(by.val, 1), "_NEWVAR")
          
          next.PSU[,c(new.var) := next_minimum(N, by.val.tail), by = c(higher.stages),
                   env = list(by.val.tail = by.val.tail)]
          if(any(is.na(next.PSU[[new.var]]))){
            if(firstStage){
              next.PUS[is.na(new.var_col), c(new.var) := head(higher.stages, 1),      ################################################################################################# @Johannes: PSU ? #####################################
                       env = list(new.var_col = new_var,
                                  higher.stages = higher.stages)]
            }else{
              next.PUS[is.na(new.var_col), c(new.var) := head(by.val.tail,1),       ################################################################################################# @Johannes: PSU ? #####################################
                       env = list(new.var_col = new_var,
                                  by.val.tail = by.val.tail)]
            }
          }
          
          # select singel PSUs and PSUs to join with
          next.PSU <- next.PSU[N == 1 | new.var == by.val.tail,
                               env = list( new.var = new.var,
                                           by.val.tail = by.val.tail)]
          dat <- merge(dat, next.PSU[, mget(c(by.val, new.var))],
                       by = c(by.val), all.x = TRUE)
          
          # sum over margins
          fpc.i_ADD <- paste0(fpc[i], "_ADD")
          dat[!is.na(new.var), c(fpc.i_ADD) := 
                sum(fpc_i[!duplicated(by.val.tail)]), by = c(new.var),
              env = list(new.var = new.var,
                         fpc_i = fpc[i],
                         by.val.tail = by.val.tail)]
          
          # assign to new group
          dat[!is.na(new.var), c(by.val.tail) := new.var, 
              env = list(new.var = new.var)]
          dat[ c(fpc[i]) := fpc_i[is.na(new.var)][1], by=c(by.val),
               env = list(new.var = new.var)]
          dat[!is.na(new.var), c(fpc[i]) := fpc.i_ADD,
              env = list(new.var = new.var,
                         fpc.i_ADD = fpc.i_ADD)]
          
          dat[, c(new.var, paste0(fpc[i], "_ADD")) := NULL]
        } 
        
        # wenn bei ingle = mean: Markiert die betroffenen Beobachtungen mit einem SINGLE_BOOT_FLAG, um später im Bootstrap-Prozess einfach den Durchschnitt anderer Replicates zu nehmen.
        else if (single.PSU == "mean") {
          # if single.PSU="mean" flag the observation as well as the all the observations in the higher group
          # Vermeidung vin NA Wetren in Bootstrap Replikaten wenn eine Schicht nur eine PSU hat
          # Die Methode ersetzt die Werte dieser einzelnen PSU durch den Durchschnitt der Bootstrap-Replikate aus den anderen PSUs.
          singles[, SINGLE_BOOT_FLAG := paste(higher.stages, .GRP, sep = "-"), # singles werden markiert -> 
                  by = c(higher.stages)]
          
          dat <- merge(dat, singles, by = c(higher.stages), all.x = TRUE)
          if (!"SINGLE_BOOT_FLAG_FINAL" %in% colnames(dat)) {
            dat[, SINGLE_BOOT_FLAG_FINAL := SINGLE_BOOT_FLAG]
          } else {
            dat[is.na(SINGLE_BOOT_FLAG_FINAL),
                SINGLE_BOOT_FLAG_FINAL := SINGLE_BOOT_FLAG]
          }
          dat[, SINGLE_BOOT_FLAG := NULL]
          
        } 
        
        else {
          message("Single PSUs detected at the following stages:\n")
          dat.print <- dat[,sum(!duplicated(clust.val)), by=c(by.val),
                           env = list(clust.val = clust.val)]
          dat.print <- dat.print[V1==1,.SD,.SDcols=c(by.val)]
        }
      } else if (method == "Rao-Wu") {
        if (single.PSU == "merge") {
          # save original PSU coding and fpc values to replace changed values
          #   bevore returning the data.table
          if (return.value == "data") {
            by.val.orig <- paste0(by.val.tail, "_ORIGINALSINGLES")
            fpc_orig <- paste0(fpc[i], "_ORIGINALSINGLES")
            
            set(dat, i = by.val.orig, value = dat[[by.val.tail]])
            set(dat, i = fpc_orig, value = dat[[fpc[i]]])
          }
          
          setkeyv(dat, higher.stages)
          next.PSU <- dat[singles, .(N = sum(!duplicated(clust.val))), by = c(by.val),
                          on = c(higher.stages),
                          env = list(clust.val = clust.val)]
          
          if(nrow(next.PSU)==1){
            warning("Only 1 single PSU in first sampling stage! Results may be unreliable. Consider manually combining strata.")
          }
          
          new.var <- paste0(tail(by.val, 1), "_NEWVAR")
          
          next.PSU[,c(new.var) := next_minimum(N, by.val.tail), by = c(higher.stages),
                   env = list(by.val.tail = by.val.tail)]
          if(any(is.na(next.PSU[[new.var]]))){
            if(firstStage){
              next.PUS[is.na(new.var_col), c(new.var) := head(higher.stages, 1),
                       env = list(new.var_col = new_var,
                                  higher.stages = higher.stages)]
            }else{
              next.PUS[is.na(new.var_col), c(new.var) := head(by.val.tail,1),
                       env = list(new.var_col = new_var,
                                  by.val.tail = by.val.tail)]
            }
          }
          
          # select singel PSUs and PSUs to join with
          next.PSU <- next.PSU[N == 1 | new.var == by.val.tail,
                               env = list( new.var = new.var,
                                           by.val.tail = by.val.tail)]
          dat <- merge(dat, next.PSU[, mget(c(by.val, new.var))],
                       by = c(by.val), all.x = TRUE)
          
          # sum over margins
          fpc.i_ADD <- paste0(fpc[i], "_ADD")  
          dat[!is.na(new.var), c(fpc.i_ADD) := 
                sum(fpc_i[!duplicated(by.val.tail)]), by = c(new.var),
              env = list(new.var = new.var,
                         fpc_i = fpc[i],
                         by.val.tail = by.val.tail)]
          
          # assign to new group
          dat[!is.na(new.var), c(by.val.tail) := new.var, 
              env = list(new.var = new.var)]
          dat[ c(fpc[i]) := fpc_i[is.na(new.var)][1], by=c(by.val),
               env = list(new.var = new.var)]
          dat[!is.na(new.var), c(fpc[i]) := fpc.i_ADD,
              env = list(new.var = new.var,
                         fpc.i_ADD = fpc.i_ADD)]
          
          dat[, c(new.var, paste0(fpc[i], "_ADD")) := NULL]
        } 
        
        # wenn bei ingle = mean: Markiert die betroffenen Beobachtungen mit einem SINGLE_BOOT_FLAG, um später im Bootstrap-Prozess einfach den Durchschnitt anderer Replicates zu nehmen.
        else if (single.PSU == "mean") {
          
          # Sicherstellen, dass die Single PSU mindestens einmal gezogen wird
          dat[!is.na(SINGLE_BOOT_FLAG_FINAL), c(bootRep) := lapply( 
            .SD,
            function(z) {
              if (all(is.na(z))) { #all(is.na(z)): Prüft, ob alle Werte in der aktuellen Spalte z NA sind.
                avg <- mean(z, na.rm = TRUE)
                #  Wenn alle Werte einer Replikat-Spalte NA sind, wird ein Standardwert (z. B. 0) verwendet -> Rao-Wu-Methode ist darauf ausgelegt, mit unsicheren Daten umzugehen
                if (is.na(avg)) avg <- 0 
                return(rep(avg, length(z)))  # Berechnet den Durchschnitt der Werte in z, ignoriert dabei NA-Werte 
              } else {
                return(z) # Wenn nicht alle Werte NA sind, wird der Originalvektor z unverändert zurückgegeben.
              }
            }
          ), by = SINGLE_BOOT_FLAG_FINAL, .SDcols = c(bootRep)] #Das bedeutet, dass alle Zeilen, die zur gleichen Single PSU gehören, zusammen verarbeitet werden.
        }
        else {
          message("Single PSUs detected at the following stages:\n")
          dat.print <- dat[,sum(!duplicated(clust.val)), by=c(by.val),
                           env = list(clust.val = clust.val)]
          dat.print <- dat.print[V1==1,.SD,.SDcols=c(by.val)]
          
          #Falls alle Werte in einer Replikat-Spalte NA sind, wird der Durchschnitt der anderen Replikate verwendet.
        }
        
      }
    }
    ################ SINGLE ENDE ################ ##############################################################################################
    
    # get Stage
    if (i == 1) { # dati: enthält alle PSU (PSU = die ursprünglichen Cluster oder Gruppen in einer komplexen Stichprobe) mit zugehörigen FPC 
      dati <- dat[,.(N=fpc_i[1], # Populationsgröße in dieser Gruppe
                     clust.val = unique(clust.val), # Cluster-ID (jede PSU)
                     f = 1, # Samplingfaktor (Startwert) = 1 da noch keine Reskalierung
                     n_prev = 1, n_draw_prev = 1, # Startwerte für vorige Stufe (irrelevant in Stufe 1)
                     sum_prev = 1),  # Zusatz für spätere Berechnungen (wird unten nicht verwendet)
                  by = c(by.val),
                  env = list(fpc_i = fpc[i],
                             clust.val = clust.val)]
      
    } else { # In Stufe 2 oder höher: Jetzt werden die Infos aus vorheriger Ziehung übernommen
      dati <- dat[,.(N=fpc_i[1], 
                     clust.val = unique(clust.val),
                     f = f[1], # Samplingfaktor aus vorheriger Stufe
                     n_prev = n_prev[1], # gezogene Anzahl auf vorheriger Stufe
                     n_draw_prev = n_draw_prev, # gezogene PSU vorher
                     sum_prev = sum_prev[1]), # vorherige Summenkomponente (für Varianzformeln)
                  by = c(by.val),
                  env = list(fpc_i = fpc[i],
                             clust.val = clust.val)]
    }
    
    # dati = Tabelle mit einer Zeile pro PSU, inklusive N, clust.val, f (sampling faktor), n_prev (Zahl gezogener PSU in vorheriger Stufe), n_draw_prev (gezogene PSU in vorheriger Stufe), sum_prev (Hilfswert), n (Anzahl PSU in Gruppe), n_draw (Anzahl PSU die neu gezogen werden sollen)
    
    deltai <- paste0("delta_", i, "_", 1:REP) #erstellt eine Variable, die eine Liste von Namen für zukünftige Delta-Werte erstellt
    dati[, n := .N, by = c(by.val)]  # n = Anzahl PSU 
    
    print(paste("dati:"))
    print(head(dati))
    print(paste("Type of dati:"))
    print(class(dati))
    
    print(paste("deltai:"))
    print(head(deltai))
    print(paste("Type of deltai:"))
    print(class(deltai))
    
    # determin number of psu to be drawn
    # dati[, n_draw := select.nstar(
    #   n[1], N[1], f[1], n_prev[1], n_draw_prev[1], sum_prev = NULL,
    #   new.method = new.method), by = c(by.val)]
    if (method == "Preston") {
      dati[,n_draw:=floor(n/2),by = c(by.val)]      # n_draw= hälfte der vorhandenen PSU
    } else {
      dati[, n_draw := n - 1, by = c(by.val)] # steht so im EZB Bericht
    }

    
    if (nrow(dati[n_draw == 0]) > 0) {  # wenn eine Gruppe kein PSU zum ziehen hat
      stop("Resampling 0 PSUs should not be possible! Please report bug in ",
           "https://github.com/statistikat/surveysd")
    }
    
    # do bootstrap for i-th stage
    if(!is.null(already.selected)){ # prüfe ob in der aktuellen stufe i schon gezogene PSUs vorhanden sind
      # if already.selected was supplied consider already selected units
      dati_selected <- already.selected[[i]]
      dati[dati_selected,c(deltai):=mget(deltai),on=c(by.val,clust.val)]
      
      dati[,c(deltai):=lapply(.SD,function(delta,n,n_draw){ # funktion OHNE zurücklegen um bootstrap stichprobe zu ziehen -> für jedes deltai (sammlung von delta werten)
        
        if (method == "Rao-Wu") {  # ***Rao-Wu
          draw.with.replacement(n[1], n_draw[1], delta=delta) # ***Rao-Wu
        } else if (method == "Preston"){
          draw.without.replacement(n[1],n_draw[1], delta=delta)
          print("n[1]:")
          print(n[1])
          
          print("n_draw[1]:")
          print(n_draw)
        }
      },n=n,n_draw=n_draw),by=c(by.val),.SDcols=c(deltai)] # überprüfe ob die Gesamtzahl der gezogenen PSUs für jede Gruppe mit dem erwarteten Wert übereinstimmt
      
      dati_check <- dati[,lapply(.SD,function(z,n_draw){
        sum(z)==n_draw[1]
      },n_draw=n_draw),by=c(by.val),.SDcols=c(deltai)]
      
      if(any(!dati_check[,.SD,.SDcols=c(deltai)])){
        stop("Wrong number of units selected! Please report bug in ",
             "https://github.com/statistikat/surveysd")
      }
    }else{ # wenn es noch keine PSUs gibt -> einfacher bootstrap ohne Wiederholung
      # do simple sampling without replacement
      
      # ***Rao-Wu
      if (method == "Rao-Wu") {
        dati[, c(deltai) := as.data.table(
          replicate(REP, draw.with.replacement(n[1], n_draw[1]),simplify = FALSE)
        ), by = c(by.val)] # ***Rao-Wu
        
      } else if (method == "Preston") {
        print("n[1]:")
        print(n[1])
        
        print("n_draw[1]:")
        print(n_draw)
        
        dati[, c(deltai) := as.data.table(
          replicate(REP, draw.without.replacement(n[1], n_draw[1]), # Funktion "replicate" wiederholt den Ziehvorgang REP mal
                    simplify = FALSE)),
          by = c(by.val)]
      }
    }
    
    
    # merge with data
    if(i>1){ # wenn nicht erste stufe
      dat[,c("f","n_prev","n_draw_prev","sum_prev"):=NULL]  # alte Variablen aus dat entfernen
    }
    
    dat <- merge(dat,dati,by=c(by.val,clust.val)) # dann basierend auf variablen by.val und clust.val die tabelle dati und dat zusammenführen
    setorder(dat,InitialOrder) # Reihenfolge ist initial order
    
    # prepare output for return.value %in% "selection"
    delta_selection <- c(delta_selection,
                         list(dat[, .SD, .SDcols = c(cluster, strata, deltai)]))  # enthält  für jede Stufe die Informationen zu den gezogenen PSUs und deren Cluster
    
    # extract information from data.table and remove again from data table
    # (less memory intensive)
    # only matrices and arrays needed for final calculation
    n.calc[, i] <- dat[, n]
    N.calc[, i] <- dat[, N]
    n_draw.calc[, i] <- dat[, n_draw]
    delta.calc[, i, ] <- as.matrix(dat[, mget(deltai)])
    
    # aktualisiere Variablen für die nächste Stufe
    # if (method == "Rao-Wu") {
    #   #dat[, f := sqrt((N - n)/(N - 1)) * (n_draw / n - 1) + 1]
    #   dat[, f := n / (n - 1)]  
    #   
    # } else if (method == "Preston") {
    #   dat[, f := n / N * f]
    # }
    dat[, f := n / N * f]    ########################################################################################## Richtig?
    dat[, n_prev := n * n_prev]
    dat[, n_draw_prev := n_draw * n_draw_prev]
    
    dat[, c("n", "N", deltai, "n_draw") := NULL] # entfernen
    
    rm(dati);gc() # entfernen,Speicher durch den Garbage Collector bereinigen
  }
  
  print("sampling fraction f: ")
  print (head(dat[,f]))
  # remove names
  dat[,c("f","n_prev","n_draw_prev","sum_prev"):=NULL]
  
  # rename delta_selection for output
  names(delta_selection) <- paste0("SamplingStage",1:stages)
  
  # calculate bootstrap replicate values
  bootRep <- paste0("bootRep", 1:REP) # für jede bootstrap wiederholung wird neue variable erstellt
  

  print("Type of n.calc:")
  print(class(n.calc))
  
  print("Type of N.calc:")
  print(class(N.calc))
  
  print("Type of n_draw.calc:")
  print(class(n_draw.calc))
  
  print("Type of delta.calc:")
  print(class(delta.calc))
  
  dat[, c(bootRep) := as.data.table(calc.replicate( # rufe function calc.replicate auf -> berechne bootstrap replicate für jede Wiederholung
    method = method,
    n = n.calc, 
    N = N.calc, 
    n_draw = n_draw.calc, 
    delta = delta.calc))] 
  
  if (single.PSU == "mean") {
    dat[!is.na(SINGLE_BOOT_FLAG_FINAL), c(bootRep) := lapply(  # wenn single.PSU mean ist und SINGLE_BOOT_FLAG_FINAL nicht na -> durchschnitt über alle bootstrap wiederholungen
      .SD,
      function(z) {
        mean(z, na.rm = TRUE)
      }
    ), by = SINGLE_BOOT_FLAG_FINAL, .SDcols = c(bootRep)]
  }
  
  print("dat: ")
  print (head(dat))
  
  setkey(dat, InitialOrder)  # reihenfolge der daten wird widerhergestellt
  if (length(removeCols) > 0) {
    dat[, c(removeCols) := NULL]
  }
  
  if ("data" %in% return.value) {
    # get original values for PSUs and fpc - if singles PSUs have been detected
    #   and merged
    if (single.PSU == "merge") {  # wenn merge: prüfe ob Spaltennamen vorhanden sind, die auf _ORIGINALSINGLES enden. Diese Spalten werden entfernt, und die ursprünglichen Namen werden wiederhergestellt, indem _ORIGINALSINGLES entfernt wird.
      c.names <- colnames(dat)
      c.names <- c.names[grepl("_ORIGINALSINGLES", c.names)]
      if (length(c.names) > 0) {
        drop.names <- gsub("_ORIGINALSINGLES", "", c.names)
        dat[, c(drop.names) := NULL]
        setnames(dat, c(c.names), drop.names)
      }
    }
    dat[, c("InitialOrder") := NULL]  # initial order entfernen
    # get original col values
    if (length(overwrite.names) > 0) {
      setnames(dat, overwrite.names.new, overwrite.names)
    }
    
    output_bootstrap <- list(data = dat)
  } else if ("replicates" %in% return.value) { #Nur die Bootstrap-Replikate (die Spalten bootRep1, bootRep2, ..., bootRep[REP]) werden zurückgegeben.
    output_bootstrap <- list(replicates = dat[, mget(bootRep)]) 
  }
  
  if("selection" %in% return.value){ # Die delta_selection-Liste (die die gezogenen PSUs enthält) wird zur Ausgabe hinzugefügt.
    if(length(return.value)>1){
      output_bootstrap <- c(output_bootstrap, list(selection = delta_selection))
    }else{
      output_bootstrap <- list(selection = delta_selection)
    }
  }
  
  if(length(output_bootstrap) == 1){
    output_bootstrap <- output_bootstrap[[1]]
  }
  return(output_bootstrap)
}

select.nstar <- function(n, N, f, n_prev, n_draw_prev, lambda_prev,
                         sum_prev = NULL, new.method) {
  if (n == 1) {
    # if only a single unit in strata
    # return missing
    # if single units are
    # not treated missing values are returned
    return(1L)
  }
  
  if (!is.null(sum_prev)) {
    n_draw <- (sum_prev) ^ 2 / ((1 - (n / N)) * n_prev * f +
                                  (sum_prev) ^ 2) * n
    n_draw <- floor(n_draw)
  } else {
    if (new.method) {
      n_draw <- (n * n_draw_prev) / (n_prev * f * (1 - n / N) + n_draw_prev)
      n_draw <- floor(n_draw)
    } else {
      n_draw <- floor(n / 2)
    }
  }
  
  return(n_draw) # n_draw = Größe der Unterstichprobe
}



# Stichprobenziehung OHNE Zurücklegen (Preston-Methode)
draw.without.replacement <- function(n, n_draw, delta = NULL) {
  # n = Geamtpopulation
  # n_draw = Ziehung einer bestimmten Anzahl von Einheiten (n_draw)
  # delta: 1= ausgewählt, 0 = nicht ausgewählt, NA = noch nicht gezogen
  
  # if no units have been selected prior
  if(is.null(delta)){
    delta <- rep(c(1.0, 0.0), c(n_draw, n - n_draw)) # wenn delta nicht übergeben wird (es wurden noch keine Einheiten erstellt), wird er hier erstellt. Alle Einheiten werden mit whkeit von 1.0 oder o.o ausgewählt
    if (length(delta) > 1) {
      delta <- sample(delta)
    }
  }else{
    # if units have already been selected
    n_selected <- sum(delta,na.rm=TRUE)  # wenn delta existiert dh. bereits Einheiten ausgewählt wurden wird delta angepasst: n_selected = Anzahl bisher ausgewählten Einheiten
    delta_new_unit <- is.na(delta) # TRUE für nicht gezogene Einheiten, FALSE für bereits gezogene Einheiten
    
    # check if only 1 unit needs to be selected
    if(sum(delta_new_unit)==1){
      # deselect one element in delta
      # so that resulting element will have a selection probability
      if(n_selected<n_draw){
        # set one entry in delta==0 to NA
        delta <- set.random2NA(delta,value2NA = 0)  # nicht gezogene Einheit (0) wird auf NA gesetzt
      }else{
        # set one entry in delta==1 to NA
        delta <- set.random2NA(delta,value2NA = 1) # bereits gezogene Einheit (1) word auf NA gesetzt
      }
      
      n_selected <- sum(delta,na.rm=TRUE)
      delta_new_unit <- is.na(delta)
    }
    
    if(n_draw<=n_selected){
      # if more entries have been selected than actually possible
      # or n_draw units have already been selected
      # due to already selected units
      # set some units with delta==1 to 0
      # such that entries with delta NA have selection probability >0
      add1 <- as.numeric(sum(delta_new_unit)>1) # Wenn noch mehr Einheiten ausgewählt werden müssen, als aktuell als NA markiert sind, werden einige der 0-Werte in delta zufällig zu 1 geändert,
      nChanges <- n_selected - n_draw + add1
      
      if(nChanges == sum(delta_new_unit==FALSE)){ 
        delta <- set.random2NA(delta, value2NA=1,
                               nChanges = n_selected - n_draw + add1) 
      }else{
        delta <- change.random.value(delta, changeVal=1,
                                     nChanges = n_selected - n_draw + add1)
      }
      delta_new_unit <- is.na(delta)
      n_selected <- sum(delta,na.rm=TRUE)
    }
    
    if((n_draw-n_selected) >= sum(delta_new_unit)){
      # if new units will all be selected
      # or more units need to be selected than possible
      # set random entries from 0 to 1
      add1 <- as.numeric(sum(delta_new_unit)>1)
      nChanges <- n_draw-n_selected-sum(delta_new_unit) + add1
      delta <- change.random.value(delta, changeVal=0,
                                   nChanges = nChanges)
      n_selected <- sum(delta,na.rm=TRUE)
      delta_new_unit <- is.na(delta)
    }
    # sicherstellen dass genau die Anzahl an Einheiten die benötigt wird in der Stichprobe sind
    # Zu viele Einheiten: von "ausgewäht" zu "nicht ausgewählt"
    # Wenn zu wenig: es werden weitere Einheiten aus denen, dei noch nicht ausgewählt wurden hinzugefügt
    
    n_draw <- n_draw - n_selected
    delta_rest <- rep(c(1.0, 0.0), c(n_draw, sum(delta_new_unit)-n_draw))
    if (length(delta_rest) > 1) {
      delta_rest <- sample(delta_rest)
    }
    delta[delta_new_unit==TRUE] <- delta_rest
    
  }
  return(delta)
}


# Stichprobenziehung MIT Zurücklegen (Rao-Wu-Methode)
draw.with.replacement <- function(n, n_draw, delta = NULL) { 
  # Funktion für Stichprobenziehung MIT Zurücklegen nach dem Rao-Wu-Verfahren
  # Argumente:
  # n = Gesamtpopulation
  # n_draw = Anzahl der zu ziehenden Einheiten
  # delta = aktueller Status, welche Einheiten ausgewählt wurden (1 = ausgewählt, 0 = nicht ausgewählt)
  
  # Initialisiere `r_hi_star` (Häufigkeiten der Ziehungen)
  if (is.null(delta)) {
    delta <- rep(0, n)  # Setze alle Einheiten auf 0 (noch nicht gezogen)
  }
  # Anzahl der bereits ausgewählten Einheiten
  n_selected <- sum(delta, na.rm = TRUE)
  
  # 1. Ziehe zusätzliche Einheiten, falls weniger als n_draw ausgewählt sind
  if (n_selected < n_draw) {
    n_needed <- n_draw - n_selected # Berechnung der fehlenden Anzahl von Ziehungen
    additional_draws <- sample(1:n, size = n_needed, replace = TRUE)  # Ziehungen mit Zurücklegen
    for (i in additional_draws) {
      delta[i] <- delta[i] + 1     # Aktualisiere delta: Erhöhe die Häufigkeit der gezogenen Einheiten
    }
  }
  
  # 2. Reduziere überschüssige Einheiten, falls mehr als n_draw ausgewählt sind
  if (n_selected > n_draw) {
    over_selected <- n_selected - n_draw # Berechnung der überschüssigen Ziehungen
    for (i in which(delta > 0)[1:over_selected]) {
      delta[i] <- delta[i] - 1     # Entferne überschüssige Ziehungen gleichmäßig
    }
  }
  
  # 3. Behandle NA-Werte in delta, falls vorhanden
  if (any(is.na(delta))) {
    delta[is.na(delta)] <- 0     # Setze NA-Werte auf 0 (nicht gezogen)
  }
  
  return(delta)
}



# change.random.value <- function(delta, changeVal = 0, nChanges = 1, method = "Preston") {
#   
#   if (method == "Preston") {
#     # Logik für die Preston-Methode (ohne Zurücklegen)
#     changeVal2 <- fifelse(changeVal == 0, 1, 0)  # Wechsel zwischen 0 und 1
#     set2NA <- which(delta == changeVal & !is.na(delta))  # Finde Positionen mit changeVal
#     
#     if (length(set2NA) > 1) {
#       # Zufällige Auswahl von nChanges Positionen
#       set2NA <- sample(set2NA, min(nChanges, length(set2NA)))
#     }
#     delta[set2NA] <- changeVal2  # Werte ändern
#     return(delta)
#     
#   } else if (method == "Rao-Wu") {
#     # Logik für die Rao-Wu-Methode (mit Zurücklegen)
#     set2Change <- which(!is.na(delta))  # Alle Positionen, die nicht NA sind
#     
#     if (length(set2Change) > 1) {
#       # Zufällige Auswahl von nChanges Positionen
#       set2Change <- sample(set2Change, min(nChanges, length(set2Change)))
#     }
#     delta[set2Change] <- delta[set2Change] + 1  # Häufigkeit erhöhen
#     return(delta)
#     
#   } else {
#     stop("Invalid method. Please use 'Preston' or 'Rao-Wu'.")  # Fehler für ungültige Methode
#   }
# }

change.random.value <- function(delta,changeVal=0,nChanges=1){
  
  changeVal2 <- fifelse(changeVal==0,1,0) # Wechsel zwischen 0 und 1
  set2NA <- which(delta==changeVal & !is.na(delta))
  if(length(set2NA)>1){
    # Zufällige Auswahl von nChanges Positionen
    set2NA <- sample(set2NA,min(nChanges,length(set2NA)))
  }
  delta[set2NA] <- changeVal2
  return(delta)
}



set.random2NA <- function(delta,value2NA=0,nChanges=1){
  
  set2NA <- which(delta==value2NA & !is.na(delta))
  if(length(set2NA)>1){
    set2NA <- sample(set2NA,min(nChanges,length(set2NA)))
  }
  delta[set2NA] <- NA
  return(delta)
}


calc.replicate <- function(n, N, n_draw, delta , method = "Preston") {
  
  # n = n.calc, N = N.calc, n_draw = n_draw.calc, delta = delta.calc
  # n/n.calc, N/N.calc, n_draw/n_draw.calc sind Matrizen
  # delta_1_1, delta_1_2, etc. sind Resampling-Indikatoren -> 0 = nicht ausgewählt, 1 = ausgewählt
  
  p <- ncol(n)  #  Anzahl der Spalten in der Matrix n,
  # n_draw <- trunc(n/2)
  # n_draw <- floor(n/(2-rowCumprods(n/N))-1)
  # n_draw[n_draw<1] <- 1
  dimdelta <- dim(delta)
  for (i in 1:p) {
    if (method == "Preston") {
      if (i == 1) { # 
        lambda <- sqrt(n_draw[, 1] *
                         (1 - n[, 1] / N[, 1]) / (n[, 1] - n_draw[, 1]))
        rep_out <- 1 - lambda + lambda * n[, i] / n_draw[, i] * delta[, i, ] # rep_out ist das Ergebnis der Berechnung für die Replikate -> Gewichtung der delta-Werte für den ersten Schritt.
      } else if (i == 2) {
        lambda <- (1 - n[, i] / N[, i]) / (n[, i] - n_draw[, i])
        lambda <- sqrt((n[, i - 1] / N[, i - 1]) * n_draw[, i] * lambda)
        rep_out <- rep_out + lambda *                                        # rep_out wird dann um die Berechnung für die zweite Stufe ergänzt, wobei auch die Abweichungen von der ersten Stufe (delta[, i - 1, ]) berücksichtigt werden.
          (sqrt(n[, i - 1] / n_draw[, i - 1]) * delta[, i - 1, ]) *
          (n[, i] / n_draw[, i] * delta[, i, ] - 1)
      } else {
        lambda <- (1 - n[, i] / N[, i]) / (n[, i] - n_draw[, i])
        lambda <- sqrt(rowProds(n[, 1:(i - 1)] / N[, 1:(i - 1)]) *
                         n_draw[, i] * lambda)
        prod_val <- matrix(0, ncol = dimdelta[3], nrow = dimdelta[1])
        for (r in 1:dimdelta[3]) {
          prod_val[, r] <- rowProds(sqrt(n[, 1:(i - 1)] / n_draw[, 1:(i - 1)]) *
                                      delta[, 1:(i - 1), r])
          
          print("prod_val:")
          print(prod_val)
        }
        # rep_out <- rep_out + lambda*rowProds(sqrt(n[,1:(i-1)]/
        #   n_draw[,1:(i-1)])*delta[,1:(i-1),]) * (n[,i]/n_draw[,i]*delta[,i,]-1)
        rep_out <- rep_out + lambda * prod_val * (n[, i] / n_draw[, i] *      # Schließlich wird rep_out für jede zusätzliche Stufe unter Berücksichtigung der vorherigen Stufen aktualisiert.
                                                    delta[, i, ] - 1)
      }
      print("n:")
      print(head(n))
      print("N:")
      print(head(N))
      
      print("delta:")
      print(head(delta))
      
      print("lambda:")
      print(lambda)
      
      
    } else if (method == "Rao-Wu") {
      # darf für merhstufige designs angewendet werden solange sampling fraction in der ersten stufe klein ist
      # Die Varianz, die durch die Auswahl der PSUs entsteht, dominiert die Gesamtvarianz.Die zusätzliche Unsicherheit durch die Auswahl innerhalb der PSUs (also auf Stufe 2 oder 3) ist vergleichsweise gering.
      
      # Für den ersten Schritt (FPC-Korrektur)
      n_h <- n[, i]               # Anzahl gezogener PSUs (Primary Sampling Units)
      m_h <- n_h - 1              # Für Resampling einen weglassen
      f_h <- n_h / N[, i]         # Stichprobenanteil
      print(f_h)
      
      # if (any(f_h > 0.1)) {  
      #   stop("Sampling Fraction is too big for Rao-Wu, choose preston instead") ################################################################################################# @Johannes: Welche Sampling Fraction ? (Recherche sagt maximales f von 0.05 - 0.1) #####################################
      # }
      
      w_hi <- N[, i] / n_h        # Designgewichte (Inverse der Auswahlwahrscheinlichkeit)
      
      # Skalierungsfaktor λ, lambda bestimmt die Varianzkorrektur
      lambda <- sqrt(m_h * (1 - f_h) / (n_h - 1))
      
      # r_hi_star: Erzeuge Replikate durch Ziehen von (n_h - 1) Einheiten ohne Zurücklegen
      r_hi_star <- delta[, i, ]
      
      # rep_out = Matrix der replizierten Gewichtungen
      rep_out <- (1 - lambda + lambda * (n_h / m_h) * r_hi_star) * w_hi
      # adjustment <- sweep(r_hi_star, 1, lambda * n_h / m_h, FUN = "*") # Jede Zeile der Matrix r_hi_star wird elementweise mit dem entsprechenden Stratum-Faktor multipliziert.
      # adjustment <- sweep(adjustment, 1, 1 - lambda, FUN = "+") # 1 - lamda + lambda * n/N * r_hi_star
      # rep_out <- adjustment * w_hi
      print("n:")
      print(head(n))
      print("N:")
      print(head(N))
      
      print("delta:")
      print(head(delta))
      
      print("lambda:")
      print(lambda)
      
    }
  }
  return(rep_out)
}


next_minimum <- function(N, by) {
  N_notOne <- N != 1
  by <- by[N_notOne][which.min(N[N_notOne])]
  if (length(by) > 1) {
    by <- sample(by, 1)
  }
  return(by)
}

















# ----------------------------------------------------------------------------------------------------
# ---------------------------------------------- TEST ------------------------------------------------
# ----------------------------------------------------------------------------------------------------
#library(surveysd)
# library(data.table)
# set.seed(1234)
# eusilc <- demo.eusilc(n = 1,prettyNames = TRUE)
# 
# eusilc[,new_strata:=paste(region,hsize,sep="_")]
# eusilc[,N.housholds:=uniqueN(hid),by=new_strata]
# eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=10,strata=c("new_strata"),
#                                        cluster="hid",fpc="N.households",method = "Preston")


# Lade notwendige Bibliothek
library(data.table)
set.seed(1234)
eusilc <- demo.eusilc(n = 1,prettyNames = TRUE)

eusilc[,N.households:=uniqueN(hid),by=region]
eusilc.bootstrap.PRESTON <- rescaled.bootstrap(method="Preston", eusilc,REP=10,strata="region",
                                               cluster="hid",fpc="N.households")

eusilc.bootstrap.RAOWU <- rescaled.bootstrap(method="Rao-Wu", eusilc,REP=10,strata="region",
                                             cluster="hid",fpc="N.households")


eusilc[,new_strata:=paste(region,hsize,sep="_")]
eusilc[,N.housholds:=uniqueN(hid),by=new_strata]
eusilc.bootstrap.PRESTON <- rescaled.bootstrap(method="Preston", eusilc,REP=10,strata=c("new_strata"),
                                               cluster="hid",fpc="N.households")
eusilc.bootstrap.RAOWU <- rescaled.bootstrap(method="Rao-Wu", eusilc,REP=10,strata=c("new_strata"),
                                             cluster="hid",fpc="N.households")


# ergebnisse ansehen
head(eusilc.bootstrap.PRESTON[, .SD, .SDcols = patterns("^bootRep")])
head(eusilc.bootstrap.RAOWU[, .SD, .SDcols = patterns("^bootRep")])

eusilc.bootstrap.PRESTON[, .(total_draws = rowSums(.SD)), .SDcols = patterns("^bootRep"), by = hid]
eusilc.bootstrap.RAOWU[, .(total_draws = rowSums(.SD)), .SDcols = patterns("^bootRep"), by = hid]

var_preston <- eusilc.bootstrap.PRESTON[, lapply(.SD, var), .SDcols = patterns("^bootRep")]
var_raowu <- eusilc.bootstrap.RAOWU[, lapply(.SD, var), .SDcols = patterns("^bootRep")]
print(var_preston)
print(var_raowu)


# weighted means
weighted_means_preston <- eusilc.bootstrap.PRESTON[
  , lapply(.SD, function(w) sum(w * pWeight * eqIncome) / sum(w * pWeight)),
  .SDcols = patterns("^bootRep")
]

weighted_means_raowu <- eusilc.bootstrap.RAOWU[
  , lapply(.SD, function(w) sum(w * pWeight * eqIncome) / sum(w * pWeight)),
  .SDcols = patterns("^bootRep")
]

print(weighted_means_preston)
print(weighted_means_raowu)


library(ggplot2)
# Daten für Plot vorbereiten
preston_melt <- melt(weighted_means_preston, measure.vars = patterns("^bootRep"))
raowu_melt <- melt(weighted_means_raowu, measure.vars = patterns("^bootRep"))

ggplot() +
  geom_boxplot(data = preston_melt, aes(x = "Preston", y = value), fill = "blue", alpha = 0.5) +
  geom_boxplot(data = raowu_melt, aes(x = "Rao-Wu", y = value), fill = "red", alpha = 0.5) +
  labs(title = "Vergleich der Bootstrap-Mittelwerte", y = "Gewichtetes Einkommen (€)")

# 10%-, 50%- und 90%-Quantil pro Replikation
quantiles_preston <- eusilc.bootstrap.PRESTON[
  , lapply(.SD, function(w) 
    quantile(rep(eqIncome, round(w * pWeight)), probs = c(0.1, 0.5, 0.9))),
  .SDcols = patterns("^bootRep")
]


quantiles_raowu <- eusilc.bootstrap.RAOWU[
  , lapply(.SD, function(w) 
    quantile(rep(eqIncome, round(w * pWeight)), probs = c(0.1, 0.5, 0.9))),
  .SDcols = patterns("^bootRep")
]




###
###
# small
demo.eusilc.small <- function(n = 8, prettyNames = FALSE, dropout_fraction = 0.05) {
  
  db030 <- rb030 <- povmd60 <- eqincome <- db090 <- eqIncome <- age <- hsize <-
    . <- povertyRisk <- ecoStat <- NULL
  
  data("eusilc", package = "laeken", envir = environment())
  setDT(eusilc)
  # generate yearly data for y years
  # Adjusted dropout fraction (default: 5%)
  eusilc[, year := 2010]
  eusilc.i <- copy(eusilc)
  nsamp <- round(eusilc[, uniqueN(db030)] * dropout_fraction)
  hhincome <- eusilc[!duplicated(db030)][["eqIncome"]]
  nextIDs <- (1:nsamp) + eusilc[, max(db030)]
  if (n > 1)
    for (i in 1:(n - 1)) {
      eusilc.i[db030 %in% sample(unique(eusilc.i$db030), nsamp),
               c("db030", "eqIncome") := .(nextIDs[.GRP], sample(hhincome, 1L)),
               by = db030]
      eusilc.i[, year := year + 1]
      eusilc <- rbind(eusilc, eusilc.i)
      nextIDs <- (1:nsamp) + eusilc[, max(db030)]
    }
  
  eusilc[, rb030 := as.integer(paste0(db030, "0", 1:.N)),
         by = list(year, db030)]
  eusilc[, povmd60 := as.numeric(
    eqIncome < .6 * laeken::weightedMedian(eqIncome, w = db090)
  ), by = year]
  eusilc[, age := cut(age, c(-Inf, 16, 25, 45, 65, Inf))]
  eusilc[, hsize := cut(hsize, c(0:5, Inf))]
  
  if (prettyNames) {
    data.table::setnames(eusilc, "db030", "hid")
    data.table::setnames(eusilc, "db040", "region")
    data.table::setnames(eusilc, "rb030", "pid")
    data.table::setnames(eusilc, "rb090", "gender")
    data.table::setnames(eusilc, "pb220a", "citizenship")
    data.table::setnames(eusilc, "pl030", "ecoStat")
    eusilc[, ecoStat := factor(ecoStat, labels = c(
      "full time", "part time", "unemployed",
      "education", "retired", "disabled", "domestic"
    ))]
    data.table::setnames(eusilc, "rb050", "pWeight")
    data.table::setnames(eusilc, "povmd60", "povertyRisk")
    eusilc[, povertyRisk := as.logical(povertyRisk)]
  }
  
  return(eusilc)
}

# Beispielaufruf mit 5% Dropout
eusilc_small <- demo.eusilc.small(n = 8, prettyNames = TRUE, dropout_fraction = 0.02)
eusilc_small[,N.households:=uniqueN(hid),by=region]
eusilc_small.bootstrap.PRESTON <- rescaled.bootstrap(method="Preston", eusilc_small,REP=10,strata="region",
                                                     cluster="hid",fpc="N.households")

eusilc_small.bootstrap.RAOWU <- rescaled.bootstrap(method="Rao-Wu", eusilc_small,REP=10,strata="region",
                                                   cluster="hid",fpc="N.households")


eusilc_small[,new_strata:=paste(region,hsize,sep="_")]
eusilc_small[,N.housholds:=uniqueN(hid),by=new_strata]
eusilc_small.bootstrap.PRESTON <- rescaled.bootstrap(method="Preston", eusilc_small,REP=10,strata=c("new_strata"),
                                                     cluster="hid",fpc="N.households")
eusilc_small.bootstrap.RAOWU <- rescaled.bootstrap(method="Rao-Wu", eusilc_small,REP=10,strata=c("new_strata"),
                                                   cluster="hid",fpc="N.households")


# ergebnisse ansehen
head(eusilc_small.bootstrap.PRESTON[, .SD, .SDcols = patterns("^bootRep")])
head(eusilc_small.bootstrap.RAOWU[, .SD, .SDcols = patterns("^bootRep")])

eusilc_small.bootstrap.PRESTON[, .(total_draws = rowSums(.SD)), .SDcols = patterns("^bootRep"), by = hid]
eusilc_small.bootstrap.RAOWU[, .(total_draws = rowSums(.SD)), .SDcols = patterns("^bootRep"), by = hid]

var_preston <- eusilc_small.bootstrap.PRESTON[, lapply(.SD, var), .SDcols = patterns("^bootRep")]
var_raowu <- eusilc_small.bootstrap.RAOWU[, lapply(.SD, var), .SDcols = patterns("^bootRep")]
print(var_preston)
print(var_raowu)

# weighted means
weighted_means_preston_small <- eusilc_small.bootstrap.PRESTON[
  , lapply(.SD, function(w) sum(w * pWeight * eqIncome) / sum(w * pWeight)),
  .SDcols = patterns("^bootRep")
]

weighted_means_raowu_small <- eusilc_small.bootstrap.RAOWU[
  , lapply(.SD, function(w) sum(w * pWeight * eqIncome) / sum(w * pWeight)),
  .SDcols = patterns("^bootRep")
]

print(weighted_means_preston_small)
print(weighted_means_raowu_small)


library(ggplot2)
# Daten für Plot vorbereiten
preston_melt_small <- melt(weighted_means_preston_small, measure.vars = patterns("^bootRep"))
raowu_melt_small <- melt(weighted_means_raowu_small, measure.vars = patterns("^bootRep"))

ggplot() +
  geom_boxplot(data = preston_melt_small, aes(x = "Preston", y = value), fill = "blue", alpha = 0.5) +
  geom_boxplot(data = raowu_melt_small, aes(x = "Rao-Wu", y = value), fill = "red", alpha = 0.5) +
  labs(title = "Vergleich der Bootstrap-Mittelwerte", y = "Gewichtetes Einkommen (€)")

# 10%-, 50%- und 90%-Quantil pro Replikation
quantiles_preston <- eusilc_small.bootstrap.PRESTON[
  , lapply(.SD, function(w) 
    quantile(rep(eqIncome, round(w * pWeight)), probs = c(0.1, 0.5, 0.9))),
  .SDcols = patterns("^bootRep")
]


quantiles_raowu <- eusilc_small.bootstrap.RAOWU[
  , lapply(.SD, function(w) 
    quantile(rep(eqIncome, round(w * pWeight)), probs = c(0.1, 0.5, 0.9))),
  .SDcols = patterns("^bootRep")
]

# ----------------------------------------------------------------------------------------------------




#############################################################################################################################################################















# 
# 
# # TEST -----------------------------------------------------------------------------------------------------------------
# set.seed(123) # Für reproduzierbare Ergebnisse
# 
# # Anzahl der Beobachtungen
# n_rows <- 10
# n_cols <- 3
# n_replicates <- 5
# 
# # Erstellen der Matrizen und Vektoren
# n <- matrix(sample(5:20, n_rows * n_cols, replace = TRUE), nrow = n_rows, ncol = n_cols)  # Stichprobenanzahl
# N <- matrix(n * sample(20:50, n_cols, replace = TRUE), nrow = n_rows, ncol = n_cols)     # Grundgesamtheiten
# n_draw <- matrix(sample(1:5, n_rows * n_cols, replace = TRUE), nrow = n_rows, ncol = n_cols)  # Anzahl der gezogenen Einheiten
# # Generieren von delta für die Preston-Methode (Binär: 0 oder 1)
# delta_preston <- array(sample(0:1, n_rows * n_cols * n_replicates, replace = TRUE), dim = c(n_rows, n_cols, n_replicates))  # Resampling-Indikatoren
# 
# # Generieren von delta für die Rao-Wu-Methode (Ganzzahlig: Häufigkeit der Ziehungen)
# delta_rao_wu <- array(0, dim = c(n_rows, n_cols, n_replicates))  # Delta initialisieren
# set.seed(123)  # Für Reproduzierbarkeit
# for (r in 1:n_replicates) {
#   for (i in 1:n_cols) {
#     for (j in 1:n_rows) {
#       # Häufigkeit der Ziehungen mit Zurücklegen
#       delta_rao_wu[j, i, r] <- sum(sample(1:n[j, i], n_draw[j, i], replace = TRUE) == 1)
#     }
#   }
# }
# 
# 
# # Zeilen repräsentieren die verschiedenen Einheiten oder Strata in der Stichprobe
# # Spalten repräsentieren die verschiedenen Stufen des Stichprobenverfahrens. Ein mehrstufiges Stichprobenverfahren kann folgende Schritte beinhalten: 
# # Stufe 1 (Primärstufe): Auswahl von Primäreinheiten (z. B. Dörfer, Schulen).
# # Stufe 2 (Sekundärstufe): Auswahl von Sekundäreinheiten innerhalb der Primäreinheiten (z. B. Haushalte in einem Dorf, Klassen in einer Schule).
# # Stufe 3 (Tertiärstufe): Auswahl von Tertiäreinheiten innerhalb der Sekundäreinheiten (z. B. Personen in einem Haushalt, Schüler in einer Klasse).
# # n[i, j] = Stichprobenanzahl =  Beispiel: n[2, 1] könnte die Anzahl der Dörfer sein, die in Stratum 2 auf der Primärstufe gezogen wurden.
# # N[i, j] = Grundgesamtheiten = Beispiel: N[3, 2] könnte die Anzahl der verfügbaren Haushalte in Stratum 3 auf der Sekundärstufe sein.
# # delta -> Dimensionen (Zeilen, Spalten, Replikate) = Beispiel: delta[1, 2, 3] gibt an, ob eine Einheit aus Stratum 1 auf Stufe 2 im dritten Replikat ausgewählt wurde
# 
# 
# # Testen der Funktion mit der "Preston"-Methode
# result_preston <- calc.replicate(n, N, n_draw, delta_preston, method = "Preston")
# print("Ergebnis mit der Preston-Methode:")
# print(result_preston)
# 
# # Testen der Funktion mit der "Rao-Wu"-Methode
# result_rao_wu <- calc.replicate(n, N, n_draw, delta_rao_wu, method = "Rao-Wu")
# print("Ergebnis mit der Rao-Wu-Methode:")
# print(result_rao_wu)
# # TEST -----------------------------------------------------------------------------------------------------------------





################################################################################################################
##################################### HELPER ###################################################################
################################################################################################################
rowProds <- function(x, na.rm = FALSE) {
  n <- nrow(x)
  y <- double(length = n)
  for (ii in seq_len(n)) {
    y[ii] <- prod(x[ii, , drop = TRUE], na.rm = na.rm)
  }
  y
}


dt.eval <- function(..., env = parent.frame()) {
  
  expressions <- paste0(...)
  if (length(expressions) > 1) {
    return(lapply(
      expressions,
      function(z) {
        eval(parse(text = z), envir = env)
      }
    ))
  } else {
    return(eval(parse(text = paste0(...)), envir = env))
  }
}

dt.eval2 <- function(...) {
  return(eval(parse(text = paste0(...)), envir = parent.frame()))
}

getEllipsis <- function(element, default, ell) {
  ifelse(is.null(ell[[element]]), default, ell[[element]])
}
getEllipsis2 <- function(element, default, ell) {
  
  if (element %in% names(ell)) {
    return(ell[[element]])
  }else{
    return(default)
  }
}

# helpfunction to create contingency tables
makeCalibTable <- function(dat, weights, period, vars) {
  # make contingency table
  formTab <- paste(weights, "~", paste(c(period, vars), collapse = "+"))
  varsTab <- xtabs(formTab, data = dat)
  return(list(varsTab))
}

paste_ <- function(a, b) {
  paste(a, b, sep = ".")
}
paste_c <- function(a, b) {
  paste(a, b, sep = ",")
}
paste_addarg <- function(a, b) {
  
  a <- tstrsplit(a, ",")
  
  return(paste(a[[1]], b, paste(a[2:length(a)], collapse = ","), sep = ","))
}

# helpfunctions for point estimates
weightedRatioNat <- function(x, w, N) {
  weightedRatio(x, w) / N * 100
}
weightedRatioR <- function(x, w) {
  sum(w[x == 1], na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE) * 100
}
weightedSumR <- function(x, w) {
  sum(as.numeric(x) * w, na.rm = TRUE)
}

povmd <- function(x, w) {
  md <- laeken::weightedMedian(x, w) * 0.6
  pmd60 <- x < md
  return(as.integer(pmd60))
}

# helpfunction for quantile calcultion with missings
quantileNA <- function(x, probs, p.names, np = length(probs)) {
  
  if (any(is.na(x))) {
    out <- rep(NA_real_, np)
  } else {
    out <- quantile(x, probs = probs)
  }
  names(out) <- p.names
  return(out)
}

randomInsert <- function(x, y, n = 20) {
  if (length(x) < 20 | length(y) < 20) {
    stop("n must be smaller than length(x) and length(y)")
  }
  
  x.indices <- sample(length(x), n)
  y.values <- sample(y, n)
  x[x.indices] <- y.values
  return(x)
}

generateRandomName <- function(nchar = 20, existingNames) {
  
  newName <- paste(sample(c(letters, LETTERS), nchar), collapse = "")
  while (newName %in% existingNames) {
    newName <- paste(sample(c(letters, LETTERS), nchar), collapse = "")
  }
  
  
  return(newName)
}

# input checking function
check.input <- function(input, input.name, input.length=NULL, 
                        input.type = NULL, decimal.possible = NULL, 
                        c.names = NULL, dat = NULL, dat.column.type = NULL){
  
  if(!is.null(input.length)){
    if(length(input) != 1){
      stop(paste(input.name,"must have length",input.length))
    }
  }
  
  if(!is.null(input.type)){
    if(!is(input,input.type)){
      stop(paste(input.name,"must be of type",input.type))
    }
  }
  
  if(!is.null(decimal.possible)){
    if(input %%1 !=0 & decimal.possible==FALSE){
      stop(paste(input.name,"cannot have a decimal part"))
    }
  }
  
  if(!is.null(c.names)){
    if(!input %in% c.names){
      stop(paste(input," is not a column in dat"))
    }
  }
  
  if(!is.null(dat) & !is.null(dat.column.type)){
    if(!is(dat[[input]],dat.column.type)){
      stop(paste(input.name,"must be a",input.type,"column in dat"))
    }
  }
  
  return(NULL)
}
################################################################################################################
##################################### HELPER ###################################################################
################################################################################################################
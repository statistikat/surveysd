###################################################
# test rescaled bootstrap
#
#
#

library(matrixStats)
library(data.table)

REP <- 10
p <- 2


draw.without.replacement <- function(n,N){
  n_draw <- trunc(n/2)
  # n_draw <- trunc(n/(2-n/N))
  delta <- rep(c(1,0),c(n_draw,n-n_draw))
  delta <- sample(delta)
  return(delta)
}

calc.replicate <- function(n,N,delta){
  p <- ncol(n)
  ndraw <- trunc(n/2)
  # ndraw <- trunc(n/(2-n/N))
  dimdelta <- dim(delta)
  for(i in 1:p){
    if(i==1){
      lambda <- sqrt(ndraw[,1]*(1-n[,1]/N[,1])/(n[,1]-ndraw[,1]))
      rep_out <- 1-lambda+lambda*n[,i]/ndraw[,i]*delta[,i,]
    }else if(i==2){
      lambda <- (1-n[,i]/N[,i])/(n[,i]-ndraw[,i])
      lambda <- sqrt((n[,i-1]/N[,i-1])*ndraw[,i]*lambda)
      rep_out <- rep_out + lambda*(sqrt(n[,i-1]/ndraw[,i-1])*delta[,i-1,]) * (n[,i]/ndraw[,i]*delta[,i,]-1)
    }else{
      lambda <- (1-n[,i]/N[,i])/(n[,i]-ndraw[,i])
      lambda <- sqrt(rowProds(n[,1:(i-1)]/N[,1:(i-1)])*ndraw[,i]*lambda)
      prod_val <- matrix(0,ncol=dimdelta[3],nrow=dimdelta[1])
      for(r in 1:dimdelta[3]){
        prod_val[,r] <- rowProds(sqrt(n[,1:(i-1)]/ndraw[,1:(i-1)])*delta[,1:(i-1),r])
      }
      # rep_out <- rep_out + lambda*rowProds(sqrt(n[,1:(i-1)]/ndraw[,1:(i-1)])*delta[,1:(i-1),]) * (n[,i]/ndraw[,i]*delta[,i,]-1)
      rep_out <- rep_out + lambda*prod_val * (n[,i]/ndraw[,i]*delta[,i,]-1)
    }
  }
  return(rep_out)
}


# grundgesamtheit
s <- 10
N1 <- sample(3:10,s,replace=TRUE)
n1 <- N1 - sample(0:1,s,replace=TRUE)
id <- rep(1:s,times=n1)
N1 <- rep(N1,times=n1)
n1 <- rep(n1,times=n1)

dat <- data.table(id,N1,n1)
dat[,id1:=1:.N,by=id]

N2 <- sample(10:30,length(N1),replace=TRUE)
n2 <- N2 - sample(2:8,length(N1),replace=TRUE)

dat <- dat[,.(id2=1:n2[.GRP],N2=N2[.GRP],n2=n2[.GRP],N1=N1[1],n1=n1[1]),by=list(id,id1)]

n.row <- nrow(dat)
delta <- array(0,dim=c(n.row,p,REP))
n <- as.matrix(dat[,.(n1,n2)])
N <- as.matrix(dat[,.(N1,N2)])

dati <- unique(dat,by=c("id","id1"))
c.names <- paste0("V",1:REP)
dati[,c(c.names):=as.data.table(replicate(REP,draw.without.replacement(n1[1],N1[1])),simplify=FALSE),by=c("id")]
deltai <- dati[,mget(c("id","id1",c.names))]

dat <- merge(dat,deltai,by=c("id","id1"))
delta[,1,] <- as.matrix(dat[,mget(c.names)])
dat[,c(c.names):=NULL]

dati <- unique(dat,by=c("id","id1","id2"))
dati[,c(c.names):=as.data.table(replicate(REP,draw.without.replacement(n2[1],N2[1])),simplify=FALSE),by=c("id","id1")]
deltai <- dati[,mget(c("id","id1","id2",c.names))]

dat <- merge(dat,deltai,by=c("id","id1","id2"))
delta[,2,] <- as.matrix(dat[,mget(c.names)])
dat[,c(c.names):=NULL]

b.rep <- calc.replicate(n,N,delta)
colSums(b.rep)
any(b.rep<0)

hist(b.rep[,1])
           

#####################################
# teste an laeken
#

library(laeken)
data("eusilc")
eusilc <- data.table(eusilc)

eusilc[!duplicated(db030),N.bdl := .N+sample(1:50,1),by=db040]

eusilc[,N.db030 := .N+sample(0:4,1),by=db030]







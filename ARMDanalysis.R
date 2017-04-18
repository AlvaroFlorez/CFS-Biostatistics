library(Surrogate)
library(pan)
library(tictoc)
#### Add NA to the dataset
AddNA <- function(Data){
### function to add NA values for the true and surrogate endpoints in order to perform the MI procedure
  ni <- table(Data[,1],Data[,3])/2
  Cluster <- unique(Data[,1])
  nmax <- max(ni)
  n.add <- nmax-ni
  idmax <- max(Data[,2])
  ID.add <- NULL
  Endpoint.add <- NULL
  Cluster.add <- rep(Cluster,colSums(t(n.add))*2)
  Treat.add <- NULL
  for(i in 1:length(Cluster)){
#    if(Cluster[i]!=43){
    nn <- sum(n.add[i,])
    if(nn > 0){
    Treat.add <- c(Treat.add,rep(c(rep(-1,n.add[i,1]),rep(1,n.add[i,2])),2))
    ID.add <- c(ID.add,rep((idmax+1):(idmax+sum(n.add[i,])),2))
    Endpoint.add <- c(Endpoint.add,c(rep(-1,sum(n.add[i,])),c(rep(1,sum(n.add[i,])))))
    idmax <- max(ID.add)#}
    }
  }
  Outcome.add <- rep(NA,sum(n.add)*2)
  Data.add <- data.frame(Cluster.add,ID.add,Treat.add,Endpoint.add,Outcome.add)
  names(Data.add) <- names(Data)
  Data <- rbind(Data,Data.add)
  Data <- Data[order(Data[,1],Data$Endpoint,Data$Treat,Data$ID),]
  return(Data)
}
#### MI procedure using PAN package
MIdata <- function(Data,K,Seed,Burn,Iter){
#### Data: ARMD data
#### K: number of imputations
#### Seed: seed value
#### Burn: Burn-in period
#### Iter: number of iterations to get imputed values
  ID <- unique(Data$ID)
  Y <- cbind(Data$Outcome[Data$Endpoint==-1],Data$Outcome[Data$Endpoint==1])
  Center <- Data$Center[Data$Endpoint==1]
  n <- dim(Y)[1]/length(unique(Center))
  Endpoint <- c(rep(-1,dim(Y)[1]),rep(1,dim(Y)[1]))
  X <- cbind(1, Data$Treat[Data$Endpoint==1])
  Z <- X
  xcol=zcol=1:2
  prior <- list(a=2,Binv=diag(1,4),c=2*4,Dinv=diag(1,4))
  CompTime <- rep(0,K)
  tic()
  Results <-  pan(Y,Center,X,xcol,zcol,prior,seed=Seed,iter=Burn+Iter)
  exectime<- toc()
  CompTime[1] <- exectime$toc - exectime$tic
  
  Imputation <- cbind(1,Center,ID,Endpoint,X[,2],as.vector(Results$y))
  for(i in 2:K){
    Seed <- Seed + 1
    tic()
    Results <-  pan(Y,Center,X,xcol,zcol,prior,seed=Seed,iter=Iter,start=Results$last)
    exectime<- toc()
    CompTime[k] <- exectime$toc - exectime$tic
    Imputation <- rbind(Imputation,cbind(i,Center,ID,Endpoint,X[,2],as.vector(Results$y)))
  }
  
  Imputation <- as.data.frame(Imputation)
  names(Imputation) <- c("Imp","Center","ID","Endpoint","Treat","Outcome")
  Imputation <- Imputation[order(Imputation$Imp,Imputation$Center,Imputation$Endpoint,Imputation$Treat),]
  write(CompTime,"ARMD.MIpanCompTime.txt", sep="\t")
  write.table(Imputation,"ARMD.MIpan.txt", sep="\t",row.names=FALSE) 
  return(Imputation)
}  

R2TrialFun <- function(D){
  A <- matrix(c(D[1,4], D[2,4]),2,1)
  B <- matrix(c(D[1,1], D[1,2], D[1,2], D[2,2]),2,2)
  C <- D[4,4]
  R2.trial <- t(A)%*%solve(B)%*%A/C
  return(R2.trial)
}
R2IndFun <- function(Sigma){
  R2.ind <- Sigma[1,2]^2/(Sigma[1,1]*Sigma[2,2])
  return(R2.ind)
}
### estimation using closed-form solution for Data
ClosedFS <-function(correction,Data,REML=FALSE){
  ## correction: type of correction used for a possible non-pd D solution (=1 correction 1 and =2 correction 2)
  ## Data: Dataset with outcome in column 7, treatment effect in column ? .....
  ## REML: Restricted maximum likelihood estimator (TRUE) or Maximum likelihood estimator (FALSE)
  N <- length(unique(Data[,1]))
  n <-unique(table(Data[,1]))/2
  X <- cbind(1,Data[1:n,3])
  Z <- t(bdiag(X,X))
  Y <- matrix(Data[,5],N,2*n,byrow=TRUE)
  ZZ <- tcrossprod(Z)
  #tA.A <- 1/N
  ### estimation of fixed effects
  B <- solve(ZZ,tcrossprod(Z,Y))
  BetaH <- apply(B,1,mean)
  
  
  Hx = X%*%solve(crossprod(X),t(X))
  E <- Y- crossprod(matrix(1,1,N),BetaH%*%Z)
  ES = E[,c(1:n)]
  ET = E[,c((n+1):(2*n))]
  
  Sum.R11 <- sum(diag(crossprod(diag(1,n) - Hx,crossprod(ES,ES))))
  Sum.R22 <- sum(diag(crossprod(diag(1,n) - Hx,crossprod(ET,ET))))
  Sum.R12 <- sum(diag(crossprod(diag(1,n) - Hx,crossprod(ES,ET))))
  
  SigmaH <- matrix(c(Sum.R11,Sum.R12,Sum.R12,Sum.R22),2,2)/(N*(n-2))
  #### estimation of the variance matrix of the random effects
  if(REML==FALSE){
    GH <- (1/N)*tcrossprod(tcrossprod(B,diag(1,N)-matrix(1/N,N,N)),B)  
  }
  else{
    GH <- (1/(N-1))*tcrossprod(tcrossprod(B,diag(1,N)-matrix(1/N,N,N)),B)  
    
  }  
  DH <- GH - solve(ZZ,kronecker(SigmaH,diag(1,2)))
  
  
  pdm.DH <- sum(eigen(DH)$values>-1e-7) == 4
  
  if(pdm.DH == 0){
    if(correction==1){
      Sigma2H <- kronecker(SigmaH,diag(1,2))
      G <- chol(ZZ)
      Ginv <- solve(G)
      Delta <- Sigma2H + tcrossprod(tcrossprod(G,DH),G)
      EigenV <- eigen(solve(Sigma2H,Delta))
      E <- diag(EigenV$values)
      diag(E)[diag(E)<=1] <- 0.0001
      L <- EigenV$vector
      I <- diag(1,4)
      I[E<=1] = 0
      
      DH <- Sigma2H%*%Ginv%*%L%*%(E-I)%*%solve(L)%*%t(Ginv)
      #pdm.DH <- sum(eigen(DH)$values>-1e-7) == 4
    }
    if(correction==2){
      EigenD <- eigen(DH)
      E <- diag(EigenD$values)
      diag(E)[diag(E)<=0] <- 0.0001
      L <- EigenD$vector
      DH <- L%*%E%*%t(L)
    }
  }  
  R2Trial <- R2TrialFun(DH)
  R2Ind <- R2IndFun(SigmaH)
  Estimates <- c(as.vector(BetaH),DH[upper.tri(DH,diag=TRUE)],SigmaH[upper.tri(SigmaH,diag=TRUE)],R2Trial,R2Ind,pdm.DH)
  return(Estimates)
}

#### data manipulation for ARMD dataset. Wide to long format and adding NA values for MI
data(ARMD)
ni <- table(ARMD$Center)
ARMD$n <- rep(ni,ni)
## remove the case with ni < 5
ARMD2 <- ARMD[ARMD$n > 5,]
### wide to long
ID <- rep(ARMD2$Id,2)
Center <- rep(ARMD2$Center,2)
Treat <- rep(ARMD2$Treat,2)
Outcome <- c(ARMD2$Diff24,ARMD2$Diff52)
Endpoint <- c(rep(-1,dim(ARMD2)[1]),rep(1,dim(ARMD2)[1]))
ARMDlong <- data.frame(Center,ID,Treat,Endpoint,Outcome,row.names = NULL)
ARMDlong <- ARMDlong[order(ARMDlong$Center,ARMDlong$Endpoint,ARMDlong$Treat,ARMDlong$ID),]
ARMDlong2 <- AddNA(ARMDlong)

write.table(ARMDlong2,"ARMD2.withNA.txt", sep="\t",row.names=FALSE) 


#### MI procedure using PAN package
K <- 1000
ARMD2 <- read.table("ARMD2.withNA.txt",header=TRUE)
MIARMD <- MIdata(ARMD2,K,1,2000,5000)
write.table(MIARMD,"ARMD.Mipan.txt", sep="\t",row.names=FALSE) 


#### CFS for ARMD-MI data (using PAN package)
MIARMD <- read.table("ARMD.MIpan.txt",header=TRUE)
MIARMD <- MIARMD[,c(1,2,3,5,4,6)]
Estimates <- matrix(0,K,20)
for(k in 1:K){
  Data <- MIARMD[MIARMD$Imp==k,-1]
  Estimates[k,] <- CFS(Data)
}
write.table(Estimates,"ARMDdata/MI2Estimates.CFS.txt", sep="\t",row.names=FALSE) 



### MI using trial by trial MI (performed in SAS)
MI.SAS.ARMD <- read.table("ARMDdata/ARMD.MI.SAS.txt",header=TRUE)
  K <- 1000
  EstimatesSAS <- matrix(0,K,20)
  CompTime <- rep(0,1000)

CompTime2 <- c()
tic()
for(k in 1:K){
    Data <- MI.SAS.ARMD[MI.SAS.ARMD$X_Imputation==k,-1]
    EstimatesSAS[k,] <- ClosedFS(2,Data,REML=TRUE)
  }
exectime<- toc()
CompTime2 <- exectime$toc - exectime$tic

  EstimatesSAS <- cbind(EstimatesSAS,CompTime)
  write.table(EstimatesSAS,"ARMDdata/MI1Estimates.CFS.txt", sep="\t",row.names=FALSE) 

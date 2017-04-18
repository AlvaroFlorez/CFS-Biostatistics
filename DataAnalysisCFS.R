library(Matrix)
library(tictoc)
library(data.table)
### Closed-form solution function (N cluster of the same size)
ClosedFS <-function(correction,Data,REML=FALSE){
  ## correction: type of correction used for a possible non-pd D solution (=1 correction 1 and =2 correction 2)
  ## Data: Dataset with outcome in column 7, treatment effect in column ? .....
  ## REML: Restricted maximum likelihood estimator (TRUE) or Maximum likelihood estimator (FALSE)
  N <- length(unique(Data$Trial.id))
  n <-unique(table(Data$Trial.id))/2
  X <- cbind(1,Data$Treat[1:n])
  Z <- t(bdiag(X,X))
  Y <- matrix(Data$Outcome,N,2*n,byrow=TRUE)
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
### function to compute R^2_trial
R2TrialFun <- function(D){
  A <- matrix(c(D[1,4], D[2,4]),2,1)
  B <- matrix(c(D[1,1], D[1,2], D[1,2], D[2,2]),2,2)
  C <- D[4,4]
  R2.trial <- t(A)%*%solve(B)%*%A/C
  return(R2.trial)
}
### function to compute R^2_trial
R2IndFun <- function(Sigma){
  R2.ind <- Sigma[1,2]^2/(Sigma[1,1]*Sigma[2,2])
  return(R2.ind)
}
### split sample method
SS.Est <- function(data,correction,REML){
## data: dataset
## correction: type of correction for non-PD solution
## REML: estimation method( = FALSE ML and = TRUE REML)
  N.subsamples <- length(unique(data[,2]))
  n.k <- c(by(data[,3],data[,2],function(x){unique(table(x))}))/2
  c.k <- c(by(data[,3],data[,2],function(x){length(table(x))}))
  
  equal.w <- rep(1/N.subsamples,N.subsamples)
  propC.w <- c.k/sum(c.k)
  propCN.w <- (c.k*n.k)/sum(c.k*n.k)
  Estimates <- by(data,data[,2],ClosedFS,correction=correction,REML=REML)
  Est.equal <- Reduce('+',mapply("*",Estimates,equal.w,SIMPLIFY = FALSE))
  Est.propC <- Reduce('+',mapply("*",Estimates,propC.w,SIMPLIFY = FALSE))
  Est.propCN <- Reduce('+',mapply("*",Estimates,propCN.w,SIMPLIFY = FALSE))
  
  D.equal <- matrix(0,4,4)
  D.equal[upper.tri(D.equal,diag=TRUE)] <-   Est.equal[5:14]
  R2.trial.equal <- R2TrialFun(D.equal)
  
  D.propC <- matrix(0,4,4)
  D.propC[upper.tri(D.propC,diag=TRUE)] <-   Est.propC[5:14]
  R2.trial.propC <- R2TrialFun(D.propC)
  
  D.propCN <- matrix(0,4,4)
  D.propCN[upper.tri(D.propCN,diag=TRUE)] <-   Est.propCN[5:14]
  R2.trial.propCN <- R2TrialFun(D.propCN)

  Sigma.equal <- matrix(0,2,2)
  Sigma.equal[upper.tri(Sigma.equal,diag=TRUE)] <-   Est.equal[15:17]
  R2.ind.equal <- R2IndFun(Sigma.equal)
  
  Sigma.propC <- matrix(0,2,2)
  Sigma.propC[upper.tri(Sigma.propC,diag=TRUE)] <-   Est.propC[15:17]
  R2.ind.propC <- R2IndFun(Sigma.propC)
  
  Sigma.propCN <- matrix(0,2,2)
  Sigma.propCN[upper.tri(Sigma.propCN,diag=TRUE)] <-   Est.propCN[15:17]
  R2.ind.propCN <- R2IndFun(Sigma.propCN)
  
  
  Est.equal <- c(Est.equal,R2.trial.equal,R2.ind.equal)
  Est.propC <- c(Est.propC,R2.trial.propC,R2.ind.propC)
  Est.propCN <- c(Est.propCN,R2.trial.propCN,R2.ind.propCN)
  return(c(Est.equal,Est.propC,Est.propCN))
}
### split sample method for all simulated datasets
EstimateFun <- function(N,D,alt,correction,REML,Return=FALSE){
  Case <- ifelse(alt,paste("DataCaseN",N,"D",D,'alt.txt',sep=''),paste("DataCaseN",N,"D",D,'.txt',sep=''))
  Data <- read.table(paste('datasets/',Case,sep=''),header=TRUE) 
  
  M <- max(Data$Sim.id)
  Estimates <- matrix(0,M,66)
  for(m in 1:M){
    data.m <- Data[Data$Sim.id==m,]
    Estimates[m,] = SS.Est(data.m,correction=correction,REML=REML)
  }
  Estimates.equal = Estimates[,1:22]
  Estimates.propC = Estimates[,23:44]
  Estimates.propCN = Estimates[,45:66]

  colnames(Estimates.equal) <- colnames(Estimates.propC) <- colnames(Estimates.propCN) <- c("BS0","BS1",
                                                                                            "BT0","BT1","dSS","dSa","daa",
                                                                                            "dST","daT","dTT","dSb","dab","dTb",
                                                                                            "dbb","SigmaSS","SigmaST","SigmaTT","R2Trial",
                                                                                            "R2Ind","Valid","R2Trial2","R2Ind2")
  Case <- ifelse(alt,paste("N",N,"D",D,'alt',sep=''),paste("N",N,"D",D,sep=''))
  if(REML==TRUE){
    write.table(Estimates.equal, paste("estimates/",Case,".equal",".REML.txt",sep=""), sep="\t",row.names=FALSE)
    write.table(Estimates.propC, paste("estimates/",Case,".propC",".REML.txt",sep=""), sep="\t",row.names=FALSE)
    write.table(Estimates.propCN, paste("estimates/",Case,".propCN",".REML.txt",sep=""), sep="\t",row.names=FALSE)
  }
  if(REML==FALSE){
    write.table(Estimates.equal, paste("estimates/",Case,".equal",".txt",sep=""), sep="\t",row.names=FALSE)
    write.table(Estimates.propC, paste("estimates/",Case,".propC",".txt",sep=""), sep="\t",row.names=FALSE)
    write.table(Estimates.propCN, paste("estimates/",Case,".propCN",".txt",sep=""), sep="\t",row.names=FALSE)
  }
  if(Return==TRUE){
    Results <- list(equal = Estimates.equal, propC = Estimates.propC,
                    propCN = Estimates.propCN)
    return(Results)
  }
}

EstimateFun(10,1,TRUE,2,TRUE,TRUE)
EstimateFun(10,2,TRUE,2,TRUE,FALSE)
EstimateFun(10,1,FALSE,2,TRUE,FALSE)
EstimateFun(10,2,FALSE,2,TRUE,FALSE)

EstimateFun(20,1,TRUE,2,TRUE,FALSE)
EstimateFun(20,2,TRUE,2,TRUE,FALSE)
EstimateFun(20,1,FALSE,2,TRUE,FALSE)
EstimateFun(20,2,FALSE,2,TRUE,FALSE)

EstimateFun(50,1,TRUE,2,TRUE,FALSE)
EstimateFun(50,2,TRUE,2,TRUE,FALSE)
EstimateFun(50,1,FALSE,2,TRUE,FALSE)
EstimateFun(50,2,FALSE,2,TRUE,FALSE)

#########
EstimateFun(10,1,TRUE,2,FALSE,TRUE)
EstimateFun(10,2,TRUE,2,FALSE,FALSE)
EstimateFun(10,1,FALSE,2,FALSE,FALSE)
EstimateFun(10,2,FALSE,2,FALSE,FALSE)

EstimateFun(20,1,TRUE,2,FALSE,FALSE)
EstimateFun(20,2,TRUE,2,FALSE,FALSE)
EstimateFun(20,1,FALSE,2,FALSE,FALSE)
EstimateFun(20,2,FALSE,2,FALSE,FALSE)

EstimateFun(50,1,TRUE,2,FALSE,FALSE)
EstimateFun(50,2,TRUE,2,FALSE,FALSE)
EstimateFun(50,1,FALSE,2,FALSE,FALSE)
EstimateFun(50,2,FALSE,2,FALSE,FALSE)






######### ADDITIONAL SIMULATIONS
EstimateFunAdd <- function(N,D,alt,R2,correction,REML,Return=FALSE){
  if(R2==.25){
    Case <- ifelse(alt,paste("N",N,"D",D,'alt.txt',sep=''),paste("N",N,"D",D,'.txt',sep=''))
  }else{
    Case <- ifelse(alt,paste("N",N,"D",D,'alt75.txt',sep=''),paste("N",N,"D",D,'75.txt',sep=''))
  }
  Data <- read.table(paste('Additional Simulation/datasets/',Case,sep=''),header=TRUE) 
  
  M <- max(Data$Sim.id)
  Estimates <- matrix(0,M,66)
  for(m in 1:M){
    data.m <- Data[Data$Sim.id==m,]
    Estimates[m,] = SS.Est(data.m,correction=correction,REML=REML)
  }
  Estimates.equal = Estimates[,1:22]
  Estimates.propC = Estimates[,23:44]
  Estimates.propCN = Estimates[,45:66]
  
  colnames(Estimates.equal) <- colnames(Estimates.propC) <- colnames(Estimates.propCN) <- c("BS0","BS1",
                                                                                            "BT0","BT1","dSS","dSa","daa",
                                                                                            "dST","daT","dTT","dSb","dab","dTb",
                                                                                            "dbb","SigmaSS","SigmaST","SigmaTT","R2Trial",
                                                                                            "R2Ind","Valid","R2Trial2","R2Ind2")
  Case <- ifelse(alt,paste("N",N,"D",D,'alt',sep=''),paste("N",N,"D",D,sep=''))
  if(REML==TRUE){
    write.table(Estimates.equal, paste("Additional Simulation/Estimates/",Case,".equal",".REML.txt",sep=""), sep="\t",row.names=FALSE)
    write.table(Estimates.propC, paste("Additional Simulation/Estimates/",Case,".propC",".REML.txt",sep=""), sep="\t",row.names=FALSE)
    write.table(Estimates.propCN, paste("Additional Simulation/Estimates/",Case,".propCN",".REML.txt",sep=""), sep="\t",row.names=FALSE)
  }
  if(REML==FALSE){
    write.table(Estimates.equal, paste("Additional Simulation/Estimates/",Case,".equal",".txt",sep=""), sep="\t",row.names=FALSE)
    write.table(Estimates.propC, paste("Additional Simulation/Estimates/",Case,".propC",".txt",sep=""), sep="\t",row.names=FALSE)
    write.table(Estimates.propCN, paste("Additional Simulation/Estimates/",Case,".propCN",".txt",sep=""), sep="\t",row.names=FALSE)
  }
  if(Return==TRUE){
    Results <- list(equal = Estimates.equal, propC = Estimates.propC,
                    propCN = Estimates.propCN)
    return(Results)
  }
}

EstimateFunAdd(10,1,TRUE,0.25,2,TRUE,TRUE)
EstimateFunAdd(10,2,TRUE,0.25,2,TRUE,FALSE)
EstimateFunAdd(10,1,FALSE,0.25,2,TRUE,FALSE)
EstimateFunAdd(10,2,FALSE,0.25,2,TRUE,FALSE)

EstimateFunAdd(20,1,TRUE,0.25,2,TRUE,TRUE)
EstimateFunAdd(20,2,TRUE,0.25,2,TRUE,FALSE)
EstimateFunAdd(20,1,FALSE,0.25,2,TRUE,FALSE)
EstimateFunAdd(20,2,FALSE,0.25,2,TRUE,FALSE)

EstimateFunAdd(50,1,TRUE,0.25,2,TRUE,TRUE)
EstimateFunAdd(50,2,TRUE,0.25,2,TRUE,FALSE)
EstimateFunAdd(50,1,FALSE,0.25,2,TRUE,FALSE)
EstimateFunAdd(50,2,FALSE,0.25,2,TRUE,FALSE)

#########
EstimateFunAdd(10,1,TRUE,0.25,2,FALSE,TRUE)
EstimateFunAdd(10,2,TRUE,0.25,2,FALSE,FALSE)
EstimateFunAdd(10,1,FALSE,0.25,2,FALSE,FALSE)
EstimateFunAdd(10,2,FALSE,0.25,2,FALSE,FALSE)

EstimateFunAdd(20,1,TRUE,0.25,2,FALSE,TRUE)
EstimateFunAdd(20,2,TRUE,0.25,2,FALSE,FALSE)
EstimateFunAdd(20,1,FALSE,0.25,2,FALSE,FALSE)
EstimateFunAdd(20,2,FALSE,0.25,2,FALSE,FALSE)

EstimateFunAdd(50,1,TRUE,0.25,2,FALSE,TRUE)
EstimateFunAdd(50,2,TRUE,0.25,2,FALSE,FALSE)
EstimateFunAdd(50,1,FALSE,0.25,2,FALSE,FALSE)
EstimateFunAdd(50,2,FALSE,0.25,2,FALSE,FALSE)

EstimateFunAdd(10,1,TRUE,0.75,2,TRUE,TRUE)
EstimateFunAdd(10,2,TRUE,0.75,2,TRUE,FALSE)
EstimateFunAdd(10,1,FALSE,0.75,2,TRUE,FALSE)
EstimateFunAdd(10,2,FALSE,0.75,2,TRUE,FALSE)

EstimateFunAdd(20,1,TRUE,0.75,2,TRUE,TRUE)
EstimateFunAdd(20,2,TRUE,0.75,2,TRUE,FALSE)
EstimateFunAdd(20,1,FALSE,0.75,2,TRUE,FALSE)
EstimateFunAdd(20,2,FALSE,0.75,2,TRUE,FALSE)

EstimateFunAdd(50,1,TRUE,0.75,2,TRUE,TRUE)
EstimateFunAdd(50,2,TRUE,0.75,2,TRUE,FALSE)
EstimateFunAdd(50,1,FALSE,0.75,2,TRUE,FALSE)
EstimateFunAdd(50,2,FALSE,0.75,2,TRUE,FALSE)

#########
EstimateFunAdd(10,1,TRUE,0.75,2,FALSE,TRUE)
EstimateFunAdd(10,2,TRUE,0.75,2,FALSE,FALSE)
EstimateFunAdd(10,1,FALSE,0.75,2,FALSE,FALSE)
EstimateFunAdd(10,2,FALSE,0.75,2,FALSE,FALSE)

EstimateFunAdd(20,1,TRUE,0.75,2,FALSE,TRUE)
EstimateFunAdd(20,2,TRUE,0.75,2,FALSE,FALSE)
EstimateFunAdd(20,1,FALSE,0.75,2,FALSE,FALSE)
EstimateFunAdd(20,2,FALSE,0.75,2,FALSE,FALSE)

EstimateFunAdd(50,1,TRUE,0.75,2,FALSE,TRUE)
EstimateFunAdd(50,2,TRUE,0.75,2,FALSE,FALSE)
EstimateFunAdd(50,1,FALSE,0.75,2,FALSE,FALSE)
EstimateFunAdd(50,2,FALSE,0.75,2,FALSE,FALSE)



##### CASE N=20
Data <- read.table("datasets/DataCaseN20D1.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN20D1.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN20D1")
gc()
Data <- read.table("datasets/DataCaseN20D2.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN20D2.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN20D2")
gc()
Data <- read.table("datasets/DataCaseN20D1alt.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN20D1alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN20D1alt")
gc()
Data <- read.table("datasets/DataCaseN20D2alt.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN20D2alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN20D2alt")
gc()
##### CASE N=10 
Data <- read.table("datasets/DataCaseN10D2.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN10D2.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN10D2")
gc()
Data <- read.table("datasets/DataCaseN10D1alt.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN10D1alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN10D1alt")
gc()
Data <- read.table("datasets/DataCaseN10D2alt.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN10D2alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN10D2alt")
### N= 50
Data <- read.table("datasets/DataCaseN50D1.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN50D1.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN50D1")
gc()
Data <- read.table("datasets/DataCaseN50D2.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN50D2.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN50D2")
gc()
Data <- read.table("datasets/DataCaseN50D1alt.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN50D1alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN50D1alt")
gc()
Data <- read.table("datasets/DataCaseN50D2alt.txt",header=TRUE)
Data.info <- read.table("datasets/DataCaseN50D2alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"EstimatesCaseN50D2alt")


#### Additional sim
EstimateFun <- function(Data,Data.info,Filename,Return=FALSE){
  N <- Data.info[1,1]
  M <- dim(Data.info)[1]
  
  Estimates.equal <- as.vector(NULL)
  Estimates.propC <- as.vector(NULL)
  Estimates.propCN <- as.vector(NULL)
  
  Validcase.equal <- NA
  Validcase.propC <- NA
  Validcase.propCN <- NA
  for(m in 1:M){
    data.m <- Data[Data[,1]==m ,]
    data.info.m <- Data.info[m,-1]
    
    ### computation of the weights
    N.subsample <- data.info.m$K
    n.k <- t(unique(t(data.info.m[1:N])))
    c.k <- as.vector(table(t(data.info.m[1:N])))
    c.k <- c.k[order(c.k)]
    equal.w <- rep(1/N.subsample,N.subsample)
    propC.w <- c.k/sum(c.k)
    propCN.w <- (c.k*n.k)/sum(c.k*n.k)
    ### estimation  
    BetaH.equal <- rep(0,4) 
    BetaH.propC <- rep(0,4)
    BetaH.propCN <- rep(0,4)
    DH.equal <- matrix(0,4,4)
    DH.propC <- matrix(0,4,4)
    DH.propCN <- matrix(0,4,4)
    SigmaH.equal <- matrix(0,2,2)
    SigmaH.propC <- matrix(0,2,2)
    SigmaH.propCN <- matrix(0,2,2)
    
    #pdm.valid <- 0
    for(i in 1:N.subsample){
      Solution <- ClosedFS(data.m[data.m[,2]==i,])
      #pdm.valid[i] <-Solution$pdm
      
      BetaH.equal <- Solution$FE*equal.w[i] + BetaH.equal
      DH.equal <- Solution$D*equal.w[i] + DH.equal
      SigmaH.equal <- Solution$Sigma*equal.w[i] + SigmaH.equal
      
      BetaH.propC <- Solution$FE*propC.w[i] + BetaH.propC
      DH.propC <- Solution$D*propC.w[i] + DH.propC
      SigmaH.propC <- Solution$Sigma*propC.w[i] + SigmaH.propC
      
      BetaH.propCN <- Solution$FE*propCN.w[i] + BetaH.propCN
      DH.propCN <- Solution$D*propCN.w[i] + DH.propCN
      SigmaH.propCN <- Solution$Sigma*propCN.w[i] + SigmaH.propCN
    }
    #Valid.case[m] <- sum(pdm.valid)/N.subsample
    
    Validcase.equal[m] <- sum(eigen(DH.equal)$values < 0) == 0 
    Validcase.propC[m] <- sum(eigen(DH.propC)$values < 0) == 0
    Validcase.propCN[m] <- sum(eigen(DH.propCN)$values < 0) == 0
    ## computation of the R^2_trial  
    R2.Trial.equal <- R2TrialFun(DH.equal)
    R2.Trial.propC <- R2TrialFun(DH.propC)
    R2.Trial.propCN <- R2TrialFun(DH.propCN)
    
    R2.Ind.equal <- R2IndFun(SigmaH.equal)
    R2.Ind.propC <- R2IndFun(SigmaH.propC)
    R2.Ind.propCN <- R2IndFun(SigmaH.propCN)
    
    
    Estimates.equal.m <- c(as.vector(BetaH.equal), DH.equal[lower.tri(DH.equal,diag=TRUE)], 
                           SigmaH.equal[lower.tri(SigmaH.equal,diag=TRUE)],R2.Trial.equal,R2.Ind.equal)
    Estimates.propC.m <- c(as.vector(BetaH.propC), DH.propC[lower.tri(DH.propC,diag=TRUE)], 
                           SigmaH.propC[lower.tri(SigmaH.propC,diag=TRUE)],R2.Trial.propC,R2.Ind.propC)
    Estimates.propCN.m <- c(as.vector(BetaH.propCN), DH.propCN[lower.tri(DH.propCN,diag=TRUE)], 
                            SigmaH.propCN[lower.tri(SigmaH.propCN,diag=TRUE)],R2.Trial.propCN,R2.Ind.propCN)
    
    Estimates.equal <- rbind(Estimates.equal,Estimates.equal.m)
    Estimates.propC <- rbind(Estimates.propC,Estimates.propC.m)
    Estimates.propCN <- rbind(Estimates.propCN,Estimates.propCN.m)
  }
  Estimates.equal <- cbind(Estimates.equal,Validcase.equal)
  Estimates.propC <- cbind(Estimates.propC,Validcase.propC)
  Estimates.propCN <- cbind(Estimates.propCN,Validcase.propCN)
  rownames(Estimates.equal) <- rownames(Estimates.propC) <- rownames(Estimates.propCN) <- NULL
  colnames(Estimates.equal) <- colnames(Estimates.propC) <- colnames(Estimates.propCN) <- c("BS0","BS1",
                                                                                            "BT0","BT1","dSS","dSa","dST","dSb","daa","daT","dab","dTT","dTb",
                                                                                            "dbb","SigmaSS","SigmaST","SigmaTT","R2Trial","R2Ind","Valid")
  
  write.table(Estimates.equal, paste("Additional Simulation/Estimates/",Filename,".equal",".txt",sep=""), sep="\t",row.names=FALSE)
  write.table(Estimates.propC, paste("Additional Simulation/Estimates/",Filename,".propC",".txt",sep=""), sep="\t",row.names=FALSE)
  write.table(Estimates.propCN, paste("Additional Simulation/Estimates/",Filename,".propCN",".txt",sep=""), sep="\t",row.names=FALSE)
  
  if(Return==TRUE){
    Results <- list(equal = as.data.frame(Estimates.equal), propC = as.data.frame(Estimates.propC),
                    propCN = as.data.frame(Estimates.propCN))
    return(Results)
  }
}

##### CASE N=10 (R=0.25)
Data <- read.table("Additional Simulation/Datasets/N10D1.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/datasets/N10D1.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N10D1")
gc()
Data <- read.table("Additional Simulation/Datasets/N10D2.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N10D2.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N10D2")
gc()
Data <- read.table("Additional Simulation/Datasets/N10D1alt.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N10D1alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N10D1alt")
gc()
Data <- read.table("Additional Simulation/Datasets/N10D2alt.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N10D2alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N10D2alt") ###CHECK!

##### CASE N=20 (R=0.25)
Data <- read.table("Additional Simulation/Datasets/N20D1.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/datasets/N20D1.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N20D1")
g()
Data <- read.table("Additional Simulation/Datasets/N20D2.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N20D2.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N20D2")
gc()
Data <- read.table("Additional Simulation/Datasets/N20D1alt.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N20D1alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N20D1alt")
gc()
Data <- read.table("Additional Simulation/Datasets/N20D2alt.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N20D2alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N20D2alt")


##### CASE N=50 (R=0.25)
Data <- read.table("Additional Simulation/Datasets/N50D1.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/datasets/N50D1.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N50D1")
g()
Data <- read.table("Additional Simulation/Datasets/N50D2.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N50D2.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N50D2")
gc()
Data <- read.table("Additional Simulation/Datasets/N50D1alt.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N50D1alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N50D1alt")
gc()
Data <- read.table("Additional Simulation/Datasets/N50D2alt.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N50D2alt.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N50D2alt")


##### CASE N=10 (R=0.75)
Data <- read.table("Additional Simulation/Datasets/N10D175.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/datasets/N10D175.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N10D175")
gc()
Data <- read.table("Additional Simulation/Datasets/N10D275.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N10D275.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N10D275")
gc()
Data <- read.table("Additional Simulation/Datasets/N10D1alt75.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N10D1alt75.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N10D1alt75")
gc()
Data <- read.table("Additional Simulation/Datasets/N10D2alt75.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N10D2alt75.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N10D2alt75")

##### CASE N=20 (R=0.75)
Data <- read.table("Additional Simulation/Datasets/N20D175.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/datasets/N20D175.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N20D175")
g()
Data <- read.table("Additional Simulation/Datasets/N20D275.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N20D275.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N20D275")
gc()
Data <- read.table("Additional Simulation/Datasets/N20D1alt75.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N20D1alt75.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N20D1alt75")
gc()
Data <- read.table("Additional Simulation/Datasets/N20D2alt75.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N20D2alt75.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N20D2alt75")


##### CASE N=50 (R=0.25)
Data <- read.table("Additional Simulation/Datasets/N50D175.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/datasets/N50D175.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N50D175")
gC()###HERE
Data <- read.table("Additional Simulation/Datasets/N50D275.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N50D275.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N50D275")
gc()
Data <- read.table("Additional Simulation/Datasets/N50D1alt75.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N50D1alt75.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N50D1alt75")
gc()
Data <- read.table("Additional Simulation/Datasets/N50D2alt75.txt",header=TRUE)
Data.info <- read.table("Additional Simulation/Datasets/N50D2alt75.info.txt",header=TRUE)
EstimateFun(Data,Data.info,"N50D2alt75")
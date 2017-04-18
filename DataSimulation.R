library(MASS)
library(gtools)
##### Target parameters #######
## fixed effecs
Beta <- c(450,300,500,500)
## Sigma matrix 
Sigma <- matrix(c(3, c(sqrt(.5*9)), c(sqrt(.5*9)), 3), ncol=2)*100
## D matrix target
D <- matrix(rep(0, times=16), ncol=4)
diag(D) <- c(1, 1, 1, 1)
D[1, 3] <- D[3, 1] <- c(.4)
D[2, 4] <- D[4, 2] <- c(sqrt(.5))
D <- D *1000
###  mean and standard deviation of number of individuals per trial
Mean.n <- 20
Sd.n <- 5
### number of simulations
M <- 1000
#### function to simulate the data
Sim.data <- function(M,N.trial,Beta,D,Sigma,Mean.n,Sd.n,Seed,Filename){
    ## M: number of simulations
    ## N.trial: number of clusters
    ## Beta: Fixed effects
    ## D: random effects variance matrix
    ## Sigma: error variance matrix
    
    Seed <- Seed-1
    Sim.data.all <- c(NULL)
    Cluster.sizes <- c(NULL)
    for(m in 1:M){
      Sim.id <- m
      Seed <- Seed+1
  ## simulation of C_K and n_k
  ### Case N.trial=10 two possibilities
        if(N.trial==10){
          set.seed(Seed)
          K = 2
          c.k <- c(5,5)
          n.k <- c()
          while(length(unique(n.k)) < K){
            n.k <- round(rnorm(length(c.k),Mean.n,Sd.n))
            n.k[n.k < 5] <- 5
          }
        }
### Case N.trial=20 more possibilities
        if(N.trial==20){
          set.seed(Seed)
          K <- sample(2:4,1)
            if(K !=4){
              Min <- 5
              Max <- N.trial-5*(K-1)
          
              all.ck <- combinations((Max-Min+1), K, v=Min:Max,repeats.allowed=TRUE)
              possible.ck <- all.ck[apply(all.ck,1,sum) == N.trial,]
              c.k <- possible.ck[sample(1:dim(possible.ck)[1],1),]
              rm(all.ck)
              rm(possible.ck)
              
            }
            else{c.k=c(5,5,5,5)}
            n.k <- c()
            while(length(unique(n.k)) < K){
              n.k <- round(rnorm(length(c.k),Mean.n,Sd.n))
              n.k[n.k < 5] <- 5

            }
            gc()
        }

### Case N.trial=50 more possibilities
        if(N.trial==50){
          set.seed(Seed)
            K <- sample(2:10,1)
            Min <- 5
            Max <- N.trial-5*(K-1)
            if(K < 10){
              all.ck <- combinations((Max-Min+1), K, v=Min:Max,repeats.allowed=TRUE)
              possible.ck <- all.ck[apply(all.ck,1,sum) == N.trial,]
              c.k <- possible.ck[sample(1:dim(possible.ck)[1],1),]
              rm(all.ck)
              rm(possible.ck)
            }else{c.k <- rep(5,10)}
              n.k <- c()
              while(length(unique(n.k)) < K){
                n.k <- round(rnorm(length(c.k),Mean.n,Sd.n))
                n.k[n.k < 5] <- 5
            }
            
            gc()
        }
        
        N.ind <- sum(n.k*c.k)
        n.trial <- rep(n.k,c.k)
        
        Cluster.sizes <- rbind(Cluster.sizes,c(N.trial,n.trial,K,N.ind))
        Subsample.id <- rep(1:K,times=2*n.k*c.k)
    
  ## data generation
        set.seed(Seed)
        btw.err <- mvrnorm(N.trial, c(0, 0, 0, 0), D)
  
        pat.indicator <- c(0,cumsum(rep(n.k,c.k)))+1
        Sim.data <- as.vector(NULL)
        for (i in 1:N.trial){
          Trial <- i
          pat.id <- rep(pat.indicator[i]:(pat.indicator[i+1]-1),2)
          n1 <- ceiling(n.trial[i]/2)
          n2 <- n.trial[i]-n1
          Treat <- c(rep(-1,n1),rep(1,n2))
          wth.err <- mvrnorm(n.trial[i], c(0, 0), Sigma)
          Surr <- Beta[1] + btw.err[i,1] + (Beta[2] + btw.err[i,2])*Treat + wth.err[,1]  
          True <- Beta[3] + btw.err[i,3] + (Beta[4] + btw.err[i,4])*Treat + wth.err[,2]  
          endpoint <- c(rep(-1,n.trial[i]),rep(1,n.trial[i]))
          data.trial <- cbind(Trial, pat.id,endpoint,rep(Treat,2), c(Surr,True))
          Sim.data <- rbind(Sim.data, data.trial)
        }
        Sim.data <- cbind(Sim.id,Subsample.id,Sim.data)
        Sim.data.all <- rbind(Sim.data.all,Sim.data)
      }
      colnames(Sim.data.all) <- c("Sim.id","Subsample.id","Trial.id","Pat.id","Endpoint","Treat","Outcome")
      colnames(Cluster.sizes) <- c("N.trial",as.character(1:N.trial),"K","N.ind")
      
      write.table(Sim.data.all, paste("datasets/",Filename,".txt",sep=""), sep="\t",row.names=FALSE)
      write.table(Cluster.sizes, paste("datasets/",Filename,".info",".txt",sep=""), sep="\t",row.names=FALSE)
      Results <- list(data=as.data.frame(Sim.data.all),n.trials=as.data.frame(Cluster.sizes))
  }
  


Sim.data(M,10,Beta,D,Sigma,Mean.n,Sd.n,10,"DataCaseN10D1")
Sim.data(M,20,Beta,D,Sigma,Mean.n,Sd.n,20,"DataCaseN20D1")
Sim.data(M,50,Beta,D,Sigma,20,5,50,"DataCaseN50D1")
Sim.data(M,10,Beta,D*0.1,Sigma,Mean.n,Sd.n,11,"DataCaseN10D2")
Sim.data(M,20,Beta,D*0.1,Sigma,Mean.n,Sd.n,21,"DataCaseN20D2")
Sim.data(M,50,Beta,D*0.1,Sigma,20,5,51,"DataCaseN50D2")

#### other case mean(nk)=10 
Sim.data(M,10,Beta,D,Sigma,10,Sd.n,100,"DataCaseN10D1alt")
Sim.data(M,10,Beta,D*0.1,Sigma,10,Sd.n,110,"DataCaseN10D2alt")
Sim.data(M,20,Beta,D,Sigma,10,Sd.n,200,"DataCaseN20D1alt")
Sim.data(M,20,Beta,D*0.1,Sigma,10,Sd.n,220,"DataCaseN20D2alt")
Sim.data(M,50,Beta,D,Sigma,10,5,500,"DataCaseN50D1alt")
Sim.data(M,50,Beta,D*0.1,Sigma,10,5,550,"DataCaseN50D2alt")
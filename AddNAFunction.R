### function to add NA values (MI porpose) to the simulated datasets
AddNAFunction <- function(Data,Data.info,Case){
  M <- max(Data$Sim.id)
  N <- Data.info$N.trial[1]
  for(m in 1:M){
    Data.m <- Data[Data$Sim.id==m,]
    N.subsample <- Data.info$K[m]
    n.max <- max(Data.info[m,2:(N+1)])
    max.id <- max(Data$Pat.id)
    n.k <- t(unique(t(Data.info[m,2:(N+1)])))
    c.k <- as.vector(table(t(Data.info[m,2:(N+1)])))
    c.k <- c.k[order(c.k)]         
    add.n <- rep(n.max-n.k,c.k)
    id.indicator <- max.id
    add.Data <- NULL
    for(i in 1:N){
        if(n.max%%2==0){
          add.n2 <- ceiling(add.n[i]/2)
          add.n1 <- add.n[i] - add.n2
        }else{
          add.n1 <- ceiling(add.n[i]/2)
          add.n2 <- add.n[i] - add.n1
        }
        add.Endpoint <- c(rep(-1,add.n[i]),rep(1,add.n[i]))
        add.Treat <- c(rep(-1,add.n1),rep(1,add.n2))
        add.Sim.id <- m
        add.Subsample.id <- m
        add.Trial.id <- i
        add.Pat.id <- rep((id.indicator+1):(id.indicator+add.n[i]),2)
        add.Outcome <- rep(99999,add.n[i]*2)
        add.DataSub <- cbind(add.Sim.id,add.Subsample.id,add.Trial.id,add.Pat.id,add.Endpoint,add.Treat,add.Outcome)
        if(add.n[i] !=0){
          add.Data <- rbind(add.Data,add.DataSub)
        }
        id.indicator <- id.indicator+add.n[i]
      }
    
    colnames(add.Data) <- names(Data.m) 
    Data.m <- rbind(Data.m, add.Data)
    Data.m <- Data.m[order(Data.m$Trial.id),]
    
    write.table(Data.m, paste("datasets/DataCase",Case,"/","DataCase",Case,".",m,".NA.txt",sep=""), sep="\t",row.names=FALSE)
    }
  }

source('./dataGen.R', echo = TRUE)
source('./statistics.R', echo = TRUE)
K <- 3
repTime <- 2000

doit <- function() {
    # J
    tmp <- lapply(1:K,function(i){
        tmp <- rep(0,n[i]*K)
        dim(tmp) <- c(n[i],K)
        tmp[,i] <- 1/sqrt(n[i])
        tmp
    })
    NEW.J <- do.call(rbind,tmp)
    # C
    tmp <- t(NEW.J)%*%rep(1,sum(n))
    C <- eigen(diag(K)-tmp%*%t(tmp)/sum(n))$vectors[,-K]
   
    #SNR 
    tmpMu <- do.call(cbind,mu)
    for(i in 1:K){
        tmpMu[,i] <- tmpMu[,i]*sqrt(n[i])
    }
    tmpMuF <- tmpMu%*%C
        
    tmpSigEig <- eigen(Sigma)$values#[-c(1,2)]
    tmpCon <- sqrt(SNR*sqrt(sum(tmpSigEig^2))/sum(tmpMuF^2))
    for(i in 1:length(mu)){
        mu[[i]] <- mu[[i]]*tmpCon
    }
    
    

    pb <- txtProgressBar(style = 3)
    
    jCX <- rep(0, repTime)
    jSC <- rep(0, repTime)
    jHBWW <- rep(0,repTime)
    jZGZ <- rep(0,repTime)
    jAsy <- rep(0,repTime)
    for (myIterator in 1:repTime) {
        X <- lapply(1:K, function(k) {
            dataGen(n[k], p, mu[[k]], Sigma)
        })
        CXt <- CXtest(n, p, K, X)
        jCX[myIterator] <- (CXt <= 0.05)
        
        SCt <- SCtest(n, p, K, X, NEW.J = NEW.J)
        jSC[myIterator] <- (SCt <= 0.05)
        
        HBWWt<- HBWWtest(n,p,K,X)
        jHBWW[myIterator] <- (HBWWt <= 0.05)
        
        ZGZt<- ZGZtest(n,p, K, X, NEW.J= NEW.J)
        jZGZ[myIterator] <- (ZGZt <= 0.05)
        
        myAsy <- newTest(n,p, K, X, NEW.J= NEW.J, C=C)
        jAsy[myIterator] <- (myAsy <= 0.05)
        
        setTxtProgressBar(pb, myIterator / repTime)
    }
    close(pb)
    return(data.frame(
        SNR=SNR,
        SC = mean(jSC),
        CX = mean(jCX),
        HBWW = mean(jHBWW),
        ZGZ = mean(jZGZ),
        Asy = mean(jAsy)
    ))
}

######################################################
tmpList<-NULL
for(SNR in seq(0,1)){
    K = 3
    n <- c(15, 15, 15)
    p = 200
    Sigma <- diag(p)
    Sigma[1,1]<- 3*p
    Sigma[2,2]<- 2*p
    Sigma[3,3]<- 1*p
    
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"1p50.csv",row.names=FALSE)
 
#########################################################################

tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(10, 10, 10)
    p = 75
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"1p75.csv",row.names=FALSE)
 
#########################################################################

tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(10, 10, 10)
    p = 100
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"1p100.csv",row.names=FALSE)

 
#########################################################################



##########  table 2 #################################3
######################################################
tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(25, 25, 25)
    p = 100
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"2p100.csv",row.names=FALSE)
 
#########################################################################

tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(25, 25, 25)
    p = 150
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"2p150.csv",row.names=FALSE)
 
#########################################################################

tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(25, 25, 25)
    p = 200
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"2p200.csv",row.names=FALSE)

 
#########################################################################

################## table 3
######################################################
tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(10, 10, 10)
    p = 50
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"3p50.csv",row.names=FALSE)
 
#########################################################################

tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(10, 10, 10)
    p = 75
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"3p75.csv",row.names=FALSE)
 
#########################################################################

tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(10, 10, 10)
    p = 100
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"3p100.csv",row.names=FALSE)

 
#########################################################################



##########  table 4 #################################3
######################################################
tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(25, 25, 25)
    p = 100
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"4p100.csv",row.names=FALSE)
 
#########################################################################

tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(25, 25, 25)
    p = 150
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"4p150.csv",row.names=FALSE)
 
#########################################################################

tmpList<-NULL
for(SNR in seq(0,10)){
    B = 100
    # sample number
    K = 3
    n <- c(25, 25, 25)
    p = 200
    Sigma <- diag(p)
    Sigma[1,1]<- 1.5*p
    Sigma[2,2]<- p
    
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"4p200.csv",row.names=FALSE)

 
#########################################################################










a<-read.csv("4p100.csv")
b<-read.csv("4p150.csv")
c<-read.csv("4p200.csv")

tmp<-merge(a,b,by.x="SNR",by.y="SNR")
tmp<-merge(tmp,c,by.x="SNR",by.y="SNR")

library(xtable)
zzz<-xtable(tmp,digits=3)
print(zzz,include.rownames=FALSE)
#########################################################################
##### plot(ecdf(jjj))
#
# TheoryCDF <- function(x){
#     tmp <- uniroot(function(t){gamma(t)-((K-1)/2)},c(1,100))$root
#     exp(-tmp*exp(-x/K))
# }
#
# curve(TheoryCDF,from=-10,to=30,add=TRUE)

source('./dataGen.R', echo = TRUE)
source('./statistics.R', echo = TRUE)
repTime <- 100

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
        
    tmpSigEig <- eigen(Sigma)$values[-c(1,2)]
    tmpCon <- sqrt(SNR*sqrt(sum(tmpSigEig^2))/sum(tmpMuF^2))
    for(i in 1:length(mu)){
        mu[[i]] <- mu[[i]]*tmpCon
    }
    
    

    pb <- txtProgressBar(style = 3)
    
    jCX <- rep(0, repTime)
    jNEW <- rep(0, repTime)
    jSC <- rep(0, repTime)
    for (myIterator in 1:repTime) {
        X <- lapply(1:K, function(k) {
            dataGen(n[k], p, mu[[k]], Sigma)
        })
        theCXstat <- CXstat(n, p, K, X)
        
        tmp <- NEWstat(n, p, K, X,NEW.J=NEW.J,C=C)
        theNEWstat <- tmp$stat
        Zinv <- tmp$Zinv
        
        tmp <- SCstat(n, p, K, X,NEW.J=NEW.J)
        theSCstat <- tmp$stat
        myGram <- tmp$myGram
        
        tmpNEWstat <- rep(0, B)
        tmpCXstat <- rep(0, B)
        tmpSCstat <- rep(0, B)
        for (xxx in 1:B) {
            theOrder <- sample.int(sum(n))
            thePer <- myPer(n, p, K, X, theOrder)
            tmpCXstat[xxx] <- CXstat(n, p, K, thePer)
            tmpNEWstat[xxx] <-
                NEWstat(n, p, K, thePer, Zinv[theOrder, theOrder],NEW.J=NEW.J,C=C)
            tmpSCstat[xxx] <-
                SCstat(n, p, K, thePer, myGram = myGram[theOrder, theOrder],NEW.J=NEW.J)
        }
        if ((sum(tmpCXstat >= theCXstat) + 1) / (B + 1) <= 0.05)
            jCX[myIterator] <- 1
        
        if ((sum(tmpNEWstat >= theNEWstat) + 1) / (B + 1) <= 0.05)
            jNEW[myIterator] <- 1
        
        if ((sum(tmpSCstat >= theSCstat) + 1) / (B + 1) <= 0.05)
            jSC[myIterator] <- 1
        
        setTxtProgressBar(pb, myIterator / repTime)
    }
    close(pb)
    return(data.frame(
        SNR=SNR,
        CX = mean(jCX),
        SC = mean(jSC),
        NEW = mean(jNEW)
    ))
}

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










a<-read.csv("1p50.csv")
b<-read.csv("1p75.csv")
c<-read.csv("1p100.csv")

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
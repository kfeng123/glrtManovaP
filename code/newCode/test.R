
K <- 3
repTime <- 2000


source('./dataGen.R', echo = TRUE)
source('./statistics.R', echo = TRUE)

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

####################################################################
# a<-read.csv("4p100.csv")
# b<-read.csv("4p150.csv")
# c<-read.csv("4p200.csv")

# tmp<-merge(a,b,by.x="SNR",by.y="SNR")
# tmp<-merge(tmp,c,by.x="SNR",by.y="SNR")

# library(xtable)
# zzz<-xtable(tmp,digits=3)
# print(zzz,include.rownames=FALSE)
#########################################################################
##### plot(ecdf(jjj))
#
# TheoryCDF <- function(x){
#     tmp <- uniroot(function(t){gamma(t)-((K-1)/2)},c(1,100))$root
#     exp(-tmp*exp(-x/K))
# }
#
# curve(TheoryCDF,from=-10,to=30,add=TRUE)

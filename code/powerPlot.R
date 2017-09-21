source('./dataGen.R', echo = TRUE)
source('./statistics.R', echo = TRUE)
repTime <- 100

rock <- function() {
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
    dim(C) <-c(2,1)
   
    #SNR 
    tmpMu <- do.call(cbind,mu)
    for(i in 1:K){
        tmpMu[,i] <- tmpMu[,i]*sqrt(n[i])
    }
    tmpMuF <- tmpMu%*%C
        
    tmpSigEig <- eigen(Sigma)$values[-1]
    tmpCon <- sqrt(SNR*sqrt(sum(tmpSigEig^2))/sum(tmpMuF^2))
    for(i in 1:length(mu)){
        mu[[i]] <- mu[[i]]*tmpCon
    }
    
    

    pb <- txtProgressBar(style = 3)
    
    jCQ <- rep(0, repTime)
    jSR <- rep(0, repTime)
    jCLX <- rep(0, repTime)
    jFZWZ <- rep(0, repTime)
    jNEW <- rep(0, repTime)
    for (myIterator in 1:repTime) {
        X <- lapply(1:K, function(k) {
            dataGen(n[k], p, mu[[k]], Sigma)
        })
        
        # Chen Qin and gram matrix
        tmp <- CQstat(n,p,X)
        theCQstat <- tmp$stat
        myGram <- tmp$myGram
        
        # Srivastava 2009
        theSRstat <- SRstat(n,p,X)$stat
        
        # CLX, it doesn't have fast algorithm
        theCLXstat <- CXstat(n, p, K, X)
        
        # FZWZ
        theFZWZstat <- FZWZstat(n,p,X)$stat
        
        # new
        tmp <- NEWstat(n, p, K, X,NEW.J=NEW.J,C=C)
        theNEWstat <- tmp$stat
        Zinv <- tmp$Zinv
        
        
        # permutation
        tmpCQstat <- rep(0, B)
        tmpSRstat <- rep(0, B)
        tmpCLXstat <- rep(0, B)
        tmpFZWZstat <- rep(0, B)
        tmpNEWstat <- rep(0, B)
        for (xxx in 1:B) {
            theOrder <- sample.int(sum(n))
            thePer <- myPer(n, p, K, X, theOrder)
            
            tmpCQstat[xxx] <-CQstat(n, p, thePer, myGram = myGram[theOrder, theOrder])
            tmpSRstat[xxx] <- SRstat(n, p, thePer)
            tmpFZWZstat[xxx] <- FZWZstat(n, p, thePer)
            tmpCLXstat[xxx] <- CXstat(n, p, K, thePer)
            tmpNEWstat[xxx] <-
                NEWstat(n, p, K, thePer, Zinv[theOrder, theOrder],NEW.J=NEW.J,C=C)
        }
        
        if ((sum(tmpCQstat >= theCQstat) + 1) / (B + 1) <= 0.05)
            jCQ[myIterator] <- 1
        
        if ((sum(tmpSRstat >= as.numeric(theSRstat)) + 1) / (B + 1) <= 0.05)
            jSR[myIterator] <- 1
        
        if ((sum(tmpFZWZstat >= theFZWZstat) + 1) / (B + 1) <= 0.05)
            jFZWZ[myIterator] <- 1
        
        if ((sum(tmpCLXstat >= theCLXstat) + 1) / (B + 1) <= 0.05)
            jCLX[myIterator] <- 1
        
        if ((sum(tmpNEWstat >= theNEWstat) + 1) / (B + 1) <= 0.05)
            jNEW[myIterator] <- 1
        
        
        setTxtProgressBar(pb, myIterator / repTime)
    }
    close(pb)
    return(data.frame(
        Rho=Rho,
        CQ = mean(jCQ),
        SR = mean(jSR),
        CLX = mean(jCLX),
        FZWZ = mean(jFZWZ),
        NEW = mean(jNEW)
    ))
}

tmpList<-NULL
for(Rho in seq(0,0.9,0.1)){
    SNR=5
    B = 100
    # sample number
    K = 2
    n <- c(20, 20)
    p = 150
    Sigma <- rep(Rho,p^2)
    dim(Sigma) <- c(p,p)
    diag(Sigma) <- rep(1,p)
    
    mu <- list(c(rep(1, p/2),rep(-1,p/2)),rep(0, p))
    
    (tmpFram<-rock())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"tu.csv",row.names=FALSE)


tCQ <- data.frame(
    rho=theOut$Rho,
    power=theOut$CQ,
    method="CQ"
)
tSR <- data.frame(
    rho=theOut$Rho,
    power=theOut$SR,
    method="SR"
)
tCLX <- data.frame(
    rho=theOut$Rho,
    power=theOut$CLX,
    method="CLX"
    
)
tFZWZ <- data.frame(
    rho=theOut$Rho,
    power=theOut$FZWZ,
    method="FZWZ"
)
tNEW <- data.frame(
    rho=theOut$Rho,
    power=theOut$NEW,
    method="NEW"
)
df=rbind(tCQ,tSR,tCLX,tFZWZ,tNEW)

library(ggplot2)
library(latex2exp)
myPlot <- ggplot(df,aes(rho,power,colour=method))+
    geom_point()+
    geom_line()+
    ylab("Empirical power")+
    xlab(TeX('$\\rho$'))+
    scale_x_continuous(breaks=seq(0,0.9,0.1))+
    guides(colour=guide_legend(title=NULL))+
    theme_bw()

ggsave("figure1.eps",myPlot)

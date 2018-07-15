
    n <- c(35, 35, 35)
    p = 1000
repTime <- 5000

tmpList<-NULL
for(Rho in seq(0,0.9,0.1)){
    
    SNR <- 5
    
    Sigma <- rep(Rho,p^2)
    dim(Sigma) <- c(p,p)
    diag(Sigma) <- rep(1,p)
    
    
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
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
        
    
    oSig <- eigen(Sigma, symmetric= TRUE)
    tmpSigEig <- oSig$values[-1]
    tmpCon <- sqrt(SNR*sqrt(sum(tmpSigEig^2))/sum(tmpMuF^2))
    for(i in 1:length(mu)){
        mu[[i]] <- mu[[i]]*tmpCon
    }
    
    
    
    tmpFram<-doit()
    tmpFram$Rho <-  Rho
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)


write.csv(theOut,"tu.csv",row.names=FALSE)

theOut <- read.csv('tu.csv')

tSC <- data.frame(
    rho=theOut$Rho,
    power=theOut$SC,
    method="Sc"
)
tCX <- data.frame(
    rho=theOut$Rho,
    power=theOut$CX,
    method="CX"
)
tHBWW <- data.frame(
    rho=theOut$Rho,
    power=theOut$HBWW,
    method="HBWW"
    
)
tZGZ <- data.frame(
    rho=theOut$Rho,
    power=theOut$ZGZ,
    method="ZGZ"
)
tLFD <- data.frame(
    rho=theOut$Rho,
    power=theOut$Asy,
    method="LFD"
)
df=rbind(tSC,tCX,tHBWW,tZGZ,tLFD)

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

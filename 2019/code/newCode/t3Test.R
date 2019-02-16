    n <- c(30, 30, 30)
    p = 800

tmpList<-NULL
for(SNR in c(0,1,2,3,4)){
    
    Sigma <- diag(p)
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
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
    tmpSigEig <- oSig$values#[-c(1,2)]
    tmpCon <- sqrt(SNR*sqrt(sum(tmpSigEig^2))/sum(tmpMuF^2))
    for(i in 1:length(mu)){
        mu[[i]] <- mu[[i]]*tmpCon
    }
    
    
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t3l.csv",row.names=FALSE)

tmpList<-NULL
for(SNR in c(0,1,2,3,4)){
    
    Sigma <- diag(p)
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
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
    tmpSigEig <- oSig$values#[-c(1,2)]
    tmpCon <- sqrt(SNR*sqrt(sum(tmpSigEig^2))/sum(tmpMuF^2))
    for(i in 1:length(mu)){
        mu[[i]] <- mu[[i]]*tmpCon
    }
    
    
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t3lb.csv",row.names=FALSE)



tmpList<-NULL
for(SNR in c(0,1,2,3,4)){
    
    Sigma <- diag(p)
    for(i in 1:p)for(j in 1:p){
        Sigma[i,j] <- 0.6^(abs(i-j))
    }
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    
    
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
    tmpSigEig <- oSig$values#[-c(1,2)]
    tmpCon <- sqrt(SNR*sqrt(sum(tmpSigEig^2))/sum(tmpMuF^2))
    for(i in 1:length(mu)){
        mu[[i]] <- mu[[i]]*tmpCon
    }
    
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t3r.csv",row.names=FALSE)


tmpList<-NULL
for(SNR in c(0,1,2,3,4)){
    
    Sigma <- diag(p)
    for(i in 1:p)for(j in 1:p){
        Sigma[i,j] <- 0.6^(abs(i-j))
    }
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
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
    tmpSigEig <- oSig$values#[-c(1,2)]
    tmpCon <- sqrt(SNR*sqrt(sum(tmpSigEig^2))/sum(tmpMuF^2))
    for(i in 1:length(mu)){
        mu[[i]] <- mu[[i]]*tmpCon
    }
    
    
    
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t3rb.csv",row.names=FALSE)
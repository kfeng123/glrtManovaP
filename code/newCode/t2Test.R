
    K = 3
    n <- c(20, 20, 20)
    p = 300


tmpList<-NULL
for(SNR in c(0,0.2,0.4,0.8,1.6,3.2)){
    K = 3
    n <- c(20, 20, 20)
    p = 300
    
    tmp <- rnorm(p*p)
    dim(tmp) <- c(p,p)
    myU<-svd(tmp)$u
    tmp2 <- diag(p)
    tmp2[1,1] <- 3*p
    tmp2[2,2] <- 2*p
    tmp2[3,3] <- 1*p
    Sigma <- myU %*% tmp2 %*% t(myU)
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t2l.csv",row.names=FALSE)

tmpList<-NULL
for(SNR in c(0,0.2,0.4,0.8,1.6,3.2)){
    K = 3
    n <- c(20, 20, 20)
    p = 300
    
    
    tmp <- rnorm(p*p)
    dim(tmp) <- c(p,p)
    myU<-svd(tmp)$u
    tmp2 <- diag(p)
    tmp2[1,1] <- 3*p
    tmp2[2,2] <- 2*p
    tmp2[3,3] <- 1*p
    Sigma <- myU %*% tmp2 %*% t(myU)
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t2lb.csv",row.names=FALSE)



tmpList<-NULL
for(SNR in c(0,0.2,0.4,0.8,1.6,3.2)){
    
    tmp <- rnorm(p*p)
    dim(tmp) <- c(p,p)
    myU<-svd(tmp)$u
    tmp2 <- diag(p)
    tmp2[1,1] <- p
    tmp2[2,2] <- p
    tmp3 <- diag(p)
    for(i in 1:p)for(j in 1:p){
        tmp3[i,j] <- rbinom(1,1,0.01)
    }
    Sigma <- myU %*% tmp2 %*% t(myU) + tmp3 %*% t(tmp3)
    
    
    
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
write.csv(theOut,"t2r.csv",row.names=FALSE)





tmpList<-NULL
for(SNR in c(0,0.2,0.4,0.8,1.6,3.2)){
    
    tmp <- rnorm(p*p)
    dim(tmp) <- c(p,p)
    myU<-svd(tmp)$u
    tmp2 <- diag(p)
    tmp2[1,1] <- p
    tmp2[2,2] <- p
    tmp3 <- diag(p)
    for(i in 1:p)for(j in 1:p){
        tmp3[i,j] <- rbinom(1,1,0.01)
    }
    Sigma <- myU %*% tmp2 %*% t(myU) + tmp3 %*% t(tmp3)
    
    
    
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
write.csv(theOut,"t2rb.csv",row.names=FALSE)
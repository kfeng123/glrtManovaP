tmpList<-NULL
for(SNR in seq(0,3)){
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
for(SNR in seq(0,3)){
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
for(SNR in seq(0,3)){
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
write.csv(theOut,"t2r.csv",row.names=FALSE)


tmpList<-NULL
for(SNR in seq(0,3)){
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
write.csv(theOut,"t2rb.csv",row.names=FALSE)
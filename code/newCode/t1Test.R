tmpList<-NULL
for(SNR in seq(0,3)){
    K = 3
    n <- c(20, 20, 20)
    p = 300
    Sigma <- diag(p)
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t1l.csv",row.names=FALSE)

tmpList<-NULL
for(SNR in seq(0,3)){
    K = 3
    n <- c(20, 20, 20)
    p = 300
    Sigma <- diag(p)
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t1lb.csv",row.names=FALSE)



tmpList<-NULL
for(SNR in seq(0,3)){
    K = 3
    n <- c(20, 20, 20)
    p = 300
    Sigma <- diag(p)
    for(i in 1:p)for(j in 1:p){
        Sigma[i,j] <- 0.6^(abs(i-j))
    }
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t1r.csv",row.names=FALSE)


tmpList<-NULL
for(SNR in seq(0,3)){
    K = 3
    n <- c(20, 20, 20)
    p = 300
    Sigma <- diag(p)
    for(i in 1:p)for(j in 1:p){
        Sigma[i,j] <- 0.6^(abs(i-j))
    }
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t1rb.csv",row.names=FALSE)
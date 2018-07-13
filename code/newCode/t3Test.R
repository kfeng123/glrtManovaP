    n <- c(30, 30, 30)
    p = 800

tmpList<-NULL
for(SNR in c(0,0.2,0.4,0.8,1.6,3.2)){
    
    Sigma <- diag(p)
    mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    #mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t3l.csv",row.names=FALSE)

tmpList<-NULL
for(SNR in c(0,0.2,0.4,0.8,1.6,3.2)){
    
    Sigma <- diag(p)
    #mu <- list(rep(1, p), rep(-1, p), rep(0, p))
    mu <- list(c(rep(1, p/5),rep(0,4*p/5)), c(rep(0, p/5),rep(1, p/5),rep(0, 3*p/5)), rep(0, p))
    
    (tmpFram<-doit())
    tmpList<-c(tmpList,list(tmpFram))
}

theOut<-do.call(rbind,tmpList)
write.csv(theOut,"t3lb.csv",row.names=FALSE)



tmpList<-NULL
for(SNR in c(0,0.2,0.4,0.8,1.6,3.2)){
    
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
write.csv(theOut,"t3r.csv",row.names=FALSE)


tmpList<-NULL
for(SNR in c(0,0.2,0.4,0.8,1.6,3.2)){
    
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
write.csv(theOut,"t3rb.csv",row.names=FALSE)
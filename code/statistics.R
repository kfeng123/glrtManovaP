
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




# Tony Cai and Yin Xia
CXstat <- function(n,p,K,X){
    invSigma <- solve(Sigma)
    transformedXbar <- lapply(X, function(X){
        invSigma %*% colMeans(X)
    })
    
    # temp <- sapply(1:p,function(i){
    #     mySum <-0
    #     for(j in 1:(K-1))for(l in (j+1):K){
    #         mySum <- mySum + n[j]*n[l]/(n[j]+n[l])*(transformedXbar[[j]][i]-transformedXbar[[l]][i])^2/invSigma[i,i]
    #     }
    #     mySum
    # })
    # 
    mySum <- rep(0,p)
    diagInvSigma <- diag(invSigma)
    for(j in 1:(K-1))for(l in (j+1):K){
        mySum <- mySum + n[j]*n[l]/(n[j]+n[l])*(transformedXbar[[j]]-transformedXbar[[l]])^2/diagInvSigma
    }
    
    max(mySum)
}

# Schott
SCstat <- function(n,p,K,X,myGram=NULL){
    out <- list()
    if(is.null(myGram)){
        bindX <- do.call(rbind, X)
        myGram <- bindX%*%t(bindX)
        out$myGram <- myGram
    }
    tmp1 <- sum(diag(myGram))
    tmp2 <- sum(diag(t(NEW.J)%*%myGram%*%NEW.J))
    tmp3 <- sum(myGram)/sum(n)
    traceE <-tmp1-tmp2
    traceH <- tmp2-tmp3
    out$stat <- (traceH/(K-1)-traceE/(sum(n)-K))/sqrt(sum(n)-1)
    return(out)
}


# the New stat


NEWstat <- function(n,p,K,X,Zinv=NULL){
    out <- list()
    if(is.null(Zinv)){
        bindX <- do.call(rbind, X)
        Zinv <- solve(bindX%*%t(bindX))
        out$Zinv <- Zinv
    }
    tmp <- t(C)%*%solve(t(NEW.J)%*%Zinv%*%NEW.J)%*%C
    out$stat <- eigen(tmp)$values[1]
    #out$stat <- sum(eigen(tmp)$values)
    return(out)
}
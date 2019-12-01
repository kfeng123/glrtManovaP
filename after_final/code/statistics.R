library(MASS)


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
SCstat <- function(n,p,K,X,myGram=NULL,NEW.J){
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

# Chen and Qin (2010)
CQstat <- function(n,p,X,myGram=NULL){
    out <- list()
    if(is.null(myGram)){
        bindX <- do.call(rbind, X)
        myGram <- bindX%*%t(bindX)
        out$myGram <- myGram
    }
    #sample 1
    tmp1=sum(myGram[1:n[[1]],1:n[[1]]])-sum(diag(myGram[1:n[[1]],1:n[[1]]]))
    #sample 2
    tmp2=sum(myGram[(n[[1]]+1):n[[2]],(n[[1]]+1):n[[2]]])-
        sum(diag(myGram[(n[[1]]+1):n[[2]],(n[[1]]+1):n[[2]]]))
    #cross term
    tmp3=sum(myGram[1:n[[1]],(n[[1]]+1):n[[2]]])
    out$stat <- tmp1/n[[1]]/(n[[1]]-1)+tmp2/n[[2]]/(n[[2]]-1)-2*tmp3/n[[1]]/n[[2]]
    return(out)
}

# Srivastava (2007)
SRstat <- function(n,p,X){
    out <- list()
    mean1 <- colMeans(X[[1]])
    var1 <- var(X[[1]])
    mean2 <- colMeans(X[[2]])
    var2 <- var(X[[2]])
    myVar <- var1*(n[[1]]-1)+var2*(n[[2]]-1)
    out$stat<-t(mean1-mean2) %*% ginv(myVar) %*% (mean1-mean2)
    return(out)
}

# Feng et al (2014)
FZWZstat <- function(n,p,X){
    out <- list()
    mean1 <- colMeans(X[[1]])
    var1 <- apply(X[[1]],2,var)
    mean2 <- colMeans(X[[2]])
    var2 <- apply(X[[2]],2,var)
    myVar <- var1/n[[1]]+var2/n[[2]]
    out$stat<-sum((mean1-mean2)^2/myVar)
    return(out)
}





# the New stat


NEWstat <- function(n,p,K,X,Zinv=NULL,NEW.J,C){
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
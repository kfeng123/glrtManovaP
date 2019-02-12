library(MASS)


# Tony Cai and Yin Xia
CXtest <- function(n,p,K,X){
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
    MOMG <- max(mySum)
    tmp <- MOMG - 4*log(p) + 2* log(log(p))
    return(
        1-exp(-1/sqrt(pi)*2/3*exp(-tmp/4))
    )
}

# Schott
SCtest <- function(n,p,K,X,myGram=NULL,NEW.J){
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
    stat <- (traceH/(K-1)-traceE/(sum(n)-K))/sqrt(sum(n)-1)
    ########################### 
    # eigenvalues of \hat{\Sigma}
    N <- sum(n)
    nn <- N-K
    estEigen <- eigen(
    (diag(N)- NEW.J%*%t(NEW.J)) %*% bindX %*% t(bindX) %*% (diag(N)- NEW.J%*%t(NEW.J)),
    symmetric = TRUE , only.values =TRUE
    )$values[1:nn]/nn
    trEst <- sum(estEigen)
    trEst2 <- sum(estEigen^2)-trEst^2/nn
    ###########################
    trEst2 <- trEst2*nn^2/(nn+2)/(nn-1)
    theVar <- 2/(K-1)/(sum(n)-K)*trEst2
    
    return(
        1-
    pnorm(stat/sqrt(theVar))
    )
}


# Jiang Hu, Zhidong Bai, Chen Wang, Wei Wang 2017 AISM
HBWWtest <-function(n,p,K,X){
    theStat <- 0
    Xbar <- lapply(X,colMeans)
    Xvar <- lapply(X,function(joj){
        sum(apply(joj,2,var))
    })
    for(i in 1:(K-1)) for(j in (i+1):K){
        theStat <- theStat + t(Xbar[[i]] - Xbar[[j]])%*%(Xbar[[i]] - Xbar[[j]])
    }
    for(i in 1:K){
        theStat <- theStat - (K-1)*  Xvar[[i]] /n[[i]]
    }
    
    theVar <- 0
    for(i in 1:K) {
         tmp <- scale(X[[i]],center=TRUE,scale=FALSE)
         tmp <- tmp%*% t(tmp)/(n[[i]]-1)
         
         theVar <- theVar +  2*(K-1)^2/n[[i]]/(n[[i]]-1) * (n[[i]]-1)^2/(n[[i]]+1)/(n[[i]]-2)*(
         sum(tmp^2)-Xvar[[i]]^2/(n[[i]]-1)
         )
    }
    for(i in 1:(K-1)) for(j in (i+1):K){
        theVar <- theVar + 4/n[[i]]/n[[j]]*sum(diag(
           var(X[[i]]) %*% var(X[[j]])
        ))
    }
    
    1-pnorm(theStat/sqrt(theVar))
    
}


#Jin-Ting Zhang, Jia Guo, Bu Zhou
ZGZtest <- function(n,p,K,X,NEW.J){
    tmp <- lapply(X,colMeans)
    Mhat <- do.call(rbind,tmp)
    tG <- rep(0,(K-1)* K)
    dim(tG) <- c(K-1,K)
    for(i in 1:(K-1)){
        tG[i,i] <- 1
        tG[i,i+1] <- -1
    }
    theD <- diag(1/n)
    Tn <- sum(diag(
        t(Mhat)%*% t(tG) %*% solve(tG %*% theD %*% t(tG)) %*% tG %*% Mhat
    ))
    
    # eigenvalues
    N <- sum(n)
    nn <- N-K
    bindX <- do.call(rbind, X)
    estEigen <- eigen(
    (diag(N)- NEW.J%*%t(NEW.J)) %*% bindX %*% t(bindX) %*% (diag(N)- NEW.J%*%t(NEW.J)),
    symmetric = TRUE , only.values =TRUE
    )$values[1:nn]/nn
    
    tr1 <- nn*(nn+1)/(nn-1)/(nn+2)*(sum(estEigen)^2-2/(nn+1)*sum(estEigen^2))
    tr2 <- nn^2 / (nn-1)/ (nn+2)*(sum(estEigen^2)-1/nn*sum(estEigen)^2)
    betahat <- tr2/sum(estEigen) 
    dhat <- (K-1)*tr1/tr2
    1- pchisq(Tn/betahat,df=dhat)
}



# the New stat

# generate wishart
myWishart <- function (n,df,Sigma){
    if(df >= dim(Sigma)[1]){
        return(rWishart(n,df,Sigma))
    }
    else{
        tmp <- list()
        for ( i in 1:n){
            tmp2 <- rnorm(df*(dim(Sigma)[1]))
            dim(tmp2) <- c(df,dim(Sigma)[1])
            tmp[[i]] <- t(tmp2) %*% tmp2
        }
        myRes <- do.call(c,tmp)
        dim(myRes) <- c(dim(Sigma)[1],dim(Sigma)[1],n)
        return(myRes)
    }
}

# reference distributions
jojo <- 1000
temp <- lapply(1:jojo, function(i){
    temp <- rnorm((K-1)^2)
    dim(temp) <- c(K-1,K-1)
    for (j1 in 1:(K-1))for (j2 in 1:(K-1)){
        if(j1>j2) {
            temp[j1,j2] <- temp[j2,j1]
        }
    }
    for (j in 1:(K-1)){
        temp[j,j] <- temp[j,j]* sqrt(2)
    }
    return(temp)
})
refW <- rep(0, (K-1)*(K-1)*jojo)
dim(refW) <- c(K-1 , K-1 ,jojo)
for(i in 1: jojo){
    refW[,,i] <- temp[[i]]
}
    


newTest <- function(n,p,K,X,NEW.J, C){
    
    # (X^\top X)^{-1}
    bindX <- do.call(rbind, X)
    Zinv <- solve(bindX%*%t(bindX))
        
    tmp <- t(C)%*%solve(t(NEW.J)%*%Zinv%*%NEW.J)%*%C
    Tx <- eigen(tmp)$values[1]
    
    
    # eigenvalues of \hat{\Sigma}
    N <- sum(n)
    nn <- N-K
    estEigen <- eigen(
    (diag(N)- NEW.J%*%t(NEW.J)) %*% bindX %*% t(bindX) %*% (diag(N)- NEW.J%*%t(NEW.J)),
    symmetric = TRUE , only.values =TRUE
    )$values[1:nn]/nn
    
    # detection of spiked covariance
    tau=5
    if(nn*estEigen[1]/sum(estEigen)< tau){
        # stanrdardized statistic
        trEst <- sum(estEigen)
        trEst2 <- sum(estEigen^2)-trEst^2/nn
        standardizedT <- (Tx-(trEst-nn* trEst2 / trEst )) / sqrt( trEst2 )
        # reference distribution
        tmpRef <- apply(refW,3, function(ele){
            eigen(ele,symmetric = TRUE, only.values = TRUE)$values[1]
        })
        return(mean(tmpRef > standardizedT))
    }
    else
    {
        # estimation of r
        for(i in 1:floor(sqrt(nn))){
            if(nn*estEigen[i+1]/sum(estEigen[(i+1):nn])< tau){
                myR <- i
                break
            }
            myR <- i
        }
        trEst <- 1/(1-myR/nn) * sum(estEigen[(myR+1):nn])
        trEst2 <- sum((estEigen[(myR+1):nn]-rep(trEst/nn,nn-myR))^2)
        # standardized statistic
        standardizedT <- 
            (Tx - (
            (1+myR/nn) * trEst - nn * trEst2 / trEst
        )) / sqrt( myR/nn^2 * trEst^2 + trEst2  )
        # reference distribution
        tempa <- trEst/nn/sqrt(myR*trEst^2/nn^2 + trEst2  )
        tempb <- sqrt(trEst2  ) /sqrt(myR*trEst^2/nn^2 + trEst2  )
        
        tmpRef <- apply(
            tempa *( myWishart(jojo, myR, diag(K-1)) -myR * array(rep(diag(K-1),jojo),c(K-1,K-1,jojo)))  + tempb * refW,
                        3,
                        function(ele){
                        eigen(ele,symmetric = TRUE, only.values = TRUE)$values[1]
                        })
        return(mean(tmpRef > standardizedT))
    }
}
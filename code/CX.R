dataGen <- function(n, p, mu, Sigma) {
    temp <- rnorm(n * p)
    dim(temp) <- c(n, p)
    temp2 <- eigen(Sigma, symmetric = TRUE)
    temp %*% diag(temp2$values ^ (1 / 2)) %*% t(temp2$vectors) + outer(rep(1,n),mu)
}

# sample number
K=5

B=100

n <- c(20,20,20,20,20)
p=30
mu <- list(rep(0,p),rep(0,p),rep(0,p),rep(0,p),rep(0,p))
Sigma <- diag(p)



CXstat <- function(n,p,K,X){
    transformedXbar <- lapply(X, function(X){
        solve(Sigma)%*%colMeans(X)
    })
    
    temp <- sapply(1:p,function(i){
        mySum <-0
        for(j in 1:(K-1))for(l in (j+1):K){
            mySum <- mySum + n[j]*n[l]/(n[j]+n[l])*(transformedXbar[[j]][i]-transformedXbar[[l]][i])^2/solve(Sigma)[i,i]
        }
        mySum
    })
    max(temp)
}

myPer <- function(n,p,K,X){
    tmp <- do.call(rbind,X)
    tmp <- tmp[sample.int(nrow(tmp)),]
    tmp2 <- rep(1:K,times=n)  
    split(data.frame(tmp),tmp2)
}

jjj <- rep(0,1000)
for(myIterator in 1:1000){
    X <- lapply(1:K,function(k){
        dataGen(n[k],p,mu[[k]],Sigma)
    })
    theStat <- CXstat(n,p,K,X)
    tmp <- rep(0,B)
    for(xxx in 1:B){
        tmp[xxx] <- CXstat(n,p,K,myPer(n,p,K,X))
    }
    if(mean(tmp>theStat)<=0.05)
        jjj[myIterator] <- 1
}

# plot(ecdf(jjj))
# 
# TheoryCDF <- function(x){
#     tmp <- uniroot(function(t){gamma(t)-((K-1)/2)},c(1,100))$root
#     exp(-tmp*exp(-x/K))
# }
# 
# curve(TheoryCDF,from=-10,to=30,add=TRUE)

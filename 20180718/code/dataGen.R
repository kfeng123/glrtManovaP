
dataGen <- function(n, p, mu, Sigma) {
    temp <- rnorm(n * p)
    dim(temp) <- c(n, p)
    temp2 <- eigen(Sigma, symmetric = TRUE)
    temp %*% diag(temp2$values ^ (1 / 2)) %*% t(temp2$vectors) + outer(rep(1,n),mu)
}



myPer <- function(n,p,K,X,theOrder){
    tmp <- do.call(rbind,X)
    tmp <- tmp[theOrder,]
    tmp2 <- rep(1:K,times=n)  
    split(data.frame(tmp),tmp2)
}
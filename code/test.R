source('./dataGen.R', echo = TRUE)
source('./statistics.R', echo = TRUE)


doit <- function() {
    pb <- txtProgressBar(style = 3)
    
    jCX <- rep(0, 500)
    jNEW <- rep(0, 500)
    jSC <- rep(0, 500)
    for (myIterator in 1:500) {
        X <- lapply(1:K, function(k) {
            dataGen(n[k], p, mu[[k]], Sigma)
        })
        theCXstat <- CXstat(n, p, K, X)
        
        tmp <- NEWstat(n, p, K, X)
        theNEWstat <- tmp$stat
        Zinv <- tmp$Zinv
        
        tmp <- SCstat(n, p, K, X)
        theSCstat <- tmp$stat
        myGram <- tmp$myGram
        
        tmpNEWstat <- rep(0, B)
        tmpCXstat <- rep(0, B)
        tmpSCstat <- rep(0, B)
        for (xxx in 1:B) {
            theOrder <- sample.int(sum(n))
            thePer <- myPer(n, p, K, X, theOrder)
            tmpCXstat[xxx] <- CXstat(n, p, K, thePer)
            tmpNEWstat[xxx] <-
                NEWstat(n, p, K, thePer, Zinv[theOrder, theOrder])
            tmpSCstat[xxx] <-
                SCstat(n, p, K, thePer, myGram = myGram[theOrder, theOrder])
        }
        if ((sum(tmpCXstat >= theCXstat) + 1) / (B + 1) <= 0.05)
            jCX[myIterator] <- 1
        
        if ((sum(tmpNEWstat >= theNEWstat) + 1) / (B + 1) <= 0.05)
            jNEW[myIterator] <- 1
        
        if ((sum(tmpSCstat >= theSCstat) + 1) / (B + 1) <= 0.05)
            jSC[myIterator] <- 1
        
        setTxtProgressBar(pb, myIterator / 500)
    }
    close(pb)
    return(data.frame(
        CX = mean(jCX),
        SC = mean(jSC),
        NEW = mean(jNEW)
    ))
}

# sample number
K = 3

B = 100

n <- c(10, 10, 10)
p = 80
Sigma <- diag(p)
Sigma[1,1]<-100
mu <- list(rep(0.4, p), rep(0, p), rep(0, p), rep(0, p), rep(0, p))
(doit())

mu <- list(rep(0.2, p), rep(0, p), rep(0, p), rep(0, p), rep(0, p))
(doit())
##### plot(ecdf(jjj))
#
# TheoryCDF <- function(x){
#     tmp <- uniroot(function(t){gamma(t)-((K-1)/2)},c(1,100))$root
#     exp(-tmp*exp(-x/K))
# }
#
# curve(TheoryCDF,from=-10,to=30,add=TRUE)
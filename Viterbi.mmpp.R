"Viterbi.mmpp" <-
function (object, ...) 
{

    #B_i,j is the probability of observing observation j from state i, 
    #here modelled as the density of an expovariate variable
    #
    #A_i,j is the probability of moving from i to j in one time step,
    #Since time is continuous, we use A_i,j(t) - the probability of being in 
    #state j after time t, given that we started in state i, calculated by
    #exponentiation of the infinitesimal generator
 
    library("expm")
 
    a<-0
 
    tau <- object$tau
    lambda <- object$lambda
    delta <- object$delta
    Q <- object$Q
    n <- length(tau)
    m <- length(lambda)
    T1 <- matrix(rep(NA,n*m),m)
    T2 <- matrix(rep(0,n*m),m)
    diffs = rep(NA,n-1)
    
    for (i in 1:(n-1)){
        diffs[i] <- tau[i+1]-tau[i]
    }
    
    for (s in 1:m){
        T1[s,1] <- log(delta[s])+dexp(diffs[1],lambda[s],log=TRUE)
    }
    
    temp <- rep(NA,m)
    
    for (i in 2:(n-1)){
        A <- expm(Q*diffs[i])
        for (j in 1:m){
            for(k in 1:m){
                temp[k] <- T1[k,i-1] + log(A[k,j]) + dexp(diffs[i],lambda[j],log=TRUE)
            }
            T1[j,i] <- max(temp)
            T2[j,i] <- which.max(temp)
        }
    }
    
    print("Calculated Ts")
    
    z <- rep(NA,n-1)
    z[n-1] <- which.max(T1[,n-1])
   
    for (i in (n-1):1){
        z[i-1] <- T2[z[i],i]
    } 
    
    return(z)
}


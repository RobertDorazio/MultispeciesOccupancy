library(mvtnorm)
library(msm)
library(MASS)
##library(clusterGeneration)
library(tmvtnorm)


### Define function for fitting a multivariate probit model using a modification of the PX-SGS-DA algorithm
### of Xiao Zhang (2020; Statistics in Medicine 39: 3637-3652)

MvProbitModel_MCMC = function(Z, X, niter=1, niterInterval=NA, quiet=FALSE) {


    ### Define functions used in Metropolis-Hastings sampling

    ## ... no functions
    


    ### Define function for drawing a random value from a truncated multivariate normal distribution using conditional distributions

    ## sampleTruncatedMultivariateNormalConditionally = function(z, x, mu, R) {

    ##     xnew = x
    ##     m = length(mu)

    ##     for (j in 1:m) {

    ##         Rvec = R[j, ][-j]
    ##         ind = (1:m)[-j]
    ##         Rsub = R[ind, ind]
    ##         RsubInv = chol2inv(chol(Rsub))

    ##         xnew.mean = mu[j] + t(Rvec) %*% RsubInv %*% (xnew[-j] - mu[-j])
    ##         xnew.var = R[j,j] - t(Rvec) %*% RsubInv %*% Rvec
    ##         xnew[j] = ifelse(z[j]==1,
    ##                          rtnorm(1, mean=xnew.mean, sd=sqrt(xnew.var), lower=0),
    ##                          rtnorm(1, mean=xnew.mean, sd=sqrt(xnew.var), upper=0))
    ##     }
    ##     return(xnew)
    ## }


    ## sampleTruncatedMultivariateNormalConditionally = function(z, x, mu, Sigma) {

    ##     xlower = ifelse(z==1, 0, -Inf)
    ##     xupper = ifelse(z==0, 0, Inf)
    ##     xnew = rtmvnorm(1, mean=mu, sigma=Sigma, lower=xlower, upper=xupper, algorithm='gibbs', start.value=x, burn.in.samples=100)
    ##     return(xnew)
    ## }



    ### Initialize Markov chain

    n = nrow(Z)
    m = ncol(Z)
    p = ncol(X)
    
    XprimeX = t(X) %*% X
    invOfXprimeX = chol2inv(chol(XprimeX))


    ## ... assign values of hyperparameters

    ## ... hierarchical prior of Sigma (Huang and Wand 2013)
    
    nu.Sigma = 2  # corrresponds to a uniform(-1,1) marginal prior for each correlation parameter of Sigma
    scale.Sigma = rep(10, m)  # arbitrarily high values induce a noninformative Half-t prior on each diagonal element of Sigma


    ## .... conjugate prior for B

    scale.B = 10  # arbitrarily high value induces a noninformative prior
    SigmaInv.B = (1/scale.B)^2 * diag(rep(1,p))
    SigmaConditional.B = chol2inv(chol(SigmaInv.B + XprimeX))

    
        
    ## ... initialize each column of B by fitting a probit regression model to each column of Z
    
    B = matrix(nrow=p, ncol=m)
    
    for (j in 1:m) {
        zvec = Z[,j]
        fit = glm(zvec ~ X - 1, family=binomial(link='probit'))
        B[,j] = fit$coefficients
    }

 
    
    ## ... initialize elements of R
    R = diag(rep(1,m))
    
    ## ... initialize elements of Dvec (= diagonal elements of D)
    ## ... (each element of Dvec must be positive!)
    DvecGuess = rep(1,m)
    Dvec = DvecGuess

    ## ... compute Sigma and SigmaInv using initialized values of R and Dvec
    D = diag(Dvec)
    Sigma = D %*% R %*% D
    SigmaInv = chol2inv(chol(Sigma))


    
    ## ... initialize elements of W conditional on B, Dvec, and Z while counterfactually assuming Sigma=diag(Dvec^2)
    W = matrix(nrow=n, ncol=m)
    mu = X %*% B
    
    for (i in 1:n) {
        for (j in 1:m) {
            W[i,j] = ifelse(Z[i,j]==1,
                            rtnorm(1, mean=mu[i,j], sd=D[j,j], lower=0),
                            rtnorm(1, mean=mu[i,j], sd=D[j,j], upper=0))
        }                   
    }
    
    
    
    ### Run the MCMC algorithm

    nparam.R = m * (m-1) / 2  # no of R params
    mc.vecR = matrix(nrow=niter, ncol=nparam.R)  # this matrix holds the Markov chain for R
    j.index = integer(0)
    i.index = integer(0)
    for (j in 1:(m-1)) {
        for (i in (j+1):m) {
            j.index = c(j.index, j)
            i.index = c(i.index, i)
        }
    }
    colnames(mc.vecR) = paste('R[', i.index, ',', j.index, ']', sep='')

    
    nparam.B = p*m  # no. of B params = no. of Beta params
    mc.vecBeta = matrix(nrow=niter, ncol=nparam.B)  # this matrix holds the Markov chain for Beta
    i.index = rep(1:p, m)
    j.index = rep(1:m, each=p)
    colnames(mc.vecBeta) = paste('Beta[', i.index, ',', j.index, ']', sep='')


    start.time = Sys.time()

    for (iter in 1:niter) {

        ## ... draw W | Sigma, B, Z

        mu = X %*% B
        for (i in 1:n) {
            
            ## Wvec = sampleTruncatedMultivariateNormalConditionally(Z[i,], W[i,], mu[i,], Sigma)
            wLower = ifelse(Z[i,]==1, 0, -Inf)
            wUpper = ifelse(Z[i,]==0, 0, Inf)
            ## Wvec = rtmvnorm(1, mean=mu[i,], sigma=Sigma, lower=wLower, upper=wUpper, algorithm='rejection')
            Wvec = rtmvnorm(1, mean=mu[i,], sigma=Sigma, lower=wLower, upper=wUpper, algorithm='gibbs', start.value=W[i,], burn.in.samples=1000)
            W[i,] = Wvec
        }

        

        ## ... draw B | W, Sigma
        
        B.mle = invOfXprimeX %*% t(X) %*% W

        vecB = rmvnorm(1, mean=as.vector(B.mle), sigma=kronecker(Sigma, invOfXprimeX))
        
        B = matrix(vecB, nrow=p, ncol=m)


        
        ## ... draw Sigma | B, W

        Eps = W  -  X %*% B
        crossProdMat = crossprod(Eps, Eps)  # equivalent to  t(Eps) %*% Eps  (but faster!)

        b = rgamma(m, shape=(nu.Sigma+m)/2, rate=(1/scale.Sigma^2) + nu.Sigma*diag(SigmaInv) )

        SigmaInv = rWishart(1, df=nu.Sigma+m-1 + n, Sigma=chol2inv(chol(2*nu.Sigma*diag(b) + crossProdMat)))[,, 1]
        Sigma = chol2inv(chol(SigmaInv))

        

        ## ... compute R and Beta from the drawn values of Sigma and B
        Dvec = sqrt(diag(Sigma))
        Dinv = diag(1/Dvec)
        R = Dinv %*% Sigma %*% Dinv
        Beta = B %*% Dinv
        
       

        ## ... save iteration's results to markov chain

        indLower = lower.tri(R)
        vecRlower = as.vector(R[indLower])
        mc.vecR[iter, ] = vecRlower
        mc.vecBeta[iter, ] = as.vector(Beta)


        ## ... optionally report progress of MCMC algorithm

        if (!quiet & (!is.na(niterInterval)) & iter==round(iter/niterInterval)*niterInterval) {
            current.time = Sys.time()
            elapsed.time = difftime(current.time, start.time, units='mins')
            ndigits = ifelse(elapsed.time>1, 1,2)
            cat(paste(paste('... completed iteration #', iter, 'after', round(elapsed.time,ndigits), 'minutes \n', sep=' ')))
        }
       
        
    }  # end of loop for iter


    return(list(mc.vecR=mc.vecR, mc.vecBeta=mc.vecBeta))
}

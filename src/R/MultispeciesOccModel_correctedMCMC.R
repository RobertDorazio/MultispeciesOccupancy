library(mvtnorm)
library(msm)
library(MASS)
library(tmvtnorm)



### Define function for fitting a multispecies occupancy model
### (wherein a multivariate probit regression model is used to specify occupancy states of each species)

MultispeciesOccModel_MCMC = function(Y, X, V, niter=500, niterInterval=floor(500/4), quiet=FALSE) {


    ### Define functions used in Metropolis-Hastings sampling

    ## ... log of unnormalized pdf for conditional posterior of Alpha

    logDensity.alpha = function(alpha, yvec, vmat, mu.alpha, sigma.alpha) {
        
        logPrior = sum(dnorm(alpha, mean=mu.alpha, sd=sigma.alpha, log=TRUE))
        pvec = as.vector(pnorm(vmat %*% alpha))
        logDensity = sum(yvec*log(pvec) + (1-yvec)*log(1-pvec))
        return(logDensity + logPrior)
    }
    
   ## ... dMvn() is used internally in the bayesm R package.   It is much faster than dmvnorm() in the mvtnorm R package.
    dMvn <- function(X,mu,Sigma) {
        k <- ncol(X)
        rooti <- backsolve(chol(Sigma),diag(k))
        quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
        return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))
    }

    


    

    ### Initialize Markov chain

    n = dim(Y)[1] # no. sites
    m = dim(Y)[3] # no. species
    p = ncol(X)   # no. regressors of occurrence probability (including intercept)
    q = dim(V)[3] # no. regressors of detection probability  (including intercept)
    K = rowSums(!is.na(Y[,,1]))  # vector of no. surveys per site

    ## ... Y is an array of detections               (dimensions = n, max(K), m)
    ## ... X is a matrix of regressors of occurrence (dimensions = n, p)
    ## ... V is an array of regressors of detection  (dimensions = n, max(K), q)


    
    ## ... compute quantities that are not affected by MCMC algorithm
    XprimeX = t(X) %*% X
    invOfXprimeX = chol2inv(chol(XprimeX))

    Ysum = apply(Y, c(1,3), sum, na.rm=TRUE)

    
    
    ## ... check data for adequacy
    yvec = colSums(Ysum > 0)
    if (any(yvec==0)) {
        cat('Data contain ', sum(yvec==0), ' species with no detections. \n')
        return(NULL)
    }  
    if (any(yvec==n)) {
        cat('Data contain ', sum(yvec==n), ' species that were detected at every location. \n')
        return(NULL)
    }  


    
    ## ... assign values of hyperparameters

    ## ... for hierarchical prior of Sigma (Huang and Wand 2013)
    
    nu.Sigma = 2  # corrresponds to a uniform(-1,1) marginal prior for each correlation parameter of Sigma
    scale.Sigma = rep(10, m)  # arbitrarily high values induce a noninformative Half-t prior on each diagonal element of Sigma


    ## ... for prior of Alpha
    
    mu.alpha = 0
    sigma.alpha = 1



    
    ## ... initialize occupancy values in Z by aggregating detections in Y among surveys

    Z = matrix(as.integer(Ysum>0), nrow=n, ncol=m)
    

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
    Dvec = rep(1,m)

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

    
    
    ## ... initialize jth column of Alpha by fitting a probit regression model that uses
    ##     data from sites where the jth species was detected at least once
    
    Alpha = matrix(nrow=q, ncol=m)

    for (j in 1:m) {
        ## ... compute vector of detections and matrix of regressors for jth species; then fit probit regression model
        yvec = integer(0)
        vmat = matrix(nrow=0, ncol=q)
        for (i in 1:n) {
            if (Ysum[i,j] > 0) {
                yvec = c(yvec, Y[i, 1:K[i], j])
                vmat = rbind(vmat, matrix(V[i, 1:K[i], ], nrow=K[i], ncol=q) )
            }
        }
        fit = glm(yvec ~ vmat - 1, family=binomial(link='probit'))
        Alpha[,j] = fit$coefficients
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

    
    nparam.B = p*m  # no. of B params (= no. of Beta params)
    mc.vecBeta = matrix(nrow=niter, ncol=nparam.B)  # this matrix holds the Markov chain for Beta
    i.index = rep(1:p, m)
    j.index = rep(1:m, each=p)
    colnames(mc.vecBeta) = paste('Beta[', i.index, ',', j.index, ']', sep='')

    
    nparam.Alpha = q*m  # no. of Alpha params
    mc.vecAlpha = matrix(nrow=niter, ncol=nparam.Alpha)  # this matrix holds the Markov chain for Alpha
    i.index = rep(1:q, m)
    j.index = rep(1:m, each=q)
    colnames(mc.vecAlpha) = paste('Alpha[', i.index, ',', j.index, ']', sep='')



    start.time = Sys.time()

    for (iter in 1:niter) {

        ## ... draw W | Sigma, B, Z

        mu = X %*% B
        for (i in 1:n) {
            wLower = ifelse(Z[i,]==1, 0, -Inf)
            wUpper = ifelse(Z[i,]==0, 0, Inf)

            ## ... check whether start value for W[i,] falls within the support region given current values of wLower and wUpper
            wStart = W[i, ]
            wInside = (W[i,] >= wLower) & (W[i,] <= wUpper)
            if (!all(wInside)) {
                for (j in 1:m) {
                    if (!wInside[j]) wStart[j] = (-1) * wStart[j]
                }
            }
            Wvec = rtmvnorm(1, mean=mu[i,], sigma=Sigma, lower=wLower, upper=wUpper, algorithm='gibbs', start.value=wStart,
                            burn.in.samples=1000)
            
            if (all(!is.nan(Wvec))) {  # check if sampling of truncated multivariate normal distribution failed)
                W[i,] = Wvec
            }
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



        ## ... draw Alpha | Y, Z
        ##     (draws are only needed for sites where a species is present based on its occupancy state Z)
        ##     (draw each column of Alpha independently

        for (j in 1:m) {
            ## ... To obtain proposal for Alpha_j, fit a probit regression model using data from sites
            ## ... where the jth species is present; then use a multivariate normal distribution as a Metropolis-Hastings
            ## ... proposal with mean and variance equal to the MLE and covariance matrix, respectively.

            ## ... compute vector of detections and matrix of regressors for jth species; then fit probit regression model
            yvec = integer(0)
            vmat = matrix(nrow=0, ncol=q)
            for (i in 1:n) {
                if (Z[i,j] > 0) {
                    yvec = c(yvec, Y[i, 1:K[i], j])
                    vmat = rbind(vmat, matrix(V[i, 1:K[i], ], nrow=K[i], ncol=q) )
                }
            }
            fit = glm(yvec ~ vmat - 1, family=binomial(link='probit'))

            if (fit$converged & !fit$boundary) {
  
                ## ... draw candidate using multivariate normal distribution as proposal
                meanOfProposal = fit$coefficients
                SigmaOfProposal = vcov(fit)
                alpha.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))

                logL.alpha.cand = logDensity.alpha(alpha.cand, yvec, vmat, mu.alpha, sigma.alpha)
                logL.alpha = logDensity.alpha(Alpha[,j], yvec, vmat, mu.alpha, sigma.alpha)
                ## logQ.alpha.cand = dmvnorm(alpha.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
                ## logQ.alpha = dmvnorm(Alpha[,j], mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
                logQ.alpha.cand = log( dMvn(matrix(alpha.cand,nrow=1), mu=meanOfProposal, Sigma=SigmaOfProposal) )
                logQ.alpha = log( dMvn(matrix(Alpha[,j],nrow=1), mu=meanOfProposal, Sigma=SigmaOfProposal) )
                logR = logL.alpha.cand - logL.alpha - logQ.alpha.cand + logQ.alpha
                if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
                    Alpha[,j] = alpha.cand
                }
            }
        }

        

        
        ## ... draw Z | Y, W, B, Sigma, Alpha
        ##     (draws are only needed if a species was not detected at least once i.e., Ysum[i,j] = 0)
        
        ## ... compute mean and variance for the normal distribution of  W[i,j] | W[i, -j]
        ## ... compute probNotDetected in K[i] surveys
        ## ... then use Bayes' rule to calculate Pr( W[i,j] > 0 | Ysum[i,j]=0, W[i, -j] )
        
        for (i in 1:n) {
            for (j in 1:m) {
                if (Ysum[i,j]==0) {
                    
                    vmat = matrix(V[i, 1:K[i], ], nrow=K[i], ncol=q)
                    pvec = pnorm(vmat %*% Alpha[,j])
                    ## probNotDetected = prod(1-pvec)
                    probNotDetected = exp(sum(log(1.0-pvec)))

                    ## .. calculate mean and variance for conditional distribution of  W[i, j] | W[i, -j]
                    
                    mu = t(B) %*% X[i, ]
                    ind = (1:m)[-j]
                    w_vec = W[i, ind]
                    Sigma_vec = Sigma[j, ][-j]     # off-diagonal elements of Sigma in row j

                    Sigma_sub = Sigma[ind, ind]    # principal submatrix of Sigma
                    Sigma_subInverse = chol2inv(chol(Sigma_sub))

                    wMean = mu[j] + t(Sigma_vec) %*% Sigma_subInverse %*% (w_vec - mu[ind])
                    wVar  = Sigma[j,j] - t(Sigma_vec) %*% Sigma_subInverse %*% Sigma_vec
                    
                    probPositiveW = pnorm(wMean/sqrt(wVar))   # Pr(W[i, j] > 0  | W[i, -j] )
                    if (probPositiveW > 0.9999) {
                        Z[i,j] = 1  # analytically, as psi[i,j] -> 1,  Z[i,j]  ~ Bernoulli(1)
                    }
                    else {
                        numerator = probNotDetected * probPositiveW
                        probZ = numerator / (numerator + 1 - probPositiveW)
                        Z[i,j] = rbinom(1, size=1, prob=probZ)
                    }
                }
            }
        }
        
       

        ## ... save iteration's results to markov chain

        indLower = lower.tri(R)
        vecRlower = as.vector(R[indLower])
        mc.vecR[iter, ] = vecRlower
        mc.vecBeta[iter, ] = as.vector(Beta)
        mc.vecAlpha[iter, ] = as.vector(Alpha)


        ## ... optionally report progress of MCMC algorithm

        if (!quiet &  (!is.na(niterInterval) & iter==round(iter/niterInterval)*niterInterval)) {
            current.time = Sys.time()
            elapsed.time = difftime(current.time, start.time, units='mins')
            ndigits = ifelse(elapsed.time>1, 1,2)
            cat(paste(paste('... completed iteration #', iter, 'after', round(elapsed.time,ndigits), 'minutes \n', sep=' ')))
        }
       
        
    }  # end of loop for iter


    return(list(mc.vecR=mc.vecR, mc.vecBeta=mc.vecBeta, mc.vecAlpha=mc.vecAlpha))
}

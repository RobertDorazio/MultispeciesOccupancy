library(mvtnorm)
library(msm)
library(MASS)
library(rjags)


### Define function for fitting a multispecies occupancy model with latent variables (factors) using JAGS 

MultispeciesOccModel_LV_JAGS = function(nfactors, Y, X, V, niter=500, niterToAdapt=500, quiet=FALSE) {

    r = nfactors  # no. latent factors used to approximate the correlation matrix 
    n = dim(Y)[1] # no. sites
    m = dim(Y)[3] # no. species
    p = ncol(X)   # no. regressors of occurrence probability (including intercept)
    q = dim(V)[3] # no. regressors of detection probability  (including intercept)
    K = rowSums(!is.na(Y[,,1]))  # vector of no. surveys per site
    
    ## ... Y is an array of detections               (dimensions = n, max(K), m)
    ## ... X is a matrix of regressors of occurrence (dimensions = n, p)
    ## ... V is an array of regressors of detection  (dimensions = n, max(K), q)

    
    ## ... create arguments for JAGS

    dataList = list(r=r, n=n, m=m, p=p, q=q, K=K, Y=Y, X=X, V=V,
                    mu.Beta=0, tau.Beta=1, mu.Alpha=0, tau.Alpha=1)
    
    params = c('Beta', 'Omega', 'Alpha')


    
    inits = function() {
    
    ## ... initialize each column of Beta by fitting a probit regression model to each column of Z

    Ysum = apply(Y, c(1,3), sum, na.rm=TRUE)
    Z = matrix(as.integer(Ysum>0), nrow=n, ncol=m)
    
    Beta = matrix(nrow=p, ncol=m)
    
    for (j in 1:m) {
        zvec = Z[,j]
        fit = glm(zvec ~ X - 1, family=binomial(link='probit'))
        Beta[,j] = fit$coefficients
    }

    ## ... initialize elements of U conditional on Beta and Z while counterfactually assuming Sigma=I
    U = matrix(nrow=n, ncol=m)
    mu = X %*% Beta
    
    for (i in 1:n) {
        for (j in 1:m) {
            U[i,j] = ifelse(Z[i,j]==1,
                            rtnorm(1, mean=mu[i,j], sd=1, lower=0),
                            rtnorm(1, mean=mu[i,j], sd=1, upper=0))
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
    
    list(Beta=Beta, U=U, Alpha=Alpha)
}


modelFilename = "LVmodel.txt"

  cat("
model {


# prior for Beta

for (l in 1:p) {
for (j in 1:m) {
Beta[l,j] ~ dnorm(mu.Beta, tau.Beta)
}
}


# prior for Alpha

for (l in 1:q) {
for (j in 1:m) {
Alpha[l,j] ~ dnorm(mu.Alpha, tau.Alpha)
}
}


# constraint on elements of Lambda that cannot be estimated

for (j in 1:(r-1)) {
for (l in (j+1):r) {
  Lambda[j,l] <- 0
}
}


# prior for estimable elements of Lambda

for (j in 1:r) {

  Lambda[j,j] ~ dunif(0,1)

for (l in (j+1):m) {

  Lambda[l,j] ~ dunif(-1,1)

}
}


# define transpose of Lambda
for (j in 1:m) {
for (l in 1:r) {
  tLambda[l,j] <- Lambda[j,l]
}
}



# constraint on psi to ensure Omega[j,j]=1
for (j in 1:m) {
  psi[j] <- 1 - sum(Lambda[j,1:r]^2)
}


# specify distribution of latent factors
for (i in 1:n) {
for (l in 1:r) {
  f[i,l] ~ dnorm(0,1)
}
}



# likelihood

for (i in 1:n) {

for (j in 1:m) {


# occupancy model

mu[i,j] <- X[i, ] %*% Beta[, j] + f[i, ] %*% tLambda[, j]

U[i,j] ~ dnorm(mu[i,j], 1/psi[j])

Z[i,j] ~ dinterval(U[i,j], 0)  # censored distribution of latent U[i,j]


for (k in 1:K[i]) {

# detection model

probit(probdetect[i,k,j]) <- V[i,k, ] %*% Alpha[, j]

prob.y[i,k,j] <- Z[i,j] * probdetect[i,k,j]

Y[i,k,j] ~ dbern(prob.y[i,k,j])

} # end of k loop


} # end of j loop

} # end of i loop



# calculation of Omega (= approximation of correlation matrix)

for (j in 1:m) {
for (k in 1:m) {
  PsiDiagonalMat[j,k] <- ifelse(j==k,  psi[j], 0)
}
}

Omega <- Lambda %*% tLambda + PsiDiagonalMat

}
", fill=TRUE, file=modelFilename)



## ... fit model to data using JAGS (called using functions in rjags package)

  start.time = Sys.time()
  jmod = jags.model(file=modelFilename, data=dataList, inits=inits, n.chains=1, n.adapt=niterToAdapt, quiet=quiet)
  mc = coda.samples(jmod, params, n.iter=niter, thin=1, quiet=quiet)

  end.time = Sys.time()
  elapsed.time = difftime(end.time, start.time, units='mins')
  cat(paste(paste('Markov chain computed in ', elapsed.time, sep=''), ' minutes\n', sep=''))


## ... return Markov chains of Alpha, Beta, and R as separate matrices

    mchain = as.matrix(mc, chain=FALSE)

    ## ... extract R from Markov chain

    nparam.R = m * (m-1) / 2  # no of R params
    j.index = integer(0)
    i.index = integer(0)
    for (j in 1:(m-1)) {
        for (i in (j+1):m) {
            j.index = c(j.index, j)
            i.index = c(i.index, i)
        }
    }
    OmegaNames = paste('Omega[', i.index, ',', j.index, ']', sep='')
    mc.vecR = mchain[, colnames(mchain) %in% OmegaNames]
    Rnames = paste('R[', i.index, ',', j.index, ']', sep='')
    colnames(mc.vecR) = Rnames

    
    ## ... extract Beta from Markov chain

    nparam.B = p*m  # no. of B params (= no. of Beta params)
    i.index = rep(1:p, m)
    j.index = rep(1:m, each=p)
    BetaNames = paste('Beta[', i.index, ',', j.index, ']', sep='')
    mc.vecBeta = mchain[, colnames(mchain) %in% BetaNames]

    
    ## ... extract Alpha from Markov chain
    
    nparam.Alpha = q*m  # no. of Alpha params
    i.index = rep(1:q, m)
    j.index = rep(1:m, each=q)
    AlphaNames = paste('Alpha[', i.index, ',', j.index, ']', sep='')
    mc.vecAlpha = mchain[, colnames(mchain) %in% AlphaNames]
   
    
    return(list(mc.vecR=mc.vecR, mc.vecBeta=mc.vecBeta, mc.vecAlpha=mc.vecAlpha))
}

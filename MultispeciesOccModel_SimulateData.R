library(MASS)
library(clusterGeneration)
library(mvtnorm)


### Define function for simulating a dataset using:
### ... (1) Model of multispecies occupancy states based on a multivariate probit regression model
### ... (2) Model of detections of species conditional on their occupancy states
### ... Model parameters are drawn randomly from normal distributions (see code).
### ... Model regressors are drawn randomly from uniform(-1,1) distributions and then centered and scaled.  

MultiSpeciesOccModel_SimulateData = function(n, m, p, q, Kmin, Kmax, seed=NA) {

## n = no. sites
## m = no. species
## p = no. occupancy predictors (including intercept)
## q = no. detection predictors (including intercept)
## Kmin = minimum no. surveys per site
## Kmax = maximum no. surveys per site
## seed = integer optionally use to seed random nummber generator

    if (!is.na(seed)) set.seed(seed)

    
### Simulate multispecies occupancy states using multivariate probit regression model

R.param = rcorrmatrix(m, alphad=5)  # a random, positive-definite correlation matrix

Beta.param = matrix(nrow=p, ncol=m)
Beta.param[1,] = qnorm(runif(m, 0.20, 0.80)) # defined on probit scale
if (p>1) {
    for (i in 2:p) {
        Beta.param[i,] = rnorm(m, 0, sd=0.5) # defined on probit scale
    }
}

X = matrix(1, nrow=n, ncol=1)
if (p>1) {
    Xcovs = matrix(runif(n*(p-1), min=-1, max=1), nrow=n)
    X = cbind(X, scale(Xcovs))
}


mu = X %*% Beta.param
umat = matrix(nrow=n, ncol=m)
for (i in 1:n) {
    umat[i,] = rmvnorm(1, mean=mu[i,], sigma=R.param)
}

Z = matrix(as.integer(umat > 0), nrow=nrow(umat), ncol=ncol(umat))


### Simulate detections of species conditional on their occupancy states

K = sample(Kmin:Kmax, size=n, replace=TRUE)  # vector of no. surveys per site

Alpha.param = matrix(nrow=q, ncol=m)
Alpha.param[1,] = qnorm(runif(m, 0.1, 0.7))   # defined on probit scale
if (q>1) {
    for (i in 2:q) {
        Alpha.param[i,] = rnorm(m, 0, sd=0.5) # defined on probit scale
    }
}

V = array(dim=c(n, max(K), q))

for (i in 1:n) {
    for (k in 1:K[i]) {
        V[i,k,1] =1
        if (q > 1) {
            V[i,k, 2:q] = runif(q-1, min=-1, max=1)
        }
    }
}
## ... center and scale each predictor to have zero mean and unit variance
if (q > 1) {
    Vmean = rep(NA, q)
    Vsd = rep(NA, q)
    for (l in 2:q) {
        tempvec = numeric(0)
        for (i in 1:n) {
            for (k in 1:K[i]) {
                tempvec = c(tempvec, V[i,k,l])
            }
        }
        Vmean[l] = mean(tempvec)
        Vsd[l] = sd(tempvec)
    }
    ## ... use mean and sd to center and scale predictors
    for (l in 2:q) {
        for (i in 1:n) {
            for (k in 1:K[i]) {
                V[i,k,l] = (V[i,k,l] - Vmean[l]) / Vsd[l]
            }
        }
    }
}


    
Y = array(dim=c(n, max(K), m))
for (j in 1:m) {
    for (i in 1:n) {
        for (k in 1:K[i]) {
            prob = pnorm(V[i,k, ] %*% Alpha.param[,j])
            Y[i,k,j] = rbinom(1, size=1, prob=Z[i,j]*prob)
        }
    }
}


return(list(Y=Y, X=X, V=V, Z=Z, Beta=Beta.param, Alpha=Alpha.param, R=R.param))
}

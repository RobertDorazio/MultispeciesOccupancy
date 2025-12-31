source('MultispeciesOccModel_SimulateData.R')
source('MultispeciesOccModel_correctedMCMC.R')
source('MultispeciesOccModel_LV_JAGS.R')
source('MvProbitModel_MCMC.R')
source('EstimatePosteriorStats.R')



### Assign constants for simulating data

n = 100    # no. sites
m =   5    # no. species

p = 3      # no. occupancy predictors (including intercept)
q = 2      # no. detection predictors (including intercept)
Kmin =  3  # minimum no. surveys per site
Kmax =  5  # maximum no. surveys per site



### Simulate data

dlist = MultiSpeciesOccModel_SimulateData(n, m, p, q, Kmin, Kmax, seed=NA)




### Fit MP occupancy model

CR = '\n'
cat(CR, CR, 'Fitting MP occupancy model ... ', CR)

fit = MultispeciesOccModel_MCMC(Y=dlist$Y, X=dlist$X, V=dlist$V, niter=25000, niterInterval=1000, quiet=FALSE)


## ... estimate posterior summary statistics for species occurrence parameters (Beta)

postStats = EstimatePosteriorStats(mc=fit$mc.vecBeta, burnin=15000)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))


## ... estimate posterior summary statistics for correlations in occurrence between species  (R)

postStats = EstimatePosteriorStats(mc=fit$mc.vecR, burnin=15000)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))


## ... estimate posterior summary statistics for species detection parameters (Alpha)

postStats = EstimatePosteriorStats(mc=fit$mc.vecAlpha, burnin=15000)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))




### Fit LF occupancy model

CR = '\n'
cat(CR, CR, 'Fitting LF occupancy model ... ', CR)

fit = MultispeciesOccModel_LV_JAGS(nfactors=floor(m/2), Y=dlist$Y, X=dlist$X, V=dlist$V, niter=75000, niterToAdapt=2000, quiet=FALSE)


## ... estimate posterior summary statistics for species occurrence parameters (Beta)

postStats = EstimatePosteriorStats(mc=fit$mc.vecBeta, burnin=15000, thin=10)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))


## ... estimate posterior summary statistics for correlations in occurrence between species  (R)

postStats = EstimatePosteriorStats(mc=fit$mc.vecR, burnin=15000, thin=10)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))


## ... estimate posterior summary statistics for species detection parameters (Alpha)

postStats = EstimatePosteriorStats(mc=fit$mc.vecAlpha, burnin=15000, thin=10)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))




### Fit Multivariate Probit Regression model

CR = '\n'
cat(CR, CR, 'Fitting multivariate probit regression model ... ', CR)

fit = MvProbitModel_MCMC(Z=dlist$Z, X=dlist$X, niter=25000, niterInterval=1000, quiet=FALSE)


## ... estimate posterior summary statistics for species occurrence parameters (Beta)

postStats = EstimatePosteriorStats(mc=fit$mc.vecBeta, burnin=15000)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))


## ... estimate posterior summary statistics for correlations in occurrence between species  (R)

postStats = EstimatePosteriorStats(mc=fit$mc.vecR, burnin=15000)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))

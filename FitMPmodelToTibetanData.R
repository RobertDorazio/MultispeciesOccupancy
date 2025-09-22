
source('MultispeciesOccModel_MCMC.R')
source('EstimatePosteriorStats.R')



### Read Tibetan data from rds file

dlist = readRDS('realdata.rds')

dlist$Y = dlist$Y[,, -c(4,8)]  #  omit detections of Brown bear and Pallas cat (too few detections for analysis)



### Fit MP model using MCMC algorithm

fit = MultispeciesOccModel_MCMC(Y=dlist$Y, X=dlist$X, V=dlist$V, niter=100000, niterInterval=10000, quiet=FALSE)


## ... estimate posterior summary statistics for correlations in occurrence between species  (R)

postStats = EstimatePosteriorStats(mc=fit$mc.vecR, burnin=50000, prob.quantiles=c(0.5, 0.025, 0.05, 0.95, 0.975))

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))



## ... estimate posterior summary statistics for species occurrence parameters (Beta)

postStats = EstimatePosteriorStats(mc=fit$mc.vecBeta, burnin=50000)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))



## ... estimate posterior summary statistics for species detection parameters (Alpha)

postStats = EstimatePosteriorStats(mc=fit$mc.vecAlpha, burnin=50000)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))

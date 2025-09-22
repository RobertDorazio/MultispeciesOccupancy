
source('MultispeciesOccModel_LV_JAGS.R')
source('EstimatePosteriorStats.R')



### Read Tibetan data from rds file

dlist = readRDS('realdata.rds')

dlist$Y = dlist$Y[,, -c(4,8)]  #  omit detections of Brown bear and Pallas cat (too few detections for analysis)





### Fit latent-variable model using JAGS algorithm

nfactors = 7  #  no. latent factors cannot exceed no. species 

fit = MultispeciesOccModel_LV_JAGS(nfactors, Y=dlist$Y, X=dlist$X, V=dlist$V, niter=300000, niterToAdapt=2000, quiet=FALSE)



## ... estimate posterior summary statistics for correlations in occurrence between species  (R)


postStats = EstimatePosteriorStats(mc=fit$mc.vecR, burnin=50000, thin=10, prob.quantiles=c(0.5, 0.025, 0.05, 0.95, 0.975))

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))


## ... estimate posterior summary statistics for species occurrence parameters (Beta)


postStats = EstimatePosteriorStats(mc=fit$mc.vecBeta, burnin=50000, thin=10)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))


## ... estimate posterior summary statistics for species detection parameters (Alpha)


postStats = EstimatePosteriorStats(mc=fit$mc.vecAlpha, burnin=50000, thin=10)

CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(postStats$estimate, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(postStats$MCerror,4))

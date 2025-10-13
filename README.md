
# Multi-species Occupancy Models with Correlated Species Occurrences

This repository contains R code for

1. simulating data based on the assumptions of a multivariate probit occupancy model with correlated species occurrences

2. fitting this model and its latent-factor approximation to the simulated data

3. fitting a multivariate probit regression model to the matrix of simulated species occurrences

4. fitting the two multi-species occupancy models to the Tibetan Plateau data


The following publication contains a complete description of these models and the MCMC algorithms used to fit them.


[Dorazio, Li, Xiao, and Kachel (2025)  An evaluation of multi-species occupancy models with correlated species occurrences.   Methods in Ecology and Evolution (in press)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70168)



## List of files

- *MultispeciesOccModel_SimulateData.R* - simulates data based on the assumptions of a multivariate probit (MP) occupancy model

- *MultispeciesOccModel_MCMC.R* - fits a multivariate probit (MP) occupancy model to data using an efficient MCMC algorithm

- *MultispeciesOccModel_LV_JAGS.R* - fits a latent factor (LF) occupancy model to data using JAGS

- *MvProbitModel_MCMC.R* - fits a multivariate probit regression model (MvP) to latent occurrences using an efficient MCMC algorithm

- *SimulationStudy-Driver.R* - calls functions to simulate a data set and then fit each model to this data set


- *realdata.rds*  - RDS file containing the camera-trap survey data of mammalian species on the Tibetan Plateau

- *FitMPmodelToTibetanData.R* - driver for fitting the MP model to the Tibetan Plateau data

- *FitLVmodelToTibetanData.R*  - driver for fitting the LF model to the Tibetan Plateau data


- *EstimatePosteriorStats.R* -- estimates posterior summary statistics from Markov chains



## Required R packages

`mvtnorm, tmvtnorm, msm, MASS, clusterGeneration, mcmcse, coda, rjags`


## Correspondence

Robert Dorazio (<rmdorazio@gmail.com>) and Lingyun Xiao (<Lingyun.Xiao@xjtlu.edu.cn>)


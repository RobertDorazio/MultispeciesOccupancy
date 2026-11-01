### Include packages and source code


### Include Julia code

path = pwd()
path = replace(path, "tests" => "src/Julia")

cd(path) do          
    include("MultispeciesOccModel.jl")
end

using .MultispeciesOccModel

using Dates


### Simulate a dataset

n = 400   # no. sites
m = 5      # no. species
Kmin = 3   # minimum no. surveys per site during each time period
Kmax = 5   # maximum no. surveys per site during each time period

qcovs = 2  # no. regressors of detection (including intercept)
pcovs = 3  # no. regressors of occurrence (including intercept)


simList = simulateData(n, m, Kmin, Kmax, qcovs, pcovs)



### Fit model using MCMC algorithm

dt_start = now()

niter = 30000
niterInterval = 1000

mc = runMCMC(simList.data,  niter, niterInterval, false)

dt_finish = now()

println("\n", "Total elapsed hours = ", round((dt_finish - dt_start) / Dates.Hour(1), digits=2) )



###  Analyze Markov chain:  write to csv file, plot, or estimate summaries of posterior distribution

using MCMCChains, StatsPlots, DataFrames, CSV

## ... write chain to file

df = DataFrame(mc.vecBeta, mc.vecBeta_name)
CSV.write("mc_vecBeta.csv", df)

df = DataFrame(mc.vecAlpha, mc.vecAlpha_name)
CSV.write("mc_vecAlpha.csv", df)

df = DataFrame(mc.vecR, mc.vecR_name)
CSV.write("mc_vecR.csv", df)



## ... plot chain's trace or autocorrelations

chain = Chains(mc.vecBeta, mc.vecBeta_name)
chain = Chains(mc.vecAlpha, mc.vecAlpha_name)
chain = Chains(mc.vecR, mc.vecR_name)

traceplot(chain)
autocorplot(chain)


## ... estimate summaries of posterior distribution

burnin = 15000
posteriorSeq = (burnin+1):niter

chain = Chains(mc.vecBeta[posteriorSeq, :], mc.vecBeta_name)
chain = Chains(mc.vecAlpha[posteriorSeq, :], mc.vecAlpha_name)
chain = Chains(mc.vecR[posteriorSeq, :], mc.vecR_name)


describe(chain)  # print summaries

## ... save summaries as dataframes

df = DataFrame(summarystats(chain))
posterior_stats = transform(df, names(df, Float64) .=> (x -> round.(x; digits=3)) .=> names(df, Float64)) # round df to 3 digits

df = DataFrame(quantile(chain, q=[0.50, 0.025, 0.975]))
posterior_quantiles = transform(df, names(df, Float64) .=> (x -> round.(x; digits=3)) .=> names(df, Float64)) # round df to 3 digits

println(posterior_quantiles)

CSV.write("posterior_summary.csv", df)



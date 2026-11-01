### Define the module for MultispeciesOccModel

module MultispeciesOccModel


### Declare packages used throughout the module

using StatsBase, Distributions, TruncatedMVN, Random, Missings, LinearAlgebra, GLM, Optim, ForwardDiff

# using LogDensityProblems, ADTypes
# using LogDensityProblemsAD, AdvancedHMC     # needed by Hamiltonian Monte Carlo 
# import AdvancedHMC: phasepoint, transition  # make these private functions accessible

using TimerOutputs


### Declare functions or lists with public access

export simulateData, runMCMC, dataList

export t_out, reset_timer!



### Define lists


## ... define lists to hold simulated data values and parameter values

struct paramList
    
    alpha::Array{Float64}
    beta::Array{Float64}
    R::Array{Float64}
    
    Z::Array{Int}
        
end


struct dataList

    Y::Array{Union{Int, Missing}}
        
    X::Array{Float64}
    V::Array{Union{Float64, Missing}}
            
end


struct simulated_dataAndparams

    param::paramList
    data::dataList

end




## ... define list to hold Markov chain

struct markovChainList

    ## ... values of parameters
    
    vecAlpha::Array{Float64}
    vecBeta::Array{Float64}
    vecR::Array{Float64}

    vecZ::Array{Int}

    ## ... names of parameters
    
    vecAlpha_name::Vector{String}
    vecBeta_name::Vector{String}
    vecR_name::Vector{String}
    
    vecZ_name::Vector{String}
    
end




### Include source code that contains definitions of functions

include("MultispeciesOccModel_SimulateData.jl")
include("MultispeciesOccModel_MCMC.jl")

include("MultispeciesOccModel_FullConditionals.jl")  # this file also includes definitions of lists not defined or used elsewhere
# include("MH_Sampler.jl")   # Metropolis Hastings sampler
# include("HMC_Sampler.jl")  # Hamiltonian Monte Carlo sampler

end


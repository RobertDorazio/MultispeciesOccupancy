
###  CONDITIONAL POSTERIOR DENSITY FUNCTIONS (UNNORMALIZED)


###  Log of unnormalized pdf for conditional posterior of Alpha_j (= jth column of Alpha)

struct argList_Alpha
    
    yvec::Vector{Int}
    vmat::Matrix{Float64}

end



function logDensity_Alpha(alpha, arg)

    ## ... assign values of function's arguments
    yvec = arg.yvec
    vmat = arg.vmat

    ## ... compute prior density

    logPrior = sum( logpdf.(Normal(), alpha) )
    
    
    ## ... compute logLikelihood

    pvec = cdf.(Normal(), vmat * alpha)

    logLikelihood = sum( logpdf.(Bernoulli.(pvec), yvec) )
    
    
    return logLikelihood + logPrior
    
end



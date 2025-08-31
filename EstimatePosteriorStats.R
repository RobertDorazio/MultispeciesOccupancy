library(mcmcse)
library(coda)


### Define function for estimating posterior stats and Monte Carlo errors
###                 and for estimating effective size of Markov chains used to estimate posterior means

EstimatePosteriorStats = function(mc, burnin=0, thin=0, prob.quantiles=c(.50, .025, .975)) {
    
    mc.names = colnames(mc)
    
    ## ... if needed, remove burnin from Markov chain
    if (burnin>0) {
        out = -(1:burnin)
        mc = matrix(mc[out, ], ncol=ncol(mc))
    }
    
    ## ... if needed, thin Markov chain
    if (thin>0) {
        keep = seq(from=1, to=nrow(mc), by=thin)
        mc = matrix(mc[keep, ], ncol=ncol(mc))
    }

    
    ## ...estimate posterior means and quantiles

    prob.names = paste(as.character(100*prob.quantiles), '%', sep='')
    post.names = mc.names
    post.stats = matrix(nrow=ncol(mc), ncol=1+length(prob.quantiles))
    dimnames(post.stats) = list(post.names, c('Mean', prob.names))
    post.stats.MCSE = post.stats

    for (j in 1:ncol(mc)) {
        xvec = as.vector(mc[, j])
        mc.estimate = mcse(x=xvec, method='obm')
        post.stats[j, 1] = mc.estimate$est
        post.stats.MCSE[j, 1] = mc.estimate$se
        for (k in 1:length(prob.quantiles)) {
            mc.estimate = mcse.q(x=xvec, q=prob.quantiles[k], method='obm')
            post.stats[j, 1+k] = mc.estimate$est
            post.stats.MCSE[j, 1+k] = mc.estimate$se
        }
    }


    ## ...estimate effective size of Markov chain for each parameter

    mc.obj = as.mcmc(mc)
    effectiveSize = effectiveSize(mc.obj)
    names(effectiveSize) = mc.names

    nominalSize = rep(nrow(mc), ncol(mc))
    names(nominalSize) = mc.names

    post.stats.size = cbind(effectiveSize, nominalSize)
    colnames(post.stats.size) = c('Effective', 'Nominal')

    

    list(estimate=post.stats, MCerror=post.stats.MCSE, size=post.stats.size)
}

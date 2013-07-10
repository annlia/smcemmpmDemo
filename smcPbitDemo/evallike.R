##---------------------------------------------------------------------------------------------##
## Part of the demo code for Sequential Monte Carlo EM for multivariate probit models (SMCEMMPM)
## Giusi Moffa and Jack Kuipers
## University of Regensburg
#
## Last modified: July 4, 2013
#
## Disclaimer: The code in this archive is not guaranteed to be optimised or free of bugs.
##        Please report any issues to the authors (giusi.moffa@ur.de, jack.kuipers@ur.de).
##---------------------------------------------------------------------------------------------##

### Define a function to evaluate the likelihood of the estimates estbe and sigem
evallike <- function(obs, X, p, estbe, sigem, smcM = 1000){
  require(mvtnorm)
	### evaluate number of observations
	nobs <- length(obs[1,])
	### define constraints
	best <- matrix(0, nrow = p, ncol = nobs)
	for(i in 1:nobs) best[,i] <- X[[i]]%*%estbe
	regProb <- rep(NA, nobs)
	for(l in 1:nobs){
		estobs <- obseval(obs[,l], best[,l], sigem, smcM)
		regProb[l] <- estobs$prob
	}
	sum(log(regProb)) ### this is the log-likelihood value (which is returned)
}

likeParallel <- function(n,...) evallike(...)
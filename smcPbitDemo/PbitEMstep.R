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

### One standard (resampling from scratch at every iteration) step of the SMCEM procedure for a Probit Model 
### with observations obs of dimension p and covariates X (ncov for each observation)
### smcM is the number of particles used for the SMC sampler
### Require "obseval.R"
#######################################################################################

PbitEMstep <- function(obs, X, p, ncov, estbe = rep(0.5,ncov+1), sigem = diag(p), smcM = 1000){
	### evaluate number of observations
	nobs <- length(obs[1,])
	### define constraints
	best <- matrix(0, nrow = p, ncol = nobs)
	for(i in 1:nobs) best[,i] <- X[[i]]%*%estbe
# initialisations for estimation step
	estsq <- vector("list",nobs)
	regProb <- rep(NA, nobs)
	sumSquare <- vector("list",nobs)
	regMean <- matrix(NA, nrow=p, ncol=nobs)
	sqAvg <- diag(0,p)
# estimation step	
	for(l in 1:nobs){
		estobs <- obseval(obs[,l], best[,l], sigem, smcM)
		tmpSum <- t(t(estobs$samp)*estobs$W)
		totSum <- apply(tmpSum, 1, sum)
		regProb[l] <- estobs$prob
		estsq[[l]] <- (tmpSum)%*%t(estobs$samp)
		tmpMix <- totSum%*%t(best[,l])
		sumSquare[[l]] <- estsq[[l]]-tmpMix-t(tmpMix)+best[,l]%*%t(best[,l])
		regMean[,l] <- -totSum + best[,l]
		sqAvg <- sqAvg+estsq[[l]]
	}
	sqAvg <- sqAvg/nobs
	loglike <- sum(log(regProb)) ## log-likelihood (before parameter update)
	Qlike <- -nobs*1*(log(abs(det(sigem)))+sum(diag(solve(sigem)%*%sqAvg)))/2
	
	Mres <- Mstep(sigem, estbe, X, regMean, sumSquare, p, Qlike)
	sigem <- Mres$newSig
	estbe <- Mres$newBe

	return(list(loglike=loglike, estbe=estbe, sigem=sigem))
}

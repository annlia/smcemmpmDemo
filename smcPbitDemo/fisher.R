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

#######################################
### Function to	evaluate fisher     ###
### information for Probit models,  ###
###  according to Louis Formula     ###
#######################################

fisher <- function(be, sig, X, V, weight){
invsig <- solve(sig)
beN <- length(be)
begrad <- rep(0, beN)
p <- length(sig[,1])
sigN <- p*(p-1)/2 ### off diagonal elements

### build matrices for differentials of off diagonal elements
### are indices consistent with vector of covariances?
emptym <- matrix(0, nrow = p, ncol = p)
F <- list()
sp <- sigN
for(j in 1:(p-1)){
	sj <- sum(0:(p-j))
	dpj <- sp-sj
	for(i in (j+1):p){
		F <- c(F, list(emptym))
		m <- i-j+dpj
		F[[m]][i,j] <- 1
		F[[m]][j,i] <- 1 
	}
}

if(!constrained||useinvariance){ ### when no correlation form is imposed
  Fdiag <- vector("list",3) ### differentials for diagonal elements (escluding the first which is assumed fixed)
  for(j in 2:p) {
    Fdiag[[j-1]] <- emptym
    Fdiag[[j-1]][j,j] <- 1
  }
  F <- c(Fdiag,F)
  sigN <- sigN+p-1 ### modification for the case when only first element of covariance matrix is fixed 
  ### size of off-diagonal elements plus p-1 (one diagonal element is fixed)
}

### evaluate gradient and variance estimate (by MC) for a given obs, 
### for which the covariates X are passed to the function
gradsig <- rep(0, sigN) ### variable for estimated gradient
dimH <- beN+sigN
estvarquad <- matrix(0, nrow=dimH, ncol=dimH) ### variable for estimated variance
quadval <- matrix(0, nrow=dimH, ncol=dimH)
M <- length(weight)

amgisFamgis <- vector("list",sigN)
traceamgisF <- rep(0,sigN)	
for(k in 1:sigN){
	tempmatrix<-invsig%*%F[[k]]
	traceamgisF[k]<-sum(diag(tempmatrix))
	amgisFamgis[[k]]<-tempmatrix%*%invsig
}
	
	for(i in 1:M){
		tmpGrBe <- -(t(X)%*%invsig%*%V[,i])
		begrad <- begrad + weight[i]*tmpGrBe ### gradient of \beta for sample i
		tmpGrSig <- rep(0, sigN)
		for(j in 1:sigN) tmpGrSig[j] <- -traceamgisF[j] + t(V[,i])%*%amgisFamgis[[j]]%*%V[,i]
		gradsig <- gradsig + 0.5*weight[i]*tmpGrSig ### gradient of \sigma for sample i
		tmpGrad <- c(tmpGrBe, 0.5*tmpGrSig)
		quadval <- quadval + weight[i]*tmpGrad%*%t(tmpGrad)
	}
	gradvec <- c(begrad,gradsig) ### estimate of complete (for \beta and \sigma) gradient vector expected value
	
	weightedsampl<-t(t(V)*weight)
	mom1weighted <- apply(weightedsampl, 1, sum)
	mom2weighted <- weightedsampl%*%t(V)
	
grSq <- gradvec%*%t(gradvec)
estvarquad <- quadval - grSq ### complete variance estimate

### evaluate estimate for the Hessian for a given obs 
### for which the covariates X are passed to the function	

	beHess <- -t(X)%*%invsig%*%X
	besigHess <- matrix(0, nrow=beN, ncol=sigN)
	sigHess <- matrix(0, nrow=sigN, ncol=sigN)
	
	for(k in 1:sigN){
		besigHess[,k] <- t(X)%*%amgisFamgis[[k]]%*%mom1weighted
		for(j in 1:k){
			ph <- amgisFamgis[[k]]%*%F[[j]]
			sigHess[k,j] <- 0.5*sum(diag(ph)) - sum(diag(ph%*%invsig%*%mom2weighted))
		}
	}
	sigHess <- sigHess+t(sigHess) - diag(diag(sigHess))
	
Hess <- matrix(0, ncol = dimH, nrow = dimH)
Hess[1:beN,1:beN] <- beHess
Hess[(beN+1):dimH,(beN+1):dimH] <- sigHess
Hess[1:beN,(beN+1):dimH] <- besigHess
Hess[(beN+1):dimH,1:beN] <- t(besigHess)

		
fish <- list()
fish$info <- -Hess-estvarquad
fish$Hess <- Hess
fish$quad <- estvarquad
fish$quval <- quadval
fish$gradsq <- grSq
return(fish)
}

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

smcEM <- function(obs, bound, prev.bound, s, prev.s, prevSamp, sf=1){
	prob <- prevSamp$prob
	bigW <- prevSamp$W
	x <- prevSamp$samp
	N <- length(bigW)
	scal.fac <- sf
	
####################################################
### It might be useful not to accept the whole move 
### to avoid degeneration in the smcEM sample
### Iterate until real target is reached...	
#sigem <- .1*sqAvg+.9*sigem
#estbe <- .1*beAvg+.9*estbe
############################	

	
	A <- 1-2*obs
	th <- 0.8*N
  mydim <- length(obs) # problem dimension
  chs <- t(chol(s))

	# Constraints to just do rejection sampling instead 
	mov.dist <- A*(pmin(qnorm(1-(.25)^(1/mydim)),A*bound))
  bound.bool <- mov.dist == bound
	not.reached <- which(!bound.bool)

	if(all(bound.bool)){ ## the region is obvious

		xtemp <- matrix(rnorm(mydim*N), nrow = mydim, ncol = N)
		xtemp <- chs%*%xtemp ## draw initial sample
		accept <- which(apply(A*xtemp>A*mov.dist, 2, all))
		xy <- xtemp[,accept]

		mmm<-1
		while(ncol(xy)<N){
			xtemp <- matrix(rnorm(mydim*N), nrow = mydim, ncol = N)
			xtemp <- chs%*%xtemp
			accept <- which(apply(A*xtemp>A*mov.dist, 2, all))
			xy <- cbind(xy,xtemp[,accept])
			mmm<-mmm+1
		}	

    ### Estimate of the probability of the region
		Pest <- ncol(xy)/(mmm*N)
		pos.left <- seq(1:N)
		x<-xy[,pos.left]
		bigW <- rep(1/N, length = N)
		scal.fac<-1
	}

else{
  ### rescale to avoid the situation of practically having to move to a bigger region
	scale.bound <- diag(bound/prev.bound)
	if(any(!diag(scale.bound)>0)) stop("Re-scaling is failing, need resampling from scratch") 
  ### !!! modify so that if rescaling is not possible sampling is performed from scratch again
  ### only for this iteration and this particular observation
	prev.s <- scale.bound%*%prev.s%*%scale.bound
	for(j in 1:N) x[,j] <- scale.bound%*%x[,j]

	pos.left <- which(as.logical(bigW)) ### samples with weight different than zero
#	ESS <- 1/(sum(bigW^2)) # Effective sample size
	Pest <- ((det(prev.s)/det(s))^(.5))*prob
	s.inv <- solve(s)
	prev.sinv <- solve(prev.s)

	quadfvec<-rep(0,N)
	for(j in 1:N){
		quadfvec[j]<-x[,j] %*% s.inv %*% x[,j]
	}

	smallw <- rep(0, length = N) # reset weights
	for(j in pos.left){
		if (all(A * x[,j] > A * bound)){
			numw <- quadfvec[j]
			denw <- x[,j] %*% prev.sinv %*% x[,j]
			smallw[j] <- exp(.5*(denw-numw))
		}
	}
	Zratio <- bigW %*% smallw
	if(!Zratio) stop("Sample degenerated")
	Pest <- Pest*Zratio
	bigW <- (bigW * smallw)/Zratio # Weight normalization
	pos.left <- which(as.logical(bigW))	
	ESS <- 1/(sum(bigW^2)) # Effective sample size

	if(ESS < th){ # Condition for resampling
		resampled<-resample.efficient(bigW)
		x<-x[,resampled]
		quadfvec<-quadfvec[resampled]
	#	smallw <- rep(1, N)
		bigW <- rep(1/N, N)
		pos.left <- seq(1:N)
	}

	for(h in 1:mh){
		mhacc <- 0 ### initialise mean acceptance probability for current generation of particles
		altprop<-chs%*%matrix(rnorm(mydim*N),nrow=mydim,byrow=T)+x
		for(j in pos.left){
			quadf <- quadfvec[j]
			if( all(A * altprop[,j] > A * bound)){
				propqf <- altprop[,j] %*% s.inv %*% altprop[,j]
				alpha <- exp(.5*(quadf - propqf))
				if( runif(1) < min(alpha,1) ){
					 x[,j] <- altprop[,j]
					quadfvec[j]<-propqf
				}
			mhacc <- mhacc + bigW[j]*min(alpha,1)
			}
		}
	}

	log.fac <- log(scal.fac) + sf.ga*(mhacc - alpha.wanted) ### scaling factor adapting step
	scal.fac <- exp(log.fac)
}
	
	res <- list()
	res$prob <- Pest
	res$W <- bigW
	res$samp <- x
	res$sf <- scal.fac
	return(res)
}

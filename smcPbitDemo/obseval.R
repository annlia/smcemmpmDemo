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

## OBJECT: funtion to sample from multivariate normal
## with covariance matrix s, constrained to the region
## defined by obs (vector \gamma of p component) 
## and bound (vector b of p component), 
## Z_i > b_i if \gamma_i = 0 and Z_i < b_i if \gamma_i = 1

obseval <- function(obs, bound, s, smcM=1000){

	countit <- 0

	mydim <- length(obs) # problem dimension
	N <- smcM # Sample size
	s.inv <- solve(s) # Inverse of covariance matrix
	chs <- t(chol(s)) # Cholesky decomposition

	# Define domain parameters, matrix of contraints
	A <- 1-2*obs ## values to define the direction of the box

	bound.bool <- rep(FALSE, mydim) # vector that keeps track of boundaries reached
	not.reached <- seq(1:mydim) # indices of boundaries not reached

	# Initialize parameters for the moving domain
	# Starting constraints from boundaries taking care of high dimensional issues 
	# 0.25 below would leave a quarter of particles if covariance matrix were identity
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
	} 

	else {	### initializations ###
		smallw <- rep(1, length = N) # Initialise weights
		pos.left <- seq(1:N) # Initialise vector of indices of weight left after SMC
		bigW <- rep(1/N, length = N) # Normalized weights
		ESS <- 1/(sum(bigW^2)) # Effective sample size, keep values along directions
		
		### Define working variables (steps) defining moving targets ###
		step <- rep(0.2, mydim)
		min.step <- step/10
		nu.step <- 1 ## 0.6
		min.nustep <- nu.step/10
		stopNu <- 83

		th <- 0.8*N ### resampling and adaptive threshold
		th.ga <- 2/N ### 0 No adaptation, otherwise something like: 2/N
		
		scal.fac <- 1 ### MH kernel scaling factor
		log.fac <- 0 ### work on log to ensure scaling factor stays positive
		alpha.wanted <- 0.6 ### target acceptance probability for control of adaptive MH
		sf.ga <- 7 ### Stepsize to adapt the scaling factor
		mh <- 4 ### number of mcmc steps within smc
		################################################################################
		
		### Define initial distribution, student-t ###
		nu <- 3 	# student t degrees of freedom
		### note, nu=1 coincides with Cauchy distribution, no variance
		Pest <- ((nu/2)^(mydim/2)*gamma(nu/2))/(gamma((nu+mydim)/2))
		### initialise probability region estimate, taking into account ###
		### that the first sample comes from a normalized distribution ###
		##################################################################
        sdf<-1

		xtemp <- t( rmvt(N, sigma = sdf*s, df = nu)) ## draw initial sample
		accept <- which(apply(A*xtemp>A*mov.dist, 2, all))
		xy <- xtemp[,accept]

		mmm<-1
		while(ncol(xy)<N){
			xtemp <- t( rmvt(N, sigma = sdf*s, df = nu)) ## draw initial sample
			accept <- which(apply(A*xtemp>A*mov.dist, 2, all))
			xy <- cbind(xy,xtemp[,accept])
			mmm<-mmm+1
		}	

### Estimate of the probability of the starting region
		Pstart <- ncol(xy)/(mmm*N)
		Pest <- Pstart*Pest

		x<-xy[,pos.left]

		quadfvec<-rep(0,N)

		for(j in pos.left){
			quadfvec[j]<-x[,j] %*% s.inv %*% x[,j]
		}

        
#print('Moving')
		
		## run the SMC update for the t distribution and moving the support
		ni.inv <- 1/nu
		esp <- (nu + mydim)/2
		repeat{
			for(i in not.reached){ ### Move along one constraint at a time
				mov.dist[i] <- A[i]*min( A[i]*mov.dist[i] + step[i], A[i]*bound[i]) # Move box
				smallw <- rep(0, N) ## reset weights
				for(j in pos.left){
					if (all(A * x[,j] > A * mov.dist)) smallw[j] <- 1
				}
				Zratio <- bigW %*% smallw
				if(!Zratio) { 
					print(paste("Error, step Degenerated zero", countit))
					break }			
				Pest <-Pest*Zratio
				bigW <- (bigW * smallw)/Zratio # Weight normalization
				pos.left <- which(as.logical(bigW))
				ESS <- 1/(sum(bigW^2))
				step[i] <- max(step[i] + th.ga*(ESS-th), min.step[i])
				if(ESS < th){ # Condition for resampling
					#x <- resamplegiusi(mydim, N, x, bigW, pos.left)
					resampled<-resample.efficient(bigW)
					x<-x[,resampled]
					quadfvec<-quadfvec[resampled]
					bigW <- rep(1/N, N)
					pos.left <- seq(1:N)
				}
				mhacc <- 0 ### initialise mean acceptance probability for current generation of particles

				altprop<-chs%*%matrix(rnorm(mydim*N),nrow=mydim,byrow=T)+x

				for(j in pos.left){ ### MH step
					if( all(A * altprop[,j] > A * mov.dist)){
						propqf <- altprop[,j] %*% s.inv %*% altprop[,j]
						alpha <- ((1 + ni.inv*(quadfvec[j]))/(1 + ni.inv*propqf))^(esp)
						if( runif(1) < min(alpha,1) ){
							x[,j] <- altprop[,j]
							quadfvec[j]<-propqf
						}
						mhacc <- mhacc + bigW[j]*min(alpha,1) 
						### evaluate mean acceptance probability for current generation of particles 
					}
				}
				log.fac <- log.fac + sf.ga*(mhacc - alpha.wanted) ### scaling factor adapting step
				scal.fac <- exp(log.fac)
			countit <- countit + 1
			}
			if(all(bound.bool)) break
			bound.bool <- mov.dist == bound
			not.reached <- which(!bound.bool)
		}

#print('Resampling')

		resampled<-resample.efficient(bigW)
		x<-x[,resampled]
		quadfvec<-quadfvec[resampled]

		smallw <- rep(1, N)
		bigW <- rep(1/N, N)
		pos.left <- seq(1:N)

#print('Changing nu')

		# Run SMC to move from student to normal
		# note, pos.left in this case should always be 1:N
		while(nu<stopNu){
			smallw <- rep(0,N)
			prev.free <- nu
			nu <- nu+nu.step # Increase degree of freedom
			for(j in pos.left){
				ni.inv <- 1/prev.free
				esp <- (prev.free + mydim)/2
				numw <- (1 + ni.inv*quadfvec[j])^(esp)
				ni.inv <- 1/nu
				esp <- (nu + mydim)/2
				denw <- (1 + ni.inv*quadfvec[j])^(esp)
				smallw[j] <- numw/denw
			}
			Zratio <- bigW %*% smallw
			if(!Zratio) { print(paste("Error, nu Degenerated zero", countit))
				break }
			Pest <- Pest*Zratio
			bigW <- (bigW * smallw)/Zratio # Weight normalization
			pos.left <- which(as.logical(bigW))		
			ESS <- 1/(sum(bigW^2)) # Effective sample size
			nu.step <- max(nu.step + th.ga*(ESS-th), min.nustep)

			if(ESS < th){ # Condition for resampling
				resampled<-resample.efficient(bigW)
				x<-x[,resampled]
				quadfvec<-quadfvec[resampled]
				smallw <- rep(1, N)
				bigW <- rep(1/N, N)
				pos.left <- seq(1:N)
			}
			mhacc <- 0 ### initialise mean acceptance probability for current generation of particles

			altprop<-chs%*%matrix(rnorm(mydim*N),nrow=mydim,byrow=T)+x

			for(j in pos.left){
				if( all(A * altprop[,j] > A * mov.dist)){
					propqf <- (altprop[,j] %*% s.inv %*% altprop[,j])
					alpha <- ((1 + ni.inv*quadfvec[j])/(1 + ni.inv*propqf))^(esp)
					if( runif(1) < min(alpha,1) ){
						x[,j] <- altprop[,j]
						quadfvec[j]<-propqf
					}
					mhacc <- mhacc + bigW[j]*min(alpha,1)
				}
			}
			log.fac <- log.fac + sf.ga*(mhacc - alpha.wanted) ### scaling factor adapting step
			scal.fac <- exp(log.fac)
			countit <- countit + 1
		}

		ni.inv <- 1/nu
		esp <- (nu + mydim)/2
		smallw <- rep(0,N)
		for(j in pos.left){
			smallw[j] <- (exp(-.5*quadfvec[j]))*(1 + ni.inv*quadfvec[j])^(esp)
		}
		Zratio <- bigW %*% smallw
		if(!Zratio) { print(paste("Error, last Degenerated zero", countit))
				break }
		Pest <-Pest*Zratio
		bigW <- (bigW * smallw)/Zratio
		altprop<-chs%*%matrix(rnorm(mydim*N),nrow=mydim,byrow=T)+x
		for(j in pos.left){
			if( all(A * altprop[,j] > A * mov.dist)){
				propqf <- (altprop[,j] %*% s.inv %*% altprop[,j])
				alpha <- exp(.5*(quadfvec[j]-propqf))
				if( runif(1) < min(alpha,1) ){
					x[,j] <- altprop[,j] 
					quadfvec[j]<-propqf
				}
			}
		}
		countit <- countit + 1

	}


	res <- list()
	res$prob <- Pest
	res$W <- bigW
	res$samp <- x
	res$count <- countit
	res$sf <- 1
	return(res)
}

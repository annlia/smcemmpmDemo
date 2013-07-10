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

### pairwise estimate of the latent variable correlation matrix
### the matrix is evetually bended to make sure it is positive defined
### the result will be used to initialise the covariance matrix of an EM algorithm

### the function takes as inputs the mean vector and the correlation matrix
### of the binary observations

siginit <- function(avg.g, cor.g){
	nc <- length(avg.g)
	pairest <- diag(nc)
	tol <- 1e-9

	# Evaluate normal distribution quantiles
	qz <- qnorm(avg.g)
	#
	for(i in 2:nc){
		for(j in 1:(i-1)){
			# Check that necessary constraints on the correlations ar satisfied
			# and if not apply threshold
			m <- max( -((avg.g[i]*avg.g[j])/((1-avg.g[i])*(1-avg.g[j])))^(1/2), -( ((1-avg.g[i])*(1-avg.g[j]))/(avg.g[i]*avg.g[j]) )^(1/2) )
			M <- min( ((avg.g[i]*avg.g[j])/((1-avg.g[i])*(1-avg.g[j])))^(1/2), ( ((1-avg.g[i])*(1-avg.g[j]) )/(avg.g[i]*avg.g[j]))^(1/2) )
			del <- (M-m)/10^9
			if(cor.g[i,j]<m){
			cor.g[i,j] <- m+del
			cor.g[j,i] <- cor.g[i,j]
			} else {
				if(cor.g[i,j]>M){
				cor.g[i,j] <- M-del
				cor.g[j,i] <- cor.g[i,j]
				}
			}
			#
			# when we are sure the correlations satisfy the necessary constraints, 
			# we can solve the equation
			# \Phi(qz[i], qz[j], rho[i,j] = corr[i,j]*(p[i]*p[j]*(1-p[i])*(1-p[j]))^(1/2) + p[i]*p[j]
			# Bisection, define the interval [a,b] as [-1,1]
			d <- cor.g[i,j]*sqrt( avg.g[i]*avg.g[j]*(1-avg.g[i])*(1-avg.g[j]) ) + avg.g[i]*avg.g[j]
			a <- -(1-10^(-7))			
			b <- (1-10^(-7))
			count <- 1
			sigma <- diag(2)
			repeat{
				sigma[upper.tri(sigma)] <- a
				sigma[lower.tri(sigma)] <- a
				l1 <- pmvnorm(lower = -Inf, upper= c(qz[i],qz[j]), mean <- c(0,0), sigma = sigma)
				sigma[upper.tri(sigma)] <- b
				sigma[lower.tri(sigma)] <- b
				l2 <- pmvnorm(lower = -Inf, upper= c(qz[i],qz[j]), mean <- c(0,0), sigma = sigma)
				sigma[upper.tri(sigma)] <- (a+b)/2
				sigma[lower.tri(sigma)] <- (a+b)/2
				l3 <- pmvnorm(lower = -Inf, upper= c(qz[i],qz[j]), mean <-  c(0,0), sigma = sigma)
				#
				if ((l1[1]-d)*(l3[1]-d) < 0){
					v1 <- l1[1]
					v2 <- l3[1]
					b <- (a+b)/2
					} else {
						v1 <- l3[1]
						v2 <- l2[1]
						a <- (a+b)/2
						}
				count <- count+1 # we define a counter to avoid an infinite loop
				if( abs(v1-v2) < 10^(-12) || (count > 10^6)) break
				}
			rho <- (a+b)/2
			pairest[i,j] <- rho
			pairest[j,i] <- rho
			}
		}
	#
	# matrix bending
	peig <- eigen(pairest, symmetric = TRUE)
	if(!all(peig$values > 0)){
		cat("bending active \n")
		rtol <- tol * peig$values[1]
		if(min(peig$values) < rtol) {
			vals <- peig$values
			vals[vals < rtol] <- rtol
			srev <- peig$vectors %*% (vals * t(peig$vectors))
			dimnames(srev) <- dimnames(pairest)
			pairest <- srev
			}
		}
return(pairest)
}
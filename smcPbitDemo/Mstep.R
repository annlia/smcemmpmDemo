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

Mstep <- function(sqAvg, beAvg, X, regMean, sumSquare, p, Qlikeold){
	convergencelimit <- 10^(-12)
	looplimit <- 30
	Qlike <- Qlikeold+10*convergencelimit
	countit <- 0

	dm1<-1 ### inverse of overall scale factor 
	if(constrained||useinvariance){
		Omega<-diag(1,p)
		### better starting point is previous value
		Omega<-sqAvg
		if(!constrained){	### better starting point is previous value for diagonal scaling
			dm1<-1/sqrt(sqAvg[1,1])
		}
	}	

	while(((Qlike-Qlikeold)>convergencelimit)&&(countit<looplimit)){
		Qlikeold <- Qlike
		countit <- countit+1
		### Update \beta ###
		tmpInv <- solve(sqAvg)
		tmpqf <- 0
		tmpXmuQf <- 0
		for(l in 1:nobs){
			tempy <- t(X[[l]])%*%tmpInv
			tmpqf <- tmpqf + tempy%*%X[[l]]
			tmpXmuQf <- tmpXmuQf + tempy%*%regMean[,l]
		}
		beAvg <- solve(tmpqf)%*%tmpXmuQf/dm1
		### Update \Sigma
		sqAvg <- diag(0,p)
		s1 <- diag(0,p)
		s2 <- diag(0,p)
		s3 <- diag(0,p)
		for(l in 1:nobs){
			s1 <- s1+sumSquare[[l]]
			tempy <- X[[l]]%*%beAvg
			s2 <- s2-tempy%*%t(regMean[,l])
			s3 <- s3+tempy%*%t(tempy)
		}
		s1<-s1/nobs
		s2<-s2/nobs
		s3<-s3/nobs
		sqAvg <- s1+dm1*s2+dm1*t(s2)+dm1*dm1*s3
		if(!constrained&&!useinvariance){
			Qlike <- -nobs*1*(log(abs(det(sqAvg)))+p)/2
		}
		if(constrained||useinvariance){

			countitdm1 <- 0
			testdm1 <- 10*convergencelimit

			while((testdm1>convergencelimit)&&(countitdm1<looplimit))
			{
				countitdm1 <- countitdm1 + 1
				dm1tmp <- dm1

				countitconstrained <- 0
				testz<-10*convergencelimit

				while((testz>convergencelimit)&&(countitconstrained<looplimit)){
					countitconstrained<-countitconstrained+1
					Omegatemp<-Omega
					if(!useinvariance){
						temp <- solve(Omegatemp*Omegatemp)
						temp2 <- diag(diag(1,p)-sqAvg)
						temp3 <- temp%*%temp2
						A<-diag(c(temp3))
					}
					if(useinvariance){
						temp <- 1/(Omegatemp[1,1]^2)
 						temp2 <- 1-sqAvg[1,1]
						temp3 <- temp*temp2
						A<-diag(0,p)
						A[1,1]<-temp3
					}
					Omega<-sqAvg+Omegatemp%*%A%*%Omegatemp
					testy<-Omega-Omegatemp
					testz<-sqrt(sum(testy*testy))/p
				} 

				if(testz>convergencelimit){
					print("No convergence!") # Project instead
					if(!useinvariance){
						sqdiag <- sqrt(diag(1/diag(sqAvg)))
					}
					if(useinvariance){
						sqdiag <- diag(1/sqrt(sqAvg[1,1]))
					}
					Omega <- sqdiag%*%sqAvg%*%sqdiag
				}

				Qlike <- -nobs*1*(log(abs(det(Omega)))+sum(diag(solve(Omega)%*%sqAvg)))/2

				if(!constrained){
					tempy <- solve(Omega)
					dm1 <- -(A[1,1]-p+sum(diag(tempy%*%s1)))/(sum(diag(tempy%*%s2)))
					sqAvg <- s1+dm1*s2+dm1*t(s2)+dm1*dm1*s3
				}

				testdm1 <- abs(dm1-dm1tmp)

			}
			sqAvg <- Omega/(dm1^2)
		}
	}
	if(doProjection){
	  sqdiag <- 1/sqrt(sqAvg[1,1])
	  sqAvg <- sqdiag*sqdiag*sqAvg ### project covariance so first element is 1
	  beAvg <- sqdiag*beAvg ### rescale beta too
	}

	if(!doProjection&&!constrained){
	  sqdiag <- 1/sqrt(sqAvg[1,1])
	}
	
	res <- list()
	res$newSig <- sqAvg
	res$newBe <- beAvg
	res
}

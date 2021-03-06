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

### From withing the working folder
rm(list=ls(all.names=TRUE))
library("mvtnorm")
pkgfolder <- "../smcPbitDemo/"
datafile <- "sixcities.RData" ### file with data to model
tosave <- TRUE
savefile <- "smcSixCdemo11062013.RData" ### file to save R workspace after processing
### From within the working directory
source(paste0(pkgfolder, "resample.R"))
source(paste0(pkgfolder, "siginit.R"))
source(paste0(pkgfolder, "obseval.R"))
source(paste0(pkgfolder, "PbitEMstart.R"))
source(paste0(pkgfolder, "PbitEMloop.R"))
source(paste0(pkgfolder, "Mstep.R"))
source(paste0(pkgfolder, "smcEM.R"))
source(paste0(pkgfolder, "evallike.R"))
source(paste0(pkgfolder, "fisher.R"))
load(datafile) ### Read data

doProjection <- FALSE # projection to sigma_11=1 form is performed after M step - projection is not needed if constrained is TRUE
constrained <- TRUE # Maximisation is performed constrained - either in correlation form or with sigma_11=1 depending on useinvariance
useinvariance <- FALSE # Maximisation uses invariant Q tilde function instead of standard Q
fixM <- FALSE # keep the number of particles constant or not
nr <- 40
M <- 1000
refineM <- 4*M
startM <- 100
nc <- length(citygam[,1])
nobs <- length(citygam[1,])
nX <- length(cityX[[1]][1,])-1
avg.gam <- apply(citygam, 1, mean)  ### mean vector
cor.gam <- cor(t(citygam))
set.seed(101)

pwsig <- siginit(avg.gam, cor.gam) ### covariance matrix pairwise estimate
updsig <- pwsig
obs <- c(t(citygam))
covmat <- matrix(0, nrow = nc*nobs, ncol=nX)
covmataux <- matrix(0, nrow= nobs, ncol = nX)
for(i in 1:nc){
  for(j in 1:nobs) covmataux[j,] <- cityX[[j]][i,2:(nX+1)]
  covmat[((i-1)*nobs+1):(i*nobs),] <- covmataux
}
indprobit <- glm(obs ~ covmat, family=binomial(link="probit"), na.action=na.pass)
indbe <- coef(indprobit) ### coefficients from fitting independent probit models
updbe <- indbe

loglike <- rep(NA, nr)
regcoeff <- list()
regcoeff[[1]] <- updbe
covest <- list()
covest[[1]] <- updsig

res <- PbitEMstart(citygam, cityX, nc, nX, updbe, updsig, M)
loglike[1] <- res$loglike
### record sigma and beta at previous step for sequential sampling
prevBe <- updbe
prevSig <- updsig
### update sigma and beta with new estiamte
updbe <- res$estbe
updsig <- res$sigem
covest[[2]] <- updsig
regcoeff[[2]] <- updbe

### sequential monte carlo EM ###

### parameters for kernel adaptation, sf defined in smcEM
log.fac <- 0 ### factor for MH kernel scale adaptation
alpha.wanted <- 0.3 ### target acceptance probability for control of adaptive MH
sf.ga <- 7 ### stepsize to adapt the scaling factor
mh <- 4 ### number of mcmc steps within smc

for (k in 2:nr){
  res <- PbitEMloop(citygam, cityX, nc, nX, updbe, prevBe, updsig, prevSig, M)
  loglike[k] <- res$loglike
  ### record sigma and beta at previous step for sequential sampling
  prevBe <- updbe
  prevSig <- updsig
  updbe <- res$estbe ### update beta
  updsig <- res$sigem ### update sigma
  covest[[k+1]] <- updsig
  regcoeff[[k+1]] <- updbe
}
sixCres <- list(loglike=loglike, estbe=updbe, estsig=updsig, regcoeff=regcoeff, covest=covest)
if(tosave) save.image(file = savefile)

# ### evaluate loglikelihood
evallike(citygam, cityX, nc, updbe, updsig, refineM)
# 
### evaluate standard errors
best <- matrix(0, nrow = nc, ncol = nobs)
for(i in 1:nobs) best[,i] <- cityX[[i]]%*%updbe

if(!constrained||useinvariance) 
  dimInf <- length(updbe)+nc*(nc+1)/2-1 else dimInf <- length(updbe)+nc*(nc-1)/2
totinfo <- matrix(0, nrow=dimInf, ncol=dimInf)
totHess <- matrix(0, nrow=dimInf, ncol=dimInf)
totestquad <- matrix(0, nrow=dimInf, ncol=dimInf)
totGradSq <- matrix(0, nrow=dimInf, ncol=dimInf)
totquval <- matrix(0, nrow=dimInf, ncol=dimInf)

for(j in 1:nobs){
  samp <- obseval(citygam[,j], best[,j], updsig, refineM)
  obsInfo <- fisher(updbe, updsig, cityX[[j]], samp$samp, samp$W)
  totquval <- totquval+obsInfo$quval ### total second moment expectation
  totGradSq <- totGradSq + obsInfo$gradsq ### total expectation of gradient quadratic form
  totestquad <- totestquad + obsInfo$quad ### total variance estimate
  totHess <- totHess + obsInfo$Hess ### total hessian estimate
  totinfo <- totinfo + obsInfo$info ###
}

myse <- sqrt(diag(solve(totinfo)))

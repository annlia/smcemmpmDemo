simProbit <- function(bb, ss, nn){
  require(mvtnorm)
  ### inputs ###
  # bb: regression coefficients
  # ss: covariance matrix
  # nn: number of onservations
  
  ### outputs ###
  # X: simulated covariates 
  # L: latent variable means
  # Z: simulated latent variables
  # gam: simulated response variable   
    
  # observations are packed into a list, probably not the best,
  # but a number of connected functions should be updated to change
  # this representation
  pp <- dim(ss)[1] # dimension of response variable
  cc <- length(bb) # number of covariates (including intercept)
  simX <- vector("list", nn) ### build a list with nobs elements
  for(i in 1:nn) simX[[i]] <- matrix(c(rep(1,pp), 
                                         runif(pp*(cc-1), -.5, .5)), 
                                       nrow = pp, ncol = cc)
  ### latent variable mean
  ll <- sapply(simX, function(xx,rr) xx%*%rr, bb)
  zz <- t(rmvnorm(nn, sigma=ss)) + ll
  gam <- zz>0
  mode(gam) <- "integer"
  res <- list()
  res$L <- ll
  res$Z <- zz
  res$X <- simX
  res$gam <- gam
  res
}

### define model with response variable of size 8 (the number of covariates can be different)
library(mvtnorm)
set.seed(7)
pbitcoeff <- c(1, .3, -.3 , .2, -.2, .1, -.1)
pbitcov <- matrix(0,8,8)
pbitcov[upper.tri(pbitcov)] <- pbitcov[upper.tri(pbitcov)] <- c(rep(.1,8), rep(.2,8), rep(.3,4), rep(.4,4), rep(.5,2), rep(.6,2))
pbitcov <- pbitcov+t(pbitcov)
diag(pbitcov) <- c(1, 1.2, 1.2, 1.1, 1.1, .9, .9, .8)
# chol(pbitcov) # check positive-definitiveness
testPbit <- simProbit(pbitcoeff, pbitcov, 1000)
# apply(testPbit$gam, 1, mean)
citygam <- testPbit$gam
cityX <- testPbit$X
save(citygam, cityX, pbitcoeff, pbitcov, file="simPbit8Data1k.RData")
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

### Build the data set for the six cities example

### response variable: citygam
citygam <- matrix(NA, nrow=4, ncol=537)

## nosmoke <- c(237,10,15,4,16,2,7,3,24,3,3,2,6,2,5,11)
## smoke <- c(118,6,8,2,11,1,6,4,7,3,3,1,4,2,4,7)


### no maternal smoking

for(j in 1:237) citygam[,j] <- rep(0,4)
for(j in 238:247) citygam[,j] <- c(0,0,0,1)
for(j in 248:262) citygam[,j] <- c(0,0,1,0)
for(j in 263:266) citygam[,j] <- c(0,0,1,1)

for(j in 267:282) citygam[,j] <- c(0,1,0,0)
for(j in 283:284) citygam[,j] <- c(0,1,0,1)
for(j in 285:291) citygam[,j] <- c(0,1,1,0)
for(j in 292:294) citygam[,j] <- c(0,1,1,1)

for(j in 295:318) citygam[,j] <- c(1,0,0,0)
for(j in 319:321) citygam[,j] <- c(1,0,0,1)
for(j in 322:324) citygam[,j] <- c(1,0,1,0)
for(j in 325:326) citygam[,j] <- c(1,0,1,1)

for(j in 327:332) citygam[,j] <- c(1,1,0,0)
for(j in 333:334) citygam[,j] <- c(1,1,0,1)
for(j in 335:339) citygam[,j] <- c(1,1,1,0)
for(j in 340:350) citygam[,j] <- c(1,1,1,1)

### maternal smoking

for(j in 351:468) citygam[,j] <- rep(0,4)
for(j in 469:474) citygam[,j] <- c(0,0,0,1)
for(j in 475:482) citygam[,j] <- c(0,0,1,0)
for(j in 483:484) citygam[,j] <- c(0,0,1,1)

for(j in 485:495) citygam[,j] <- c(0,1,0,0)
for(j in 496:496) citygam[,j] <- c(0,1,0,1)
for(j in 497:502) citygam[,j] <- c(0,1,1,0)
for(j in 503:506) citygam[,j] <- c(0,1,1,1)

for(j in 507:513) citygam[,j] <- c(1,0,0,0)
for(j in 514:516) citygam[,j] <- c(1,0,0,1)
for(j in 517:519) citygam[,j] <- c(1,0,1,0)
for(j in 520:520) citygam[,j] <- c(1,0,1,1)

for(j in 521:524) citygam[,j] <- c(1,1,0,0)
for(j in 525:526) citygam[,j] <- c(1,1,0,1)
for(j in 527:530) citygam[,j] <- c(1,1,1,0)
for(j in 531:537) citygam[,j] <- c(1,1,1,1)

### there are only two kinds of covariates, corresponding to smoking and no smoking condition
citycov <- list()

ones <- rep(1,4)
zeros <- rep(0,4)
age <- c(-2,-1,0,1)

citycov$nosmoke <- matrix(c(ones, age, zeros, zeros), nrow=4)
citycov$smoke <- matrix(c(ones, age, ones, ones*age), nrow=4)

### list of covariates, for each observation
cityX <- list()
for(j in 1:350) cityX[[j]]<- citycov$nosmoke
for(j in 351:537) cityX[[j]] <- citycov$smoke

save.image("sixcities.RData")

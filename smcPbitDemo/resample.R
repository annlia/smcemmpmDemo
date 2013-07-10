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

resample.efficient <- function(weighty){

	num.samples <- length(weighty)
	n.exp.samps <- num.samples * weighty
        det.reps    <- floor(n.exp.samps)
	num.reps <- sum(det.reps)
        id <- integer(num.reps)
        current.index <- integer(1)
        for (i in 1:num.samples) 
        {
            if (det.reps[i] != 0) 
            {
                id[current.index + 1:det.reps[i]] <- i
                current.index <- current.index + det.reps[i]
            }
        }
        num.samples <- num.samples - num.reps
	weighty <- n.exp.samps-det.reps
	weighty <- weighty/sum(weighty)

if(num.samples>0){

        lbs <- seq(0, by=1/num.samples, length=num.samples)
        ubs <- lbs+1/num.samples
        uniforms <- runif(num.samples, lbs, ubs)
	#uniforms <- lbs+runif(num.samples)/num.samples

        ids       <- integer(num.samples)
        cusum     <- cumsum(weighty)

        current.index     <- 1
        for (i in 1:num.samples) 
        {
            found <- FALSE
            while (!found) 
            {
                if (uniforms[i] > cusum[current.index]) 
                {
                    current.index <- current.index + 1
                }
                else 
                {
                    found <- TRUE
                }
            }
            ids[i] <- current.index
        }

	id<-c(id,ids)
}

return(id)

}

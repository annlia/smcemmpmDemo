    ## Demo code for Sequential Monte Carlo EM for multivariate probit models (SMCEMMPM)
    #
    ## Giusi Moffa and Jack Kuipers
    ## University of Regensburg
    #
    ## Last modified: July 4, 2013
    #
    ## Disclaimer: The code in this archive is not guaranteed to be optimised or free of bugs.
    ##        Please report any issues to the authors (giusi.moffa@ur.de, jack.kuipers@ur.de).
    
****

#### Notes to the SMCEMMPM demo

The supplement `kntiSMCEMpbit.pdf` to the SMCEMMPM article contains detailed procedures, including the code, to fit a probit model to the two datasets considered in the paper.
    
* Multivariate Probit model fit of the six cities dataset, with covariance matrix in correlation form, for a more direct comparison to the existing literature.   
  The six cities analysis can be reproduced by running the file `sixCseqDemo.R` in the `smcPbitDemo` folder.
  
* Multivariate Probit model of a simulated dataset, with response variable of dimension 8.   
  The simulated data analysis can be reproduced by running the file `smcPbitDemo.R` in the `smcPbitDemo` folder.
    
### Files in the archive, `smcPbitDemo` folder
****
****

****
#### Demo data files:
****

__simPbit8Data1k.RData__ Contains simulated data from a multivariate probit model with response variable of dimension 8

__sixcities.RData__ Contains the six cities dataset

__smcemPbitDemo11062013.RData__ Contains image file with results from the analysis of the data in `simPbit8Data1k.RData` (as obtained by `smcemPbitDemo.R`)

__smcSixCdemo11062013.RData__ Contains image file with results from the analysis of the data in `sixcities.RData` (as obtained by `sixCseqDemo.R`)

 __samples__ Folder needed to save results from one iteration of the SMC EM to the next. An empty folder with this name needs to be present in the working folder before starting the analysis
****

****
#### Demo code files:
****

__buildsixcities.R__ Formats the six cities data as needed for the analysis

__evallike.R__ Evaluates the log-likelihood of given regression coefficients and covariance matrix 

__fisher.R__ Evaluates fisher information of multivariate probit models according to Louis Formula

__Mstep.R__ Performs the M step of the SMCEM procedure for multivariate probit models

__obseval.R__ Performs the initial (sampling of multivariate truncated Gaussian from scratch) E step of the SMCEM procedure for multivariate probit models

__PbitEMloop.R__ Performs one SMC EM step of the SMCEM procedure for multivariate probit models, including a call to `smcEM` and to `Mstep`

__PbitEMstart.R__ Performs the first step (sampling from scratch) of the SMC EM procedure for multivariate probit models

__PbitEMstep.R__ Performs the EM step for multivariate probit models, with sampling from scratch at each iterations (only used without the sequential update of the parameters) 

__resample.R__ Performs resampling for the SMC sampler

__siginit.R__ Evaluates an initial pairwise estimate of the covariance matrix

__simProbit.R__ Simulates data from a multivariate probit models

__sixCseqDemo.R__ Performs analysis of six cities data `sixcities.RData'

__smcEM.R__ Performs one step of the sequential SMC EM update procedure for multivariate probit models

__smcemPbitDemo.R__ Performs analysis of simulated data `simPbit8Data1k.RData'
****

****
#### Knitting files, `knitSMCEMpbit` folder
****
__knit6data.Rnw__ Source file (as a child of `knitSMCEMpbit.Rnw`) for knitting the six cities analysis example

__knitSimPbit.Rnw__ Source file (as a child of `knitSMCEMpbit.Rnw`) for knitting the simulated data analysis example

__knitSMCEMpbit.Rnw__ Knits the whole supplement file `knitSMCEMpbit.pdf`

__knitSMCEMpbit.pdf__ Supplement to SMCEMMPM article. It contains detailed procedures, including the code, to fit a probit model to the two datasets considered in the paper.
****

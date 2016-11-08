# --------------------------------------------------------------------------
# stats.R
# 
# Statistics functions script for the coastwide (single stock) 
# simulation-estimation experiments comparing single species and 
# multispecies assessment models.
# 
# Author: Samuel D N Johnson
# Date: 7 September, 2016
#
# --------------------------------------------------------------------------


## THIS ONLY WORKS FOR ADMB OUTPUT
# makeErrorDists()
# This function takes a blob object produced by a sim-est procedure
# and produces distributions of absolute and relative errors
# inputs:   blob=list object output of simEstProc
# ouputs:   blob=list object with error distributions appended
# usage:    to create output that is saved for later analysis
.makeRelErrorDists <- function ( blob )
{
  # recover control list, OM and AMs
  ctrl  <- blob$ctrl
  om    <- blob$om
  opMod <- blob$opMod
  ss    <- blob$am$ss
  ms    <- blob$am$ms

  # Get control constants
  nReps <- ctrl$nReps
  nS    <- opMod$nS
  nT    <- opMod$nT

  # First get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  ssHess <- ss$hesspd
  ssHess <- as.logical(apply ( X = ssHess, FUN = prod, MARGIN = 1 ) )
  msHess <- ms$hesspd
  ssSuccessful <- which ( ssHess )
  msSuccessful <- which ( msHess )
  success <- intersect ( ssSuccessful, msSuccessful )
  failure <- (1:nReps)[-success]

  # Create a list to hold error values
  ssErr <- list ( 
                Bmsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                Umsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                sigma2= matrix ( NA, nrow = nReps, ncol = nS ),
                tau2  = matrix ( NA, nrow = nReps, ncol = nS ),
                q     = matrix ( NA, nrow = nReps, ncol = nS ),
                mlnq  = matrix ( NA, nrow = nReps, ncol = nS ),
                dep   = matrix ( NA, nrow = nReps, ncol = nS ),
                BnT   = matrix ( NA, nrow = nReps, ncol = nS )
              )
  # Slightly difference structure for MS model
  msErr <- list ( 
                Bmsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                Umsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                sigma2= matrix ( NA, nrow=nReps, ncol=1),
                Sigma2= matrix ( NA, nrow = nReps, ncol = nS ),
                tau2  = matrix ( NA, nrow = nReps, ncol = nS ),
                q     = matrix ( NA, nrow = nReps, ncol = nS ),
                mlnq  = matrix ( NA, nrow = nReps, ncol = 1 ),
                dep   = matrix ( NA, nrow = nReps, ncol = nS ),
                BnT   = matrix ( NA, nrow = nReps, ncol = nS )
              )

  # append error lists to blob
  # Single species only has one error term
  ss$err.mle <- ssErr

  # now append Sigma2 to the err list
  ms$err.mle <- msErr

  # Fill in ss MLE relative errors
  for ( s in 1:nS )
  {
    ss$err.mle$Bmsy[,s]   <- (ss$Bmsy[,s] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ss$err.mle$Umsy[,s]   <- (ss$Umsy[,s] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]
    ss$err.mle$kappa2[,s] <- (ss$kappa2[,s] - opMod$pars$kappa2)/opMod$pars$kappa2
    ss$err.mle$tau2[,s]   <- (ss$tau2[,s] - opMod$tau2[s])/opMod$tau2[s]
    ss$err.mle$q[,s]      <- (ss$q[,s] - opMod$q[s])/opMod$q[s]
    ss$err.mle$mlnq[,s]   <- (ss$mlnq[,s] - mean(log(opMod$q)))
    ss$err.mle$dep[,s]    <- (ss$dep[,s] - om$dep[,s])/om$dep[,s]
    ss$err.mle$BnT[,s]    <- (ss$Bt[,s,nT] - om$Bt[,s,nT])/om$Bt[,s,nT]

    # Now fill in ms MLE relative errors
    # some are only estimated once (instead of nS times)
    if (s == 1)
    {
      ms$err.mle$kappa2     <- t(ms$kappa2 - opMod$pars$kappa2)/opMod$pars$sigma2
      ms$err.mle$mlnq       <- t(ms$mlnq - mean(log(opMod$q)))
    }
    # Now the rest of the pars
    ms$err.mle$Bmsy[,s]   <- (ms$Bmsy[,s] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ms$err.mle$Umsy[,s]   <- (ms$Umsy[,s] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]
    ms$err.mle$Sigma2[,s] <- (ms$Sigma2[,s] - opMod$pars$Sigma2[s])/opMod$pars$Sigma2[s]
    ms$err.mle$tau2[,s]   <- (ms$tau2[,s] - opMod$tau2[s])/opMod$tau2[s]
    ms$err.mle$q[,s]      <- (ms$q[,s] - opMod$q[s])/opMod$q[s]
    ms$err.mle$dep[,s]    <- (ms$dep[,s] - om$dep[s])/om$dep[s]
    ms$err.mle$BnT[,s]    <- (ms$Bt[,s,nT] - om$Bt[,s,nT])/om$Bt[,s,nT]
  }

  # Append these to blob
  blob$am$ss <- ss
  blob$am$ms <- ms

  # Now save the good replicates
  blob$goodReps <- success

  blob
}



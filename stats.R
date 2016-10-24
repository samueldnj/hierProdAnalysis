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
makeRelErrorDists <- function ( blob )
{
  # recover control list, OM and AMs
  ctl <- blob$ctl
  om  <- blob$om
  ss  <- blob$am$ss
  ms  <- blob$am$ms

  # Get control constants
  nReps <- ctl$nReps
  nS    <- ctl$nS
  nT    <- ctl$nT

  # First get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  ssHess <- ss$hesspd
  ssHess <- as.logical(apply ( X = ssHess, FUN = prod, MARGIN = 1 ) )
  msHess <- ms$hesspd
  ssSuccessful <- which ( ssHess )
  msSuccessful <- which ( msHess )
  success <- intersect ( ssSuccessful, msSuccessful )
  failure <- (1:nReps)[-success]


  
  # Create a list to hold error values
  err <- list ( 
                Bmsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                Fmsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                sigma2= matrix ( NA, nrow = nReps, ncol = nS ),
                tau2  = matrix ( NA, nrow = nReps, ncol = nS ),
                q     = matrix ( NA, nrow = nReps, ncol = nS ),
                mlnq  = matrix ( NA, nrow = nReps, ncol = nS ),
                slnq  = matrix ( NA, nrow = nReps, ncol = nS ),
                dep   = matrix ( NA, nrow = nReps, ncol = nS ),
                BnT   = matrix ( NA, nrow = nReps, ncol = nS )
              )
# append error lists to blob
ss$err.mle <- err
ms$err.mle <- err



# Fill in ss MLE relative errors
ss$err.mle$Bmsy     <- t( (t(ss$Bmsy) - ctl$Bmsy)/ctl$Bmsy)[success,]
ss$err.mle$Fmsy     <- t( (t(ss$Fmsy) - ctl$Fmsy)/ctl$Fmsy)[success,]
ss$err.mle$sigma2   <- ((ss$sigma2 - ctl$sigma^2)/ctl$sigma^2)[success,]
ss$err.mle$tau2     <- t( (t(ss$tau2) - (ctl$tau)^2)/(ctl$tau)^2)[success,]
ss$err.mle$q        <- t( (t(ss$q) - ctl$q)/ ctl$q )[success,]
ss$err.mle$mlnq     <- t( t(ss$mlnq) - mean(log(ctl$q)))[success,]
# ss$err.mle$slnq   <- t( t(ss$mlnq) - mean(log(ctl$q))) 
ss$err.mle$dep      <- ((ss$dep - om$dep)/om$dep)[success,]
ss$err.mle$BnT      <- ((ss$Bt[,,nT] - om$Bt[,,nT])/om$Bt[,,nT])[success,]

# Now fill in ms MLE relative errors
ms$err.mle$Bmsy     <- t( (t(ms$Bmsy) - ctl$Bmsy)/ctl$Bmsy)[success,]
ms$err.mle$Fmsy     <- t( (t(ms$Fmsy) - ctl$Fmsy)/ctl$Fmsy)[success,]
ms$err.mle$sigma2   <- ((ms$sigma2 - ctl$sigma^2)/ctl$sigma^2)[success,]
ms$err.mle$tau2     <- t( (t(ms$tau2) - (ctl$tau)^2)/(ctl$tau)^2)[success,]
ms$err.mle$q        <- t( (t(ms$q) - ctl$q)/ ctl$q )[success,]
ms$err.mle$mlnq     <- as.matrix(t( t(ms$mlnq) - mean(log(ctl$q)))[success,])
# ms$err.mle$slnq   <- t( t(ms$mlnq) - mean(log(ctl$q))) 
ms$err.mle$dep      <- ((ms$dep - om$dep)/om$dep)[success,]
ms$err.mle$BnT      <- ((ms$Bt[,,nT] - om$Bt[,,nT])/om$Bt[,,nT])[success,]

# Append these to blob
blob$am$ss <- ss
blob$am$ms <- ms

# Now save the good replicates
blob$goodReps <- success

blob
}



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

#.statTable()
# Wrapper for .simStats, produces stacked tables of stats for a group
# of simulations.
# inputs:   sims=integer vector indicating simulations in ./project/
# outputs:  statTable=table of statistics for a project/group
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
.statTableMSE <- function (sims=1,tabName = "statTable.csv")
{ 
  # call function
  tableList <- lapply ( X = sims, FUN = .simStatMSE )

  # now make the table and return
  statTable <-  do.call("rbind",tableList)
  savePath <- file.path(getwd(),"project","Statistics",tabName)
  write.csv ( statTable, file = savePath )
  statTable
}

# .simStatMSE()
# Produces a statistics table for leading pars in a simulation from
# the output produced by a runSimEst() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame of mean squared error BnT and Umsy
# usage:    in lapply to produce stats for a group of simulations
.simStatMSE <- function ( sim=1 )
{
  # First, load blob
  .loadSim(sim)
  om    <- blob$om
  opMod <- blob$opMod
  pars  <- blob$opMod$pars
  ss    <- blob$am$ss
  ms    <- blob$am$ms
  
  # Control info
  nS      <- blob$opMod$nS
  nT      <- blob$opMod$nT
  species <- blob$ctrl$speciesName
  nReps   <- blob$ctrl$nReps

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","species","kappa2True",
                  "Sigma2True", "corrMult","ssBnT","msBnT","ssUmsy","msUmsy",
                  "ssBmsy","msBmsy","ssMSY","msMSY","ssDep","msDep",
                  "ssq","msq", "msHessPD", "ssHessPD" )
  
  statTable <- matrix(NA,nrow=nS,ncol=length(colLabels))
  
  colnames(statTable)   <- colLabels
  statTable             <- as.data.frame(statTable)

  # get multiplier for shared effects
  kappaMult <- opMod$kappaMult
  if (is.null(opMod$kappaMult)) kappaMult <- 1

  # Start filling stat table
  # First, info and true pars
  statTable$scenario    <- blob$ctrl$scenarioName
  statTable$mp          <- blob$ctrl$mpLabel
  statTable$species     <- species
  statTable$kappa2True  <- opMod$pars$kappa2*kappaMult
  statTable$Sigma2True  <- opMod$pars$Sigma2
  statTable$corrMult    <- opMod$corrMult

  # Now errors
  for (s in 1:nS)
  {
    statTable[s,"ssBnT"]    <- mean ( (ss$Bt[,s,nT] - om$Bt[,s,nT])^2)
    statTable[s,"msBnT"]    <- mean ( (ms$Bt[,s,nT] - om$Bt[,s,nT])^2)
    statTable[s,"ssUmsy"]   <- mean ( (ss$Umsy[,s] - pars$Umsy[s])^2)
    statTable[s,"msUmsy"]   <- mean ( (ms$Umsy[,s] - pars$Umsy[s])^2)
    statTable[s,"ssBmsy"]   <- mean ( (ss$Bmsy[,s] - pars$Bmsy[s])^2)
    statTable[s,"msBmsy"]   <- mean ( (ms$Bmsy[,s] - pars$Bmsy[s])^2)
    statTable[s,"ssMSY"]    <- mean ( (ss$msy[,s] -  pars$Umsy[s]*pars$Bmsy[s])^2)
    statTable[s,"msMSY"]    <- mean ( (ms$msy[,s] -  pars$Umsy[s]*pars$Bmsy[s])^2)
    statTable[s,"ssDep"]    <- mean ( (ss$dep[,s] -  om$dep[,s])^2)
    statTable[s,"msDep"]    <- mean ( (ms$dep[,s] -  om$dep[,s])^2)
    statTable[s,"ssq"]      <- mean ( (ss$q[,s]   -  opMod$q[s])^2)
    statTable[s,"msq"]      <- mean ( (ms$q[,s]   -  opMod$q[s])^2)
    statTable[s,"ssHessPD"] <- mean ( ss$hesspd[,s] )
    statTable[s,"msHessPD"] <- mean ( ms$hesspd )
  }
  
  # return
  statTable
}


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
                kappa2= matrix ( NA, nrow = nReps, ncol = nS ),
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
                kappa2= matrix ( NA, nrow = nReps, ncol=1),
                Sigma2= matrix ( NA, nrow = nReps, ncol = nS ),
                tau2  = matrix ( NA, nrow = nReps, ncol = 1 ),
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
    ss$err.mle$kappa2[,s] <- (ss$kappa2[,s] - (opMod$pars$kappa2+opMod$pars$Sigma2[s]))/(opMod$pars$kappa2+opMod$pars$Sigma2[s])
    ss$err.mle$tau2[,s]   <- (ss$tau2[,s] - opMod$tau2[s])/opMod$tau2[s]
    ss$err.mle$q[,s]      <- (ss$q[,s] - opMod$q[s])/opMod$q[s]
    ss$err.mle$mlnq[,s]   <- (ss$mlnq[,s] - mean(log(opMod$q)))
    ss$err.mle$dep[,s]    <- (ss$dep[,s] - om$dep[,s])/om$dep[,s]
    ss$err.mle$BnT[,s]    <- (ss$Bt[,s,nT] - om$Bt[,s,nT])/om$Bt[,s,nT]

    # Now fill in ms MLE relative errors
    # some are only estimated once (instead of nS times)
    if (s == 1)
    {
      ms$err.mle$kappa2     <- t(ms$kappa2 - opMod$pars$kappa2)/opMod$pars$kappa2
      ms$err.mle$mlnq       <- t(ms$mlnq - mean(log(opMod$q)))
      ms$err.mle$tau2       <- (ms$tau2 - mean(opMod$tau2))/mean(opMod$tau2)
    }
    # Now the rest of the pars
    ms$err.mle$Bmsy[,s]   <- (ms$Bmsy[,s] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ms$err.mle$Umsy[,s]   <- (ms$Umsy[,s] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]
    ms$err.mle$Sigma2[,s] <- (ms$Sigma2[,s] - opMod$pars$Sigma2[s])/opMod$pars$Sigma2[s]
    
    ms$err.mle$q[,s]      <- (ms$q[,s] - opMod$q[s])/opMod$q[s]
    ms$err.mle$dep[,s]    <- (ms$dep[,s] - om$dep[s])/om$dep[s]
    ms$err.mle$BnT[,s]    <- (ms$Bt[,s,nT] - om$Bt[,s,nT])/om$Bt[,s,nT]
  }

  # Now get errors from posterior estimates
  ss$err.post <- ssErr
  ms$err.post <- msErr
  # first calculate SS posterior mean and median
  ssMCMC <- ss$mcPar
  ssMCMC <- apply ( X = ssMCMC, FUN = mean, MARGIN = c(1,2,4))
  # then MS
  msMCMC <- ms$mcPar
  msMCMC <- apply ( X = msMCMC, FUN = mean, MARGIN = c(1,2,4))

  # Then load into the blob error distributions. These contain
  # NAs, which will have to be cleared up in the plotting
  for (s in 1:nS)
  {
    ss$err.post$Bmsy[,s]   <- (ssMCMC[,s,"Bmsy"] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ss$err.post$Umsy[,s]   <- (ssMCMC[,s,"Umsy"] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]
    ss$err.post$kappa2[,s] <- (ssMCMC[,s,"kappa2"] - (opMod$pars$kappa2+opMod$pars$Sigma2[s]))/(opMod$pars$kappa2+opMod$pars$Sigma2[s])
    ss$err.post$tau2[,s]   <- (ssMCMC[,s,"tau2"] - opMod$tau2[s])/opMod$tau2[s]
    ss$err.post$q[,s]      <- (ssMCMC[,s,"q"] - opMod$q[s])/opMod$q[s]
    # ss$err.post$mlnq[,s]   <- (ss$mlnq[,s] - mean(log(opMod$q)))
    ss$err.post$dep[,s]    <- (ssMCMC[,s,"dep_bar"] - om$dep[,s])/om$dep[,s]
    ss$err.post$BnT[,s]    <- (ssMCMC[,s,"BnT"] - om$Bt[,s,nT])/om$Bt[,s,nT]

    # Now fill in ms MLE relative errors
    # some are only estimated once (instead of nS times)
    if (s == 1)
    {
      ms$err.post$kappa2     <- (msMCMC[,s,"kappa2"] - opMod$pars$kappa2)/opMod$pars$kappa2
      ms$err.post$Sigma2     <- (msMCMC[,s,"Sigma2"] - mean(opMod$pars$Sigma2))/mean(opMod$pars$Sigma2)
      ms$err.post$tau2       <- (msMCMC[,s,"tau2"] - mean(opMod$tau2))/mean(opMod$tau2)
      # ms$err.post$mlnq       <- t(ms$mlnq - mean(log(opMod$q)))
    }

    ms$err.post$Bmsy[,s]   <- (msMCMC[,s,"Bmsy"] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ms$err.post$Umsy[,s]   <- (msMCMC[,s,"Umsy"] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]
    ms$err.post$q[,s]      <- (msMCMC[,s,"q"] - opMod$q[s])/opMod$q[s]
    # ms$err.post$mlnq[,s]   <- (ms$mlnq[,s] - mean(log(opMod$q)))
    ms$err.post$dep[,s]    <- (msMCMC[,s,"dep_bar"] - om$dep[,s])/om$dep[,s]
    ms$err.post$BnT[,s]    <- (msMCMC[,s,"BnT"] - om$Bt[,s,nT])/om$Bt[,s,nT]
  }


  # Append these to blob
  blob$am$ss <- ss
  blob$am$ms <- ms

  # Now save the good replicates
  blob$goodReps <- success

  blob
}



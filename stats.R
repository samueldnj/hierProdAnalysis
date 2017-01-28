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

  # get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  success <- blob$goodReps

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","species","kappaTrue",
                  "SigmaTrue", "kappaMult", "corrMult","ssBnT","msBnT","ssUmsy","msUmsy",
                  "ssBmsy","msBmsy","ssMSY","msMSY","ssDep","msDep",
                  "ssq","msq", "msHessPD", "ssHessPD", "nReps",
                  "Umax", "tUpeak", "tUtrough" )
  
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
  statTable$kappaTrue   <- sqrt(opMod$pars$kappa2)*kappaMult
  statTable$SigmaTrue   <- sqrt(opMod$pars$Sigma2)
  statTable$kappaMult   <- kappaMult
  statTable$corrMult    <- opMod$corrMult
  statTable$nReps       <- nReps
  statTable$Umax        <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUtrough    <- opMod$tUtrough
  statTable$tUpeak      <- opMod$tUpeak

  # Now errors
  for (s in 1:nS)
  {
    statTable[s,"ssBnT"]    <- mean ( (ss$Bt[success,s,nT] - om$Bt[success,s,nT])^2, na.rm=TRUE)
    statTable[s,"msBnT"]    <- mean ( (ms$Bt[success,s,nT] - om$Bt[success,s,nT])^2, na.rm=TRUE)
    statTable[s,"ssUmsy"]   <- mean ( (ss$Umsy[success,s] - pars$Umsy[s])^2, na.rm=TRUE)
    statTable[s,"msUmsy"]   <- mean ( (ms$Umsy[success,s] - pars$Umsy[s])^2, na.rm=TRUE)
    statTable[s,"ssBmsy"]   <- mean ( (ss$Bmsy[success,s] - pars$Bmsy[s])^2, na.rm=TRUE)
    statTable[s,"msBmsy"]   <- mean ( (ms$Bmsy[success,s] - pars$Bmsy[s])^2, na.rm=TRUE)
    statTable[s,"ssMSY"]    <- mean ( (ss$msy[success,s] -  pars$Umsy[s]*pars$Bmsy[s])^2, na.rm=TRUE)
    statTable[s,"msMSY"]    <- mean ( (ms$msy[success,s] -  pars$Umsy[s]*pars$Bmsy[s])^2, na.rm=TRUE)
    statTable[s,"ssDep"]    <- mean ( (ss$dep[success,s,nT] -  om$dep[success,s,nT])^2, na.rm=TRUE)
    statTable[s,"msDep"]    <- mean ( (ms$dep[success,s,nT] -  om$dep[success,s,nT])^2, na.rm=TRUE)
    statTable[s,"ssq"]      <- mean ( (ss$q[success,s]   -  opMod$q[s])^2, na.rm=TRUE)
    statTable[s,"msq"]      <- mean ( (ms$q[success,s]   -  opMod$q[s])^2, na.rm=TRUE)
    statTable[s,"ssHessPD"] <- sum ( ss$hesspd[,s] , na.rm=TRUE)
    statTable[s,"msHessPD"] <- sum ( ms$hesspd , na.rm=TRUE)
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
                dep   = matrix ( NA, nrow = nReps, ncol = nS ),
                BnT   = matrix ( NA, nrow = nReps, ncol = nS ),
                totRE = matrix ( NA, nrow = nReps, ncol = nS )
              )
  # Slightly difference structure for MS model
  msErr <- list ( 
                Bmsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                Umsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                kappa2= matrix ( NA, nrow = nReps, ncol = 1 ),
                Sigma2= matrix ( NA, nrow = nReps, ncol = nS ),
                tau2  = matrix ( NA, nrow = nReps, ncol = nS ),
                q     = matrix ( NA, nrow = nReps, ncol = nS ),
                dep   = matrix ( NA, nrow = nReps, ncol = nS ),
                BnT   = matrix ( NA, nrow = nReps, ncol = nS ),
                totRE = matrix ( NA, nrow = nReps, ncol = nS )
              )

  # append error lists to blob
  # Single species only has one error term
  ss$err.mle <- ssErr

  # now append Sigma2 to the err list
  ms$err.mle <- msErr

  # Fill in ss MLE relative errors
  for ( s in 1:nS )
  {
    ss$err.mle$Bmsy[success,s]   <- (ss$Bmsy[success,s] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ss$err.mle$Umsy[success,s]   <- (ss$Umsy[success,s] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]
    ss$err.mle$kappa2[success,s] <- (ss$kappa2[success,s] - (opMod$pars$kappa2+opMod$pars$Sigma2[s]))/(opMod$pars$kappa2+opMod$pars$Sigma2[s])
    ss$err.mle$tau2[success,s]   <- (ss$tau2[success,s] - opMod$tau2[s])/opMod$tau2[s]
    ss$err.mle$q[success,s]      <- (ss$q[success,s] - opMod$q[s])/opMod$q[s]
    ss$err.mle$dep[success,s]    <- (ss$dep[success,s,nT] - om$dep[success,s,nT])/om$dep[success,s,nT]
    ss$err.mle$BnT[success,s]    <- (ss$Bt[success,s,nT] - om$Bt[success,s,nT])/om$Bt[success,s,nT]
    ss$err.mle$totRE[success,s]  <- ss$err.mle$kappa2[success,s]

    # Now fill in ms MLE relative errors
    # some are only estimated once (instead of nS times)
    if (s == 1)
    {
      ms$err.mle$kappa2[success,]    <- (ms$kappa2[success,] - opMod$pars$kappa2)/opMod$pars$kappa2
    }
    # Now the rest of the pars
    ms$err.mle$Sigma2[success,s] <- (ms$Sigma2[success,s] - opMod$pars$Sigma2[s])/opMod$pars$Sigma2[s]
    ms$err.mle$tau2[success,s]   <- (ms$tau2[success,s] - opMod$tau2[s])/opMod$tau2[s]
    ms$err.mle$Bmsy[success,s]   <- (ms$Bmsy[success,s] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ms$err.mle$Umsy[success,s]   <- (ms$Umsy[success,s] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]    
    ms$err.mle$q[success,s]      <- (ms$q[success,s] - opMod$q[s])/opMod$q[s]
    ms$err.mle$dep[success,s]    <- (ms$dep[success,s,nT] - om$dep[success,s,nT])/om$dep[success,s,nT]
    ms$err.mle$BnT[success,s]    <- (ms$Bt[success,s,nT] - om$Bt[success,s,nT])/om$Bt[success,s,nT]

    # calculate total RE variance
    totVarFit <- ms$kappa2 + ms$Sigma2[,s]
    totVarOM  <- opMod$pars$kappa2*(opMod$kappaMult^2) + opMod$pars$Sigma2[s]
    ms$err.mle$totRE[success,s]  <- (totVarFit[success] - totVarOM) / totVarOM
  }

  # Append these to blob
  blob$am$ss <- ss
  blob$am$ms <- ms

  # Now save the good replicates
  blob$goodReps <- success

  blob
}



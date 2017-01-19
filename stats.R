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
.statTableMSE <- function( sims=1, tabName = "statTable.csv" )
{ 
  # call function
  tableList <- lapply ( X = sims, FUN = .simStatMSE )

  # now make the table and return
  statTable <-  do.call( "rbind", tableList )
  savePath  <- file.path(getwd(),"project","Statistics",tabName)
  write.csv( statTable, file = savePath )
  statTable
}

# .simStatMSE()
# Produces a statistics table for leading pars in a simulation from
# the output produced by a runSimEst() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame of mean squared error BnT and Umsy
# usage:    in lapply to produce stats for a group of simulations
.simStatMSE <- function ( sim=1, est="MCMC" )
{
  # First, load blob
  .loadSim(sim)
  om    <- blob$om
  opMod <- blob$opMod
  pars  <- blob$opMod$pars

  # Control info
  nS      <- blob$opMod$nS
  nT      <- blob$opMod$nT
  species <- blob$ctrl$speciesName
  nReps   <- blob$ctrl$nReps
  
  # Get estimates, whether MCMC or MLE
  ss    <- blob$am$ss
  ms    <- blob$am$ms

  # Get the good replicate numbers
  ssHess <- ss$hesspd
  ssHess[is.na(ssHess)] <- FALSE
  ssHess <- as.logical(apply ( X = ssHess, FUN = prod, MARGIN = 1 ) )
  msHess <- ms$hesspd
  msHess[is.na(msHess)] <- FALSE
  ssSuccessful <- which ( ssHess )
  msSuccessful <- which ( msHess )
  success <- intersect ( ssSuccessful, msSuccessful )
  failure <- (1:nReps)[-success]


  if( est == "MCMC" )
  {
    # Set quantiles (add to control file)
    qProbs  <- c(0.025,0.5,0.975)

    # Save MCMC output
    ssMCMC  <- blob$am$ss$mcOut
    msMCMC  <- blob$am$ms$mcOut

    # Apply functions to extract quantiles and posterior means
    ssPM    <- apply ( X = ssMCMC, FUN = mean, MARGIN = c(1,2,4), na.rm=T)
    ssQuant <- apply ( X = ssMCMC, FUN = quantile, MARGIN = c(1,2,4), probs = qProbs, na.rm = T )
    msPM    <- apply ( X = msMCMC, FUN = mean, MARGIN = c(1,3), na.rm=T)
    msQuant <- apply ( X = msMCMC, FUN = quantile, MARGIN = c(1,3), probs = qProbs, na.rm=T )
  }

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
    if ( est == "MLE")
    {
      statTable[s,"ssBnT"]    <- mean ( (ss$Bt[,s,nT] - om$Bt[,s,nT])^2,na.rm=TRUE)
      statTable[s,"msBnT"]    <- mean ( (ms$Bt[,s,nT] - om$Bt[,s,nT])^2,na.rm=TRUE)
      statTable[s,"ssUmsy"]   <- mean ( (ss$Umsy[,s] - pars$Umsy[s])^2,na.rm=TRUE)
      statTable[s,"msUmsy"]   <- mean ( (ms$Umsy[,s] - pars$Umsy[s])^2,na.rm=TRUE)
      statTable[s,"ssBmsy"]   <- mean ( (ss$Bmsy[,s] - pars$Bmsy[s])^2,na.rm=TRUE)
      statTable[s,"msBmsy"]   <- mean ( (ms$Bmsy[,s] - pars$Bmsy[s])^2,na.rm=TRUE)
      statTable[s,"ssMSY"]    <- mean ( (ss$msy[,s] -  pars$Umsy[s]*pars$Bmsy[s])^2,na.rm=TRUE)
      statTable[s,"msMSY"]    <- mean ( (ms$msy[,s] -  pars$Umsy[s]*pars$Bmsy[s])^2,na.rm=TRUE)
      statTable[s,"ssDep"]    <- mean ( (ss$dep[,s] -  om$dep[,s])^2,na.rm=TRUE)
      statTable[s,"msDep"]    <- mean ( (ms$dep[,s] -  om$dep[,s])^2,na.rm=TRUE)
      statTable[s,"ssq"]      <- mean ( (ss$q[,s]   -  opMod$q[s])^2,na.rm=TRUE)
      statTable[s,"msq"]      <- mean ( (ms$q[,s]   -  opMod$q[s])^2,na.rm=TRUE)
      statTable[s,"ssHessPD"] <- mean ( ss$hesspd[,s] )
      statTable[s,"msHessPD"] <- mean ( ms$hesspd )  
    }
    if( est == "MCMC" )
    {
      # Create labels to recover MCMC table entries
      ssBmsy    <- paste("Bmsy",1,sep="")
      ssUmsy    <- paste("Umsy",1,sep="")
      ssMSY     <- paste("msy",1,sep="")
      ssqlab    <- paste("q",1,sep="")
      ssDnTlab  <- paste("DnT",1,sep="")
      ssBnTlab  <- paste("Bst_",1,"_",nT,sep="")

      msBmsy    <- paste("Bmsy",s,sep="")
      msUmsy    <- paste("Umsy",s,sep="")
      msMSY     <- paste("msy",s,sep="")
      msqlab    <- paste("q",s,sep="")
      msDnTlab  <- paste("DnT",s,sep="")
      msBnTlab  <- paste("Bst_",s,"_",nT,sep="")


      statTable[s,"ssBnT"]    <- mean ( (ssPM[success,s,ssBnTlab] - om$Bt[success,s,nT])^2,na.rm=TRUE)
      statTable[s,"msBnT"]    <- mean ( (msPM[success,msBnTlab] - om$Bt[success,s,nT])^2,na.rm=TRUE)
      statTable[s,"ssUmsy"]   <- mean ( (ssPM[success,s,ssUmsy] - pars$Umsy[s])^2,na.rm=TRUE)
      statTable[s,"msUmsy"]   <- mean ( (msPM[success,msUmsy] - pars$Umsy[s])^2,na.rm=TRUE)
      statTable[s,"ssBmsy"]   <- mean ( (ssPM[success,s,ssBmsy] - pars$Bmsy[s])^2,na.rm=TRUE)
      statTable[s,"msBmsy"]   <- mean ( (msPM[success,msBmsy] - pars$Bmsy[s])^2,na.rm=TRUE)
      statTable[s,"ssMSY"]    <- mean ( (ssPM[success,s,ssMSY] -  pars$Umsy[s]*pars$Bmsy[s])^2,na.rm=TRUE)
      statTable[s,"msMSY"]    <- mean ( (msPM[success,msMSY] -  pars$Umsy[s]*pars$Bmsy[s])^2,na.rm=TRUE)
      statTable[s,"ssDep"]    <- mean ( (ssPM[success,s,ssDnTlab] -  om$dep[,s])^2,na.rm=TRUE)
      statTable[s,"msDep"]    <- mean ( (msPM[success,msDnTlab] -  om$dep[,s])^2,na.rm=TRUE)
      statTable[s,"ssq"]      <- mean ( (ssPM[success,s,ssqlab]   -  opMod$q[s])^2,na.rm=TRUE)
      statTable[s,"msq"]      <- mean ( (msPM[success,msqlab]   -  opMod$q[s])^2,na.rm=TRUE)
      statTable[s,"ssHessPD"] <- mean ( ss$hesspd[,s], na.rm=TRUE )
      statTable[s,"msHessPD"] <- mean ( ms$hesspd, na.rm=TRUE )  
    }
    
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
  ssHess[is.na(ssHess)] <- FALSE
  ssHess <- as.logical(apply ( X = ssHess, FUN = prod, MARGIN = 1 ) )
  msHess <- ms$hesspd
  msHess[is.na(msHess)] <- FALSE
  ssSuccessful <- which ( ssHess )
  msSuccessful <- which ( msHess )
  success <- intersect ( ssSuccessful, msSuccessful )
  failure <- (1:nReps)[-success]

  # Create a list to hold error values
  ssErr <- list ( 
                Bmsy   = matrix ( NA, nrow = nReps, ncol = nS ),
                Umsy   = matrix ( NA, nrow = nReps, ncol = nS ),
                kappa2 = matrix ( NA, nrow = nReps, ncol = nS ),
                totVar = matrix ( NA, nrow = nReps, ncol = nS ),
                tau2   = matrix ( NA, nrow = nReps, ncol = nS ),
                q      = matrix ( NA, nrow = nReps, ncol = nS ),
                # mlnq   = matrix ( NA, nrow = nReps, ncol = nS ),
                dep    = matrix ( NA, nrow = nReps, ncol = nS ),
                BnT    = matrix ( NA, nrow = nReps, ncol = nS )
              )
  # Slightly difference structure for MS model
  msErr <- list ( 
                Bmsy   = matrix ( NA, nrow = nReps, ncol = nS ),
                Umsy   = matrix ( NA, nrow = nReps, ncol = nS ),
                kappa2 = matrix ( NA, nrow = nReps, ncol = 1 ),
                Sigma2 = matrix ( NA, nrow = nReps, ncol = nS ),
                totVar = matrix ( NA, nrow = nReps, ncol = nS ),
                tau2   = matrix ( NA, nrow = nReps, ncol = nS ),
                q      = matrix ( NA, nrow = nReps, ncol = nS ),
                # mlnq   = matrix ( NA, nrow = nReps, ncol = 1 ),
                dep    = matrix ( NA, nrow = nReps, ncol = nS ),
                BnT    = matrix ( NA, nrow = nReps, ncol = nS )
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
    ss$err.mle$totVar[success,s] <- (ss$kappa2[success,s] - (opMod$pars$kappa2+opMod$pars$Sigma2[s]))/(opMod$pars$kappa2+opMod$pars$Sigma2[s])
    ss$err.mle$tau2[success,s]   <- (ss$tau2[success,s] - opMod$tau2[s])/opMod$tau2[s]
    ss$err.mle$q[success,s]      <- (ss$q[success,s] - opMod$q[s])/opMod$q[s]
    ss$err.mle$dep[success,s]    <- (ss$dep[success,s] - om$dep[success,s])/om$dep[success,s]
    ss$err.mle$BnT[success,s]    <- (ss$Bt[success,s,nT] - om$Bt[success,s,nT])/om$Bt[success,s,nT]

    # Now fill in ms MLE relative errors
    # some are only estimated once (instead of nS times)
    if (s == 1)
    {
      ms$err.mle$kappa2[success,s]    <- (ms$kappa2[success,] - opMod$pars$kappa2)/opMod$pars$kappa2      
    }
    # Now the rest of the pars
    ms$err.mle$tau2[success,s]   <- (ms$tau2[success,] - opMod$tau2[s])/opMod$tau2[s]
    ms$err.mle$Sigma2[success,s] <- (ms$Sigma2[success,] - opMod$pars$Sigma2[s])/opMod$pars$Sigma2[s]  
    ms$err.mle$totVar[success,s] <- (ms$Sigma2[success,] + ms$kappa2[success,] - (opMod$pars$Sigma2[s]+opMod$pars$kappa2))/(opMod$pars$Sigma2[s]+opMod$pars$kappa2)  
    ms$err.mle$Bmsy[success,s]   <- (ms$Bmsy[success,s] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ms$err.mle$Umsy[success,s]   <- (ms$Umsy[success,s] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]
    ms$err.mle$q[success,s]      <- (ms$q[success,s] - opMod$q[s])/opMod$q[s]
    ms$err.mle$dep[success,s]    <- (ms$dep[success,s] - om$dep[s])/om$dep[s]
    ms$err.mle$BnT[success,s]    <- (ms$Bt[success,s,nT] - om$Bt[success,s,nT])/om$Bt[success,s,nT]
  }

  # Now get errors from posterior estimates
  ss$err.post <- ssErr
  ms$err.post <- msErr
  # Recover MCMC chains
  ssMCMC <- blob$am$ss$mcOut
  msMCMC <- blob$am$ms$mcOut
  # Calculate posterior means
  ssPostMean <- apply ( X = ssMCMC, FUN = mean, MARGIN = c(1,2,4), na.rm=T)
  msPostMean <- apply ( X = msMCMC, FUN = mean, MARGIN = c(1,3), na.rm=T)

  # Now save relative error distributions when restricting to successful
  # runs. These contain NAs that will have to be cleared up in the plotting
  for (s in 1:nS)
  {
    # Create labels to recover parameter columns from MCMC tables
    ssBmsy    <- paste("Bmsy",1,sep="")
    ssUmsy    <- paste("Umsy",1,sep="")
    ssqlab    <- paste("q",1,sep="")
    ssDnTlab  <- paste("DnT",1,sep="")
    ssBnTlab  <- paste("Bst_",1,"_",nT,sep="")

    msBmsy    <- paste("Bmsy",s,sep="")
    msUmsy    <- paste("Umsy",s,sep="")
    msqlab    <- paste("q",s,sep="")
    msDnTlab  <- paste("DnT",s,sep="")
    msBnTlab  <- paste("Bst_",s,"_",nT,sep="")
    # Now save SS estimates
    ss$err.post$Bmsy[success,s]   <- (ssPostMean[success,s,ssBmsy] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ss$err.post$Umsy[success,s]   <- (ssPostMean[success,s,ssUmsy] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]
    ss$err.post$kappa2[success,s] <- (ssPostMean[success,s,"kappa2"] - (opMod$pars$kappa2+opMod$pars$Sigma2[s]))/(opMod$pars$kappa2+opMod$pars$Sigma2[s])
    ss$err.post$totVar[success,s] <- (ssPostMean[success,s,"kappa2"] - (opMod$pars$kappa2+opMod$pars$Sigma2[s]))/(opMod$pars$kappa2+opMod$pars$Sigma2[s])
    ss$err.post$tau2[success,s]   <- (ssPostMean[success,s,"tau2"] - opMod$tau2[s])/opMod$tau2[s]
    ss$err.post$q[success,s]      <- (ssPostMean[success,s,ssqlab] - opMod$q[s])/opMod$q[s]
    ss$err.post$dep[success,s]    <- (ssPostMean[success,s,ssDnTlab] - om$dep[success,s])/om$dep[success,s]
    ss$err.post$BnT[success,s]    <- (ssPostMean[success,s,ssBnTlab] - om$Bt[success,s,nT])/om$Bt[success,s,nT]
    # Now fill in ms MLE relative errors
    # some are only estimated once (instead of nS times)
    if (s == 1)
    {
      ms$err.post$kappa2[success,s]     <- (msPostMean[success,"kappa2"] - opMod$pars$kappa2)/opMod$pars$kappa2
    }

    ms$err.post$Sigma2[success,s] <- (msPostMean[success,"Sigma2"] - opMod$pars$Sigma2[s])/opMod$pars$Sigma2[s]
    ms$err.post$totVar[success,s] <- (msPostMean[success,"Sigma2"] + msPostMean[success,"kappa2"] - (opMod$pars$Sigma2[s]+opMod$pars$kappa2))/(opMod$pars$Sigma2[s]+opMod$pars$kappa2)
    ms$err.post$tau2[success,s]   <- (msPostMean[success,"tau2"] - opMod$tau2[s])/opMod$tau2[s]
    ms$err.post$Bmsy[success,s]   <- (msPostMean[success,msBmsy] - opMod$pars$Bmsy[s])/opMod$pars$Bmsy[s]
    ms$err.post$Umsy[success,s]   <- (msPostMean[success,msUmsy] - opMod$pars$Umsy[s])/opMod$pars$Umsy[s]
    ms$err.post$q[success,s]      <- (msPostMean[success,msqlab] - opMod$q[s])/opMod$q[s]
    ms$err.post$dep[success,s]    <- (msPostMean[success,msDnTlab] - om$dep[success,s])/om$dep[success,s]
    ms$err.post$BnT[success,s]    <- (msPostMean[success,msBnTlab] - om$Bt[success,s,nT])/om$Bt[success,s,nT]
  }

  # Append these to blob
  blob$am$ss <- ss
  blob$am$ms <- ms

  # Now save the good replicates
  blob$goodReps <- success

  blob
}



# --------------------------------------------------------------------------
# simulation.R
# 
# Script for the simulation procedure to produce data for ADMB and TMB
# estimators.
# 
# Author: Samuel Johnson
# Date: 24 August, 2016
# 
# 
# --------------------------------------------------------------------------

# opModel()
# The operating model takes control list parameters and simulates biological 
# and observational data for each species, creating an om list object that 
# can then be used to create dat and pin files to be passed to estimators.
# inputs:     control=control list object, read in from control file
#             seed=random number seed
# ouputs:     om=list object containing the data for each species
# usage:      to prepare data for each replicate (seed value)
opModel <- function (control = ctlList, seed = 1)
{
  # First, set the seed
  if ( !is.null(seed) ) set.seed(seed)

  # Now recover multi-species parameters
  nT <- control$nT          # length of simulation
  nS <- control$nS          # number of species

  # Process error component sds
  sigma <- control$sigma    # longitudinal proc error sd
  # Include some kind of auto-regressive parameter here, or 
  # complete cov mtx for epst
  Sigma <- control$Sigma    # interspecies proc error cov mtx

  # Obs error sd
  tau <- control$tau

  # Initialise list to hold the data
  om <- list (  Bt    = matrix (NA,nrow=nS, ncol=nT ),
                Ct    = matrix (NA,nrow=nS, ncol=nT ),
                Ft    = matrix (NA,nrow=nS, ncol=nT ),
                ItTrue= matrix (NA,nrow=nS, ncol=nT ),
                epst  = matrix (NA,nrow=nS, ncol=nT ),
                zetat = matrix (NA,nrow=nS, ncol=nT ),
                deltat= matrix (NA,nrow=nS, ncol=nT ),
                It    = matrix (NA,nrow=nS, ncol=nT ),
                dep   = numeric(length = nS),
                sigma = sigma,
                Sigma = Sigma,
                tau   = tau,
                q     = control$q,
                nT    = nT,
                nS    = nS,
                msy   = control$msy,
                Fmsy  = control$Fmsy )

  # Now create epst and zetat vectors using the proc error components
  #### No (auto) correlation yet ####
  epst      <- matrix(rnorm ( n = nS*nT ), nrow = nS, ncol = nT)
  zetat     <- matrix(rnorm ( n = nS*nT ), nrow = nS, ncol = nT)

  ## Here, we (will) generate the (auto) correlated standard normal
  ## deviations for epst and zetat, then pass them to the 
  ## production model, where they're scaled and bias is removed
  ## Do we need a fillRanWalk approach??? Parametric auto-corr?
  
  # standard normal deviations for obs error (uncorrelated)
  deltat    <- matrix(rnorm ( n = nS*nT ), nrow = nS, ncol = nT) 

  # Create Ft time series
  ## REPLACE WITH A FUNCTION LATER FOR MULTIPLE Ft TRAJECTORIES ##
  FtProp <- numeric ( length = nT )
  Fmult <- control $ F_mult
  # Then get time of peak fishing mortality
  tFpeak <- control $ tFpeak
  # Compute gradients of the pieces
  Fgrad1 <- ( Fmult [ 2 ] - Fmult [ 1 ] ) / (tFpeak - 1 )
  Fgrad2 <- ( Fmult [ 3 ] - Fmult [ 2 ] ) / 
                  ( nT - tFpeak )
  # Populate Ft with multipliers of F_MSY
  FtProp [ 1:tFpeak ]         <- Fmult [ 1 ] + 0:(tFpeak - 1) * Fgrad1
  FtProp [ (tFpeak + 1):nT ]  <- Fmult [ 2 ] + 1:(nT - tFpeak) * Fgrad2

  # Overwrite Ft with actual F values
  Fst <- matrix ( FtProp, nrow = nS, ncol = nT, byrow = TRUE )
  for ( s in 1:nS ) Fst[s,] <- control$Fmsy[s] * Fst[s,]

  # Loop over species and fill list entries with biological
  # and observational data
  for ( s in 1:nS )
  {
    bio <- logProdModel ( msy = control$msy[s], Fmsy = control$Fmsy[s],
                          nT = nT, Ft = Fst[s,], 
                          epst = epst[s,], sigma = sigma, 
                          zetat = zetat[s,], Sigma = diag(Sigma)[1] )
    obs <- obsModel ( Bt = bio $ Bt, q = control$q[s], nT = nT,
                      deltat = deltat[s,], tau = tau[s] )

    om$Bt[s,]     <- bio$Bt
    om$Ct[s,]     <- bio$Ct
    om$Ft[s,]     <- bio$Ft
    om$ItTrue[s,] <- control$q[s]*bio$Bt
    om$It[s,]     <- obs$It
    om$epst[s,]   <- bio$epst
    om$zetat[s,]  <- bio$zetat
    om$deltat[s,] <- obs$deltat
    om$dep[s]     <- bio$dep
  }

  return(om)
}

# makeDataLists()
# The makeDataLists function takes an om list object produced by opModel()
# and creates data structures that will be passed to estimators
# inputs:   omList=operating model list object produced by opModel()
# outputs:  msDat=multispecies data list; 
#           msPar=multispecies initial parameter value list
#           ssDat=nS-list of single species data lists
#           ssPar=nS-list of single species init. par lists
makeDataLists <- function ( omList = om, ctl = ctlList )
{
  # Needs to create 3 SS dat and par lists, and 1 MS dat and par list.
  # First, SS:
  # Recover number of species
  nS <- omList $ nS
  # Make dat and par lists
  ssDat <- vector ( mode = "list", length = nS )
  ssPar <- vector ( mode = "list", length = nS )
  ssRep <- vector ( mode = "list", length = nS )
  for (s in 1:nS )
  {
    # Make dat list
    ssDat[[s]] <- list (  nT      = omList$nT,
                          Ct      = omList$Ct[s,],
                          It      = omList$It[s,],
                          dumm    = 999 )
    # make par list
    ssPar[[s]] <- list (  lnMSY       = log(omList$msy[s]),
                          lnFmsy      = log(omList$Fmsy[s]),
                          mMSY        = omList$msy[s],
                          sMSY        = ctl$sMSYmult[s]*omList$msy[s],
                          mFmsy       = omList$Fmsy[s],
                          sFmsy       = ctl$sFMSYmult[s]*omList$Fmsy[s],
                          alpha.sigma = ctl$alpha.sigma[s],
                          beta.sigma  = ctl$beta.sigma[s],
                          alpha.tau   = ctl$alpha.tau[s],
                          beta.tau    = ctl$beta.tau[s],
                          mlnq        = log(omList$q[s]),
                          slnq        = log(3),
                          epst        = omList$epst[s,]+omList$zetat[s,] )

  }
  # now make dat and par lists for the MS model
  msDat <- list ( nT    = omList$nT,
                  nS    = omList$nS,
                  Ct    = omList$Ct,
                  It    = omList$It,
                  dumm  = 999 )
  msPar <- list ( lnMSY       = log(omList$msy),
                  lnFmsy      = log(omList$Fmsy),
                  mMSY        = omList$msy,
                  sMSY        = ctl$sMSYmult*omList$msy,
                  mFmsy       = omList$Fmsy,
                  sFmsy       = ctl$sMSYmult*omList$Fmsy,
                  alpha.sigma = ctl$alpha.sigma,
                  beta.sigma  = ctl$beta.sigma,
                  alpha.tau   = ctl$alpha.tau,
                  beta.tau    = ctl$beta.tau,
                  mlnq        = mean(log(omList$q)),
                  slnq        = log(3),
                  epst        = omList$epst+omList$zetat )

  outlist <- list ( ssDat = ssDat, 
                    ssPar = ssPar, 
                    msDat = msDat, 
                    msPar = msPar )
  return(outlist)
}



# callProcADMB()
# Function that calls an ADMB estimator for a given set of
# dat and par list objects.
# inputs:   dat=list object containing .dat file contents in order
#           par=list object containing .pin file contents in order
#           lab=character vector label for use in output file names
#           fitTrials=number of times to modify pin before giving up
#           activeFileRoot=name of ADMB model being called
#           mcTrials=number of MCMC trials to perform
#           mcSave=thinning value for MCMC trials
# ouputs:   rep=list object containing the rep file contents
#           mcPar=data.frame of MCMC parameter distributions
#           mcBio=data.frame of MCMC biomass distributions
#           localMin=Logical indicating if local minimum was found
#           hessPosDef=logical indicating if MCMC was possible
# usage:    fitting a replicate model run to given data
callProcADMB <- function (  dat=ssDat[[1]], par=ssPar[[1]], lab=NULL, 
                            fitTrials = 3, activeFileRoot="ssProd",
                            mcTrials = 1, mcSave = 1, maxfn = 1000 )
{ 
  # Set the active filename root, paste with label and extensions
  datFile   <- paste ( activeFileRoot,lab,".dat", sep = "")
  pinFile   <- paste ( activeFileRoot,lab,".pin", sep = "")

  # No label on admb output files
  repFile   <- paste ( activeFileRoot,".rep", sep = "")
  mcParFile <- paste ( activeFileRoot,"ParMC.dat", sep = "")
  mcBioFile <- paste ( activeFileRoot,"BioMC.dat", sep = "")

  # Backup par and bio MC output for failed MCMC runs
  mcParBackup <- paste (activeFileRoot,lab,"ParMC.bak", sep = "")
  mcBioBackup <- paste (activeFileRoot,lab,"BioMC.bak", sep = "")


  # Set path, exec and paste together command call
  path <- getwd()
  exec <- file.path(path,activeFileRoot)
  procCall <- paste ( exec, " -ainp ", pinFile, " -ind ", datFile, 
                            " -mcmc ", mcTrials, " -maxfn ", maxfn,
                            " -mcsave ", mcSave, sep = "" )
  mcEval <- paste ( exec, " -ainp ", pinFile, " -ind ", datFile, 
                    " -mceval", sep = "" )

  # Pull out lnMSY and sMSY for use in refitting procedure later
  lnMSY <- par$lnMSY 
  sMSY <- par$sMSY
  # Write out the dat file
  cat ( "## ", activeFileRoot, " data file, created ", 
        format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
        sep = "", file = datFile, append = FALSE )
  lapply (  X = seq_along(dat), FUN = writeADMB, x = dat, 
            activeFile=datFile )

  for ( i in 0:fitTrials )
  {
    hessPosDef <- TRUE
    localMin <- TRUE
    # decrease prior sd if needed
    par$lnMSY <- lnMSY * (fitTrials + i)/fitTrials
    par$sMSY <- sMSY * (fitTrials - i + 1)  / (fitTrials)
    # Write to pin file 
    cat ( "## ", activeFileRoot, " initial parameter file, created ", 
          format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
          sep = "", file = pinFile, append = FALSE )
    lapply (  X = seq_along(par), FUN = writeADMB, x = par, 
              activeFile=pinFile )

    # Now run the model and read in rep and MCMC files
    system ( command = procCall, wait =TRUE, ignore.stdout = TRUE )
    system ( command = mcEval, wait =TRUE, ignore.stdout = TRUE )

    rep <- lisread ( repFile )
    mcPar <- try ( read.table ( mcParFile, header = TRUE) )
    mcBio <- try ( read.table (mcBioFile) )

    # Check if that the exit code is correct in the rep file
    if ((rep$iExit != 1 | rep$maxGrad > 1e-4) )
    {
      localMin <- FALSE
      hessPosDef <- FALSE
      if (i < fitTrials)
      {
        cat ( activeFileRoot, 
              " failed to find local min with given .pin file, repeating ",
              i+1, "th time.\n", sep = "")
      }
    }
    if ( class (mcPar) == "try-error" | class (mcBio) == "try-error" )
    {
      hessPosDef <- FALSE 
      if (i < fitTrials)
      {
        cat ( activeFileRoot, 
              " has non-positive definite hessian with given .pin file, repeating ",
              i+1, "th time.\n", sep = "")
      }
    }
    if ( localMin & hessPosDef )
    # if (hessPosDef) 
      break
  }
  if ( !localMin | !hessPosDef )
  {
    if ( !hessPosDef )
    {
      mcPar <- read.table (mcParBackup, header = TRUE )
      mcBio <- read.table (mcBioBackup)
    }
    cat(  "\n", activeFileRoot, "model had trouble fitting in", i, 
          "trials with given pin and dat file.\n",
          "Copying bad pin and dat for forensics.\n" )
    # Save the *.pin and *.dat files for forensics.
    badPinFile <- file.path(  getwd(),"badfits",
                              paste( activeFileRoot,lab,Sys.time(),".pin", sep="" ) )
    file.copy( pinFile, badPinFile, overwrite=TRUE )
    badDatFile <- file.path(  getwd(),"badfits",
                              paste( activeFileRoot,lab,Sys.time(),".dat", sep="" ) )
    file.copy( datFile, badDatFile, overwrite=TRUE )
  } else { 
    # copy Par and Bio MCMC output for later (rep too??)
    file.copy( mcParFile, mcParBackup, overwrite = TRUE )
    file.copy( mcBioFile, mcBioBackup, overwrite = TRUE )
  }

  out <- list ( rep = rep, mcPar = mcPar, mcBio = mcBio,
                localMin = localMin, hessPosDef = hessPosDef )
  out
}

# seedFit()
# seedFit takes a given seed value and control list, and
# returns a full run of the simulation-estimation procedure using
# the supplied seed values and the control list parameters
# inputs:   seed=integer value provided to the rng
#           ctl=control-list set in controlFile.txt
# ouputs:   om=list of operating model quantities
#           ssFit=nS-list of ss model outputs
#           msFit=list of ms model outputs
# usage:    to be lapplied over seed values as part of a monte-carlo trial
seedFit <- function ( seed = 1, ctl, quiet = FALSE )
{
  # Run operating model to generate data
  om <- opModel ( control = ctl, seed = seed )

  # Create data objects for EMs
  datPar <- makeDataLists ( om, ctlList)

  # Call EMs
  ssFit <- list()
  for (s in 1:ctlList$nS ) 
  {
    if (!quiet) cat ( "Fitting species ", s, " of ", ctl$nS, "\n", sep = "")
    ssFit[[s]] <- callProcADMB (  dat = datPar$ssDat[[s]], 
                                  par = datPar$ssPar[[s]],
                                  lab=s, fitTrials = ctl$fitTrials, 
                                  activeFileRoot = "ssProd",
                                  mcTrials = ctl$mcTrials, 
                                  mcSave = 1, maxfn = 5000)
  }
  if (!quiet) cat ( "Fitting ", ctlList$nS," species hierarchically.\n", 
                    sep = "")
  msFit <- callProcADMB ( dat = datPar$msDat, par = datPar$msPar,
                          fitTrials = ctl$fitTrials, 
                          activeFileRoot = "msProd",
                          mcTrials = ctl$mcTrials, 
                          mcSave = 1, maxfn = 5000 )

  cat ("Completed replicate ", seed - ctl$rSeed, " of ", ctl$nReps, ".\n", sep = "")
  # Return fits
  out <- list ( om = om,
                ssFit = ssFit,
                msFit = msFit )
  out
}

# simEstProc()
# This function runs the entire sim-est procedure, taking inputs from
# the control list and returning a list object containing all operating
# model and assessment model output. Somewhat inspired by .mgmtProc()
# in mseR
# inputs:   ctl=control list containing settings
#           parallel=logical switch to turn on parallel processing
# outputs:  blob=list object containing all sim-est details for every rep
# usage:    inside a run function wrapper  
simEstProc <- function ( ctl, quiet = TRUE )
{
  # Recover simulation controls from the control list
  nReps <- ctl$nReps
  nS    <- ctl$nS
  nT    <- ctl$nT
  rSeed <- ctl$rSeed

  # Create list object to store all simulated values
  om <- list (  Bt    = array (NA,dim = c(nReps,nS,nT) ),
                Ct    = array (NA,dim = c(nReps,nS,nT) ),
                Ft    = array (NA,dim = c(nReps,nS,nT) ),
                epst  = array (NA,dim = c(nReps,nS,nT) ),
                zetat = array (NA,dim = c(nReps,nS,nT) ),
                deltat= array (NA,dim = c(nReps,nS,nT) ),
                It    = array (NA,dim = c(nReps,nS,nT) ),
                ItTrue= array (NA,dim = c(nReps,nS,nT) ),
                dep   = matrix(NA,nrow = nReps, ncol = nS)
              )
                # sigma = numeric(length = nReps),
                # Sigma = numeric(length = nReps),
                # tau   = numeric(length = nReps),
                # q     = numeric(length = nReps),
                # nT    = numeric(length = nReps),
                # nS    = numeric(length = nReps),
                # msy   = numeric(length = nReps),
                # Fmsy  = numeric(length = nReps) )

  # Now create a list object to ss and ms assessment model outputs
  am <- list ( ss=NULL, ms=NULL)

  # single species
  am$ss <- list ( msy   = matrix ( NA, nrow = nReps, ncol = nS ),
                  Fmsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                  Bmsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                  q     = matrix ( NA, nrow = nReps, ncol = nS ),
                  sigma2= matrix ( NA, nrow = nReps, ncol = nS ),
                  tau2  = matrix ( NA, nrow = nReps, ncol = nS ),
                  dep   = matrix ( NA, nrow = nReps, ncol = nS ),
                  epst  = array (NA, dim = c(nReps,nS,nT)),
                  Bt    = array (NA, dim = c(nReps,nS,nT)),
                  Ft    = array (NA, dim = c(nReps,nS,nT)),
                  mlnq  = matrix ( NA, nrow = nReps, ncol = nS ),
                  slnq  = matrix ( NA, nrow = nReps, ncol = nS ),
                  mcPar = vector ( mode = "list", length = nReps),
                  mcBio = vector ( mode = "list", length = nReps),
                  locmin= matrix ( NA, nrow = nReps, ncol = nS),
                  hesspd= matrix ( NA, nrow = nReps, ncol = nS) )

  # make a sublist in the mcPar and mcBio lists, for each species
  # Change to a 4d array later, once the # of par estimates is settled
  for ( i in 1:nReps)
  {
    am$ss$mcPar[[i]]     <- vector(mode = "list", length = nS)
    am$ss$mcBio[[i]]     <- vector(mode = "list", length = nS)  
  }
  

  # multispecies (coastwide)
  am$ms <- list ( msy   = matrix ( NA, nrow = nReps, ncol = nS ),
                  Fmsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                  Bmsy  = matrix ( NA, nrow = nReps, ncol = nS ),
                  q     = matrix ( NA, nrow = nReps, ncol = nS ),
                  sigma2= matrix ( NA, nrow = nReps, ncol = nS ),
                  tau2  = matrix ( NA, nrow = nReps, ncol = nS ),
                  dep   = matrix ( NA, nrow = nReps, ncol = nS ),
                  epst  = array (NA, dim = c(nReps,nS,nT)),
                  Bt    = array (NA, dim = c(nReps,nS,nT)),
                  Ft    = array (NA, dim = c(nReps,nS,nT)),
                  mlnq  = vector ( "numeric", length=nReps ),
                  slnq  = vector ( "numeric", length=nReps),
                  mcPar = vector ( mode = "list", length = nReps),
                  mcBio = vector ( mode = "list", length = nReps),
                  locmin= vector ( mode = "logical", length = nReps),
                  hesspd= vector ( mode = "logical", length = nReps) )

  # make a sublist in the mcPar and mcBio lists, for each species
  # Change to a 4d array later, once the # of par estimates is settled
  for (i in 1:nReps)
  {
    am$ss$mcPar[[i]]     <- vector(mode = "list", length = nS)
    am$ss$mcBio[[i]]     <- vector(mode = "list", length = nS)  
  }
  

  # The BLOOOOOOBBBBBB
  blob <- list ( ctl = ctl, om = om, am = am)

  # Now create a vector of seed values, able to be lapplied over
  seeds <- rSeed + 1:nReps

  simEst <- lapply ( X = seeds, FUN = seedFit, ctl = blob$ctl,
                     quiet = quiet )

  cat ( "Finished fitting ", nReps, " replicates.\n",
        "Now saving output.", sep = "" )

  for ( i in 1:length(simEst) )
  {
    # Save OM values
    blob$om$Bt[i,,]         <- simEst[[i]]$om$Bt
    blob$om$Ct[i,,]         <- simEst[[i]]$om$Ct
    blob$om$epst[i,,]       <- simEst[[i]]$om$epst
    blob$om$Ft[i,,]         <- simEst[[i]]$om$Ft
    blob$om$zetat[i,,]      <- simEst[[i]]$om$zetat
    blob$om$deltat[i,,]     <- simEst[[i]]$om$Bt
    blob$om$It[i,,]         <- simEst[[i]]$om$It
    blob$om$ItTrue[i,,]     <- simEst[[i]]$om$ItTrue
    blob$om$dep[i,]         <- simEst[[i]]$om$dep

    # Save AM results
    # Loop over single species
    for ( s in 1:nS )
    {
      blob$am$ss$msy[i,s]         <- simEst[[i]]$ssFit[[s]]$rep$MSY
      blob$am$ss$Fmsy[i,s]        <- simEst[[i]]$ssFit[[s]]$rep$FMSY
      blob$am$ss$Bmsy[i,s]        <- simEst[[i]]$ssFit[[s]]$rep$BMSY
      blob$am$ss$q[i,s]           <- simEst[[i]]$ssFit[[s]]$rep$q
      blob$am$ss$sigma2[i,s]      <- simEst[[i]]$ssFit[[s]]$rep$sigma2
      blob$am$ss$tau2[i,s]        <- simEst[[i]]$ssFit[[s]]$rep$tau2
      blob$am$ss$dep[i,s]         <- simEst[[i]]$ssFit[[s]]$rep$D
      blob$am$ss$epst[i,s,]       <- simEst[[i]]$ssFit[[s]]$rep$epst
      blob$am$ss$Bt[i,s,]         <- simEst[[i]]$ssFit[[s]]$rep$Bt
      blob$am$ss$Ft[i,s,]         <- simEst[[i]]$ssFit[[s]]$rep$Ft
      blob$am$ss$mlnq[i,s]        <- simEst[[i]]$ssFit[[s]]$rep$mlnq
      blob$am$ss$slnq[i,s]        <- simEst[[i]]$ssFit[[s]]$rep$slnq
      blob$am$ss$mcPar[[i]][[s]]  <- simEst[[i]]$ssFit[[s]]$mcPar
      blob$am$ss$mcBio[[i]][[s]]  <- simEst[[i]]$ssFit[[s]]$mcBio
      blob$am$ss$locmin[i,s]      <- simEst[[i]]$ssFit[[s]]$localMin
      blob$am$ss$hesspd[i,s]      <- simEst[[i]]$ssFit[[s]]$hessPosDef
    }
    # Now multispecies
    blob$am$ms$msy[i,]            <- simEst[[i]]$msFit$rep$MSY
    blob$am$ms$Fmsy[i,]           <- simEst[[i]]$msFit$rep$FMSY
    blob$am$ms$Bmsy[i,]           <- simEst[[i]]$msFit$rep$BMSY
    blob$am$ms$q[i,]              <- simEst[[i]]$msFit$rep$q
    blob$am$ms$sigma2[i,]         <- simEst[[i]]$msFit$rep$sigma2
    blob$am$ms$tau2[i,]           <- simEst[[i]]$msFit$rep$tau2
    blob$am$ms$dep[i,]            <- simEst[[i]]$msFit$rep$D
    blob$am$ms$epst[i,,]          <- simEst[[i]]$msFit$rep$epst
    blob$am$ms$Bt[i,,]            <- simEst[[i]]$msFit$rep$Bt
    blob$am$ms$Ft[i,,]            <- simEst[[i]]$msFit$rep$Ft
    blob$am$ms$mlnq[i]            <- simEst[[i]]$msFit$rep$mlnq
    blob$am$ms$slnq[i]            <- simEst[[i]]$msFit$rep$slnq
    # Split up mcPar and mcBio for the ms model
    for ( s in 1:nS )
    {
      parIdx <- seq ( from = s, to = nrow(simEst[[i]]$msFit$mcPar), by = nS)
      bioIdx <- seq ( from = s, to = nrow(simEst[[i]]$msFit$mcBio), by = nS)
      blob$am$ms$mcPar[[i]][[s]]  <- simEst[[i]]$msFit$mcPar[parIdx,]
      blob$am$ms$mcBio[[i]][[s]]  <- simEst[[i]]$msFit$mcBio[bioIdx,]
    }
    blob$am$ms$locmin[i]          <- simEst[[i]]$msFit$localMin
    blob$am$ms$hesspd[i]          <- simEst[[i]]$msFit$hessPosDef
  }

  blob
}



# logProdModel()
# Function to simulate a logistic based surplus production model,
# assuming that population is initially at equilibrium (2*Bmsy).
# inputs:     msy=maximum sustainable yield; Fmsy=fishing mortality for MSY
#             nT=length of simulation; Ct=nT-vector of catch (opt)
# 						Ft=nT-vector of fishing mortality (opt)
#             epst=nT-vector of proc errors (could be RW, or IID)
#             sigma = process error sd (uncorr)
# outputs:    Bt=nT-vector of modeled biomass; Ct=nT-vector of catch
# 						Ut=vector of exploitation rates; 
# usage:      Bexp(eps) = indices of abundance for estimation w/ logN errors eps
logProdModel <- function ( msy = 1, Fmsy = 0.1, nT = 50, Ft = NULL, Ct = NULL,
                           epst = rnorm ( n = nT ), sigma = 0.2, 
                           zetat = rep(0,nT), Sigma = 0 )
{
  # First, initialise a vector to hold biomass
  Bt <- numeric ( length = nT )

  # Multiply errors by variance
  epst <- epst * sigma
  zetat <- zetat * Sigma

  # Compute bias correction
  sigma2 <- sigma * sigma
  Sigma2 <- Sigma * Sigma

  if ( is.null ( Ft ) & is.null ( Ct ) ) 
  {
    cat ( "You're missing harvest, stupid.\n" )
    return ()
  }

  if ( is.null ( Ct ) ) 
  {
  	Ct <- numeric ( length = nT )
  }

  # Compute BMSY
  Bmsy <- msy / Fmsy

  # Populate the vector, with a special case for t = 1
  Bt [ 1 ] <- 2 * Bmsy * exp ( (epst[1] - sigma2 / 2.) + (zetat[1] - Sigma2/2.) )

  # Loop over remaining years
  for ( t in 2:nT )
  {
    if ( Ct [ t-1 ] == 0 ) Ct [ t-1 ] <- Ft[t-1] * Bt [ t-1 ]
    Bt [ t ] <-   (	Bt [ t-1 ] +                       # prev year's B
                  	2. * Fmsy * Bt [ t-1 ] * 
                  	( 1. - Bt [ t-1 ] / Bmsy / 2. ) -  # recruitment
                  	Ct [ t-1 ]                         # Catch
                  )
    ## THIS IS BOGUS ##
    Bt[t] <- max(Bt[t], 1e-1)

    Bt[t] <- Bt[t] * exp ( epst [ t ] - sigma2 / 2. +
                            zetat[ t ] - Sigma2 / 2. ) # Process error
}  # End loop for biomass
  # Generate catch for nT
  Ct[nT] <- Ft[nT] * Bt[nT]

  # compute depletion value for comparison to model fit
  dep <- Bt [ nT ] / Bmsy / 2.

  # Return B, C, F and U series and
  # depletion value
  outList <- list ()
  outList $ Bt <- Bt
  outList $ Ct <- Ct
  outList $ Ft <- Ft
  # outList $ Ut <- Ut
  outList $ epst <- epst
  outList $ zetat <- zetat
  outList $ dep <- dep
  return ( outList )
}

# Why keep the observation model separate???? 
# Right now just for modularity.

# obsModel()
# A function which creates observations (indices of biomass) from
# given true series of biomass and observation model parameters.
# inputs:			Bt=true series of biomass; q=survey catchability parameter
#							nT=length of time series; deltat=nT-vector of standard errors
# 						tau=obs error standard deviation (uncorrelated)
# outputs:		It=nT-vector of survey observations
# usage:			creates survey observations which are supplied to the estimator
obsModel <- function ( 	Bt, q = 0.05, nT = length (Bt), 
												deltat = rnorm(nT), tau = 0.4 )
{
  deltat <- deltat*tau
  tau2 <- tau*tau
	It <- q * Bt * exp ( deltat - tau2/2.0 )
	outList <- list ( It = It, deltat = deltat )
}


#### DEPRECATED, now callProcADMB works for both SS and MS models ####
# callMSprodADMB() 
# Function that calls the ADMB ssProd estimator for a given stock.
# inputs:   dat=list object containing .dat file contents in order
#           par=list object containing .pin file contents in order
# ouputs:   rep=list object containing the rep file contents
# usage:    fitting a replicate model run
callMSProdADMB <- function (  dat = msDat, par = msPar,
                              fitTrials = 3 )
{ 
  # Set the active filename root
  activeFileRoot <- "msProd"
  datFile <- paste ( activeFileRoot, ".dat", sep = "")
  pinFile <- paste ( activeFileRoot, ".pin", sep = "")

  # Write out the dat file
  cat ( "## Single Species Production Model data file, created ", 
        format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
        sep = "", file = datFile, append = FALSE )
  lapply (  X = seq_along(dat), FUN = writeADMB, x = dat, 
            activeFile=datFile )

  # Write to pin file
  pinFile <- paste ( activeFileRoot, ".pin", sep = "" )
  cat ( "## Single Species Production Model pin file, created ", 
        format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
        sep = "", file = pinFile, append = FALSE )
  lapply (  X = seq_along(par), FUN = writeADMB, x = par, 
            activeFile=pinFile )

  # Now run the model
  path <- getwd()
  exec <- file.path(path,"msProd")
  mleCall <- paste ( exec, " -ainp ", pinFile, " -ind ", datFile, 
                        " -maxfn 5000", sep = "" )
  system ( command = mleCall, wait =TRUE, ignore.stdout = TRUE )

  repFile <- paste ( activeFileRoot, ".rep", sep = "" )
  msRep <- lisread ( repFile )
  msRep
}

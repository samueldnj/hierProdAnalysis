# --------------------------------------------------------------------------
# simulation.R
# 
# Script for the simulation procedure to produce data for ADMB and TMB
# estimators.
# 
# Author: Samuel Johnson
# Date: 24 August, 2016
# 
# --------------------------------------------------------------------------

# runSimEst()
# Wrapper function to run the entire sim-est procedure for a given set of
# control parameters, calculate errors and save results to a (given) folder
# inputs:   ctlFile=character with the name/path of the control file
#           folder=optional character name of output folder 
#                   ie saves to ./project/folder
# ouputs:   NULL
# usage:    from the console to run the procedure
runSimEst <- function ( ctlFile = "simCtlFile.txt", folder=NULL )
{ 
  # read in control file
  controlList <- .readParFile ( ctlFile )
  controlList <- .createList  ( controlList )
  # Run simEst Procedure
  blob <- .simEstProc( obj = controlList, quiet = TRUE )
  # Make error distributions
  blob <- .makeRelErrorDists ( blob )
  
  # save output to project folder
  # First, if a folder name isn't nominated, create a default sim folder
  if ( is.null(folder) )
  {
    stamp <- paste( format(Sys.time(),format="%d%m%Y%H%M%S" ),sep="" )
    folder <- paste ( "sim_",stamp, sep = "" )
  } else folder <- paste ("sim_", folder, sep = "")
  # Now paste together the path to the folder and create it
  path <- file.path (getwd(),"project",folder)
  dir.create ( path )
  cat( "\nMSG (saveSim) Created simulation folder ",folder,"in project.\n" )
  # Save blob
  .saveSim(blob=blob, name=folder, path = path)
  # Copy control file and ms par file to folder for posterity
  file.copy(from=ctlFile,to=file.path(path,ctlFile))
  file.copy(from=blob$opMod$parFile,to=file.path(path,blob$opMod$parFile))
  # Done
}


# opModel()
# The operating model takes control list parameters and simulates biological 
# and observational data for each species, creating an om list object that 
# can then be used to create dat and pin files to be passed to estimators.
# inputs:     control=control list object, read in from control file
#             seed=random number seed
# ouputs:     om=list object containing the data for each species
# usage:      to prepare data for each replicate (seed value)
.opModel <- function (obj, seed = 1)
{
  # First, set the seed
  if ( !is.null(seed) ) set.seed(seed)

  # Now recover multi-species parameters
  nT <- obj$opMod$nT          # length of simulation
  nS <- obj$opMod$nS          # number of species

  # Process error component vars
  sigma <- sqrt(obj$opMod$pars$sigma)    # longitudinal shared proc error sd
  Sigma <- sqrt(obj$opMod$pars$Sigma)    # ms proc error cov mtx
  
  # Rescale shared effects sd if desired
  if (!is.null(obj$opMod$sigmaMult)) sigma <- sigma*obj$opMod$sigmaMult

  # Correlation parameters
  rho     <- obj$opMod$rho          # longitudinal auto-correlation
  msCorr  <- obj$opMod$pars$msCorr  # ms cross correlation

  # Multiply the off-diagonal elements of msCorr if desired
  if (!is.null(obj$opMod$corrMult))
  {
    corrMult <- obj$opMod$corrMult
    msCorr <- corrMult*msCorr
    diag(msCorr) <- 1
  }

  # Obs error var
  tau <- obj$opMod$tau

  # Initialise list to hold the data
  om <- list (  Bt    = matrix (NA,nrow=nS, ncol=nT ),
                Ct    = matrix (NA,nrow=nS, ncol=nT ),
                Ft    = matrix (NA,nrow=nS, ncol=nT ),
                ItTrue= matrix (NA,nrow=nS, ncol=nT ),
                epst  = numeric (length=nT ),
                zetat = matrix (NA,nrow=nS, ncol=nT ),
                deltat= matrix (NA,nrow=nS, ncol=nT ),
                It    = matrix (NA,nrow=nS, ncol=nT ),
                dep   = numeric(length = nS),
                sigma2= sigma*sigma,
                Sigma2= Sigma*Sigma,
                msCorr= msCorr,
                rho   = rho,
                tau   = tau,
                q     = obj$opMod$q,
                nT    = nT,
                nS    = nS,
                Bmsy  = obj$opMod$pars$Bmsy,
                Fmsy  = obj$opMod$pars$Fmsy )

  # Now create epst and zetat vectors using the proc error components
  epst      <- .fillRanWalk(  z=rnorm(n = nT), s=sigma,
                              rho=rho )
  zetat     <- matrix  (rnorm ( n = nS*nT ), nrow = nS, ncol = nT)
  zetat     <- .genCorrDevs ( zetat,
                              Mcorr=msCorr,
                              Sigma=Sigma )

  # standard normal deviations for obs error (uncorrelated)
  deltat    <- matrix(rnorm ( n = nS*nT ), nrow = nS, ncol = nT) 

  # Create Ft time series
  ## REPLACE WITH A FUNCTION LATER FOR MULTIPLE correlated Ft TRAJECTORIES ##
  FtProp <- numeric ( length = nT )
  Fmult <- obj$opMod$Fmult
  # Then get time of peak fishing mortality
  tFpeak <- obj$opMod$tFpeak
  # Compute gradients of the pieces
  Fgrad1 <- ( Fmult [ 2 ] - Fmult [ 1 ] ) / (tFpeak - 1 )
  Fgrad2 <- ( Fmult [ 3 ] - Fmult [ 2 ] ) / 
                  ( nT - tFpeak )
  # Populate Ft with multipliers of F_MSY
  FtProp [ 1:tFpeak ]         <- Fmult [ 1 ] + 0:(tFpeak - 1) * Fgrad1
  FtProp [ (tFpeak + 1):nT ]  <- Fmult [ 2 ] + 1:(nT - tFpeak) * Fgrad2

  # Overwrite Ft with actual F values
  Fst <- matrix ( FtProp, nrow = nS, ncol = nT, byrow = TRUE )
  for ( s in 1:nS ) Fst[s,] <- obj$opMod$pars$Fmsy[s] * Fst[s,]

  # Loop over species and fill list entries with biological
  # and observational data
  for ( s in 1:nS )
  {
    bio <- .logProdModel ( Bmsy = obj$opMod$pars$Bmsy[s], 
                          Fmsy = obj$opMod$pars$Fmsy[s],
                          nT = nT, Ft = Fst[s,], 
                          epst = epst,
                          zetat = zetat[s,] )
    obs <- .obsModel ( Bt = bio $ Bt, q = obj$opMod$q[s], nT = nT,
                      deltat = deltat[s,], tau = tau[s] )

    om$Bt[s,]     <- bio$Bt
    om$Ct[s,]     <- bio$Ct
    om$Ft[s,]     <- bio$Ft
    om$ItTrue[s,] <- obj$opMod$q[s]*bio$Bt
    om$It[s,]     <- obs$It
    om$zetat[s,]  <- bio$zetat
    om$deltat[s,] <- obs$deltat
    om$dep[s]     <- bio$dep
  }
  om$epst         <- epst

  return(om)
}

# makeDataLists()
# The makeDataLists function takes an om list object produced by opModel()
# and creates data structures that will be passed to estimators
# inputs:   obj=list object containing control and om lists
# outputs:  msDat=multispecies data list; 
#           msPar=multispecies initial parameter value list
#           ssDat=nS-list of single species data lists
#           ssPar=nS-list of single species init. par lists
.makeDataLists <- function ( obj )
{
  # Recover om and control lists
  om <- obj$om

  # Needs to create nS SS dat and par lists, and 1 MS dat and par list.
  # First, SS:
  # Recover number of species
  nS <- obj$opMod$ nS
  nT <- obj$opMod$nT
  # Make dat and par lists
  ssDat <- vector ( mode = "list", length = nS )
  ssPar <- vector ( mode = "list", length = nS )

  # Sum catch to get starting biomass
  sumCat <- apply ( X = obj$om$Ct, MARGIN=1, FUN = sum )
  sumCat <- as.numeric(sumCat)
  maxF <- apply ( X = obj$om$Ft, MARGIN=1,FUN=max)
  maxF <- as.numeric(maxF)

  # loop over species
  for (s in 1:nS )
  {
    # Make dat list
    ssDat[[s]] <- list (  nT          = om$nT,
                          Ct          = om$Ct[s,],
                          It          = om$It[s,],
                          phz.Bmsy    = obj$assess$phz_Bmsy,
                          phz_Fmsy    = obj$assess$phz_Fmsy,
                          phz_tau     = obj$assess$phz_tau,
                          phz_sigma   = obj$assess$phz_sigma,
                          phz_q       = obj$assess$phz_q,
                          phz_mlnq    = obj$assess$phz_mlnq,
                          phz_slnq    = obj$assess$phz_slnq,
                          phz_dev     = obj$assess$phz_dev,
                          phz_AR      = obj$assess$phz_AR,
                          dumm        =  999 )
    # make par list
    ssPar[[s]] <- list (  lnBmsy      = log(sumCat[s]),
                          lnFmsy      = log(om$Fmsy[s]),
                          lnTau2      = log(om$tau[s]),
                          lnsigma2    = log(om$sigma+om$Sigma[s]),
                          # lnq         = log(obj$om$q[s]),
                          mBmsy       = obj$assess$mBmsy[s],
                          sBmsy       = obj$assess$sBmsy[s],
                          mFmsy       = obj$assess$mFmsy[s],
                          sFmsy       = obj$assess$sFmsy[s],
                          alpha_sigma = obj$assess$alpha_sigma, 
                          beta_sigma  = obj$assess$beta_sigma,
                          alpha_tau   = obj$assess$alpha_tau[s],
                          beta_tau    = obj$assess$beta_tau[s],
                          mlnq        = obj$assess$mlnq,
                          slnq        = obj$assess$slnq,
                          rho         = obj$opMod$rho,
                          epst        = rep(0,nT) )

  }
  # now make dat and par lists for the MS model
  msDat <- list ( nT          = om$nT,
                  nS          = om$nS,
                  Ct          = om$Ct,
                  It          = om$It,
                  phz_Bmsy    = obj$assess$phz_Bmsy,
                  phz_Fmsy    = obj$assess$phz_Fmsy,
                  phz_tau     = obj$assess$phz_tau,
                  phz_sigma   = obj$assess$phz_sigma,
                  phz_q       = obj$assess$phz_q,
                  phz_mlnq    = obj$assess$phz_mlnq,
                  phz_slnq    = obj$assess$phz_slnq,
                  phz_dev     = obj$assess$phz_dev,
                  phz_AR      = obj$assess$phz_AR,
                  phz_chol    = obj$assess$phz_chol,
                  dumm        = 999 )
  msPar <- list ( lnBmsy      = log(sumCat),
                  lnFmsy      = log(om$Fmsy),
                  lnTau2      = log(om$tau),
                  lnSigma2    = log(om$Sigma),
                  lnsigma2    = log(om$sigma),
                  # lnq         = log(om$q),
                  mBmsy       = obj$assess$mBmsy,
                  sBmsy       = obj$assess$sBmsy,
                  mFmsy       = obj$assess$mFmsy,
                  sFmsy       = obj$assess$sFmsy,
                  alpha_sigma = obj$assess$alpha_sigma,
                  beta_sigma  = obj$assess$beta_sigma,
                  alpha_Sigma = obj$assess$alpha_Sigma,
                  beta_Sigma  = obj$assess$beta_Sigma,
                  alpha_tau   = obj$assess$alpha_tau,
                  beta_tau    = obj$assess$beta_tau,
                  mlnq        = obj$assess$mlnq,
                  slnq        = obj$assess$slnq,
                  rho         = obj$opMod$rho,
                  c           = rep(0,(nS*(nS-1)/2)),
                  epst        = rep(0,nT),
                  zetat       = matrix(0,nS,nT) )
  # return list of dat and pin objects for running estimators
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
.callProcADMB <- function ( dat=ssDat[[1]], par=ssPar[[1]], lab=NULL, 
                            fitTrials = 3, activeFileRoot="ssProd",
                            mcTrials = 1, mcSave = 1, maxfn = 10000 )
{ 
  # Set the active filename root, paste with label and extensions
  datFile   <- paste ( activeFileRoot,lab,".dat", sep = "")
  pinFile   <- paste ( activeFileRoot,lab,".pin", sep = "")

  # No label on admb output files
  parFile   <- paste ( activeFileRoot,".par", sep = "")
  repFile   <- paste ( activeFileRoot,".rep", sep = "")
  psvFile   <- paste ( activeFileRoot,".psv", sep = "")
  mcParFile <- paste ( activeFileRoot,"ParMC.dat", sep = "")
  mcBioFile <- paste ( activeFileRoot,"BioMC.dat", sep = "")

  # Set path, exec and paste together command call
  path        <- getwd()
  exec        <- file.path(path,activeFileRoot)
  procCall    <- paste ( exec, " -ainp ", pinFile, " -ind ", datFile, 
                            " -mcmc ", mcTrials, " -maxfn ", maxfn,
                            " -mcsave ", mcSave, " -nox", sep = "" )
  procCallpar <- paste ( exec, " -ainp ", parFile, " -ind ", datFile, 
                            " -mcmc ", mcTrials, " -maxfn ", maxfn,
                            " -mcsave ", mcSave, " -nox", sep = "" )
  mcEval      <- paste (  exec, " -ainp ", pinFile, " -ind ", datFile, 
                          " -mceval", sep = "" )
  mcEvalpar   <- paste (  exec, " -ainp ", pinFile, " -ind ", datFile, 
                          " -mceval", sep = "" )

  # Pull out lnMSY and sMSY for use in refitting procedure later
  lnBmsy   <- par$lnBmsy
  sBMSY    <- par$sBMSY
  # Write out the dat file
  cat ( "## ", activeFileRoot, " data file, created ", 
        format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
        sep = "", file = datFile, append = FALSE )
  lapply (  X = seq_along(dat), FUN = writeADMB, x = dat, 
            activeFile=datFile )

  # Write to pin file 
  cat ( "## ", activeFileRoot, " initial parameter file, created ", 
        format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
        sep = "", file = pinFile, append = FALSE )
  lapply (  X = seq_along(par), FUN = writeADMB, x = par, 
            activeFile=pinFile )

  # Now run the model and read in rep and MCMC files
  prcOut  <- system ( command = procCall, wait =TRUE, ignore.stdout = FALSE, 
                      intern = TRUE, ignore.stderr=TRUE )
  evalOut <- system ( command = mcEval, wait =TRUE, 
                      ignore.stdout = FALSE, intern=TRUE, ignore.stderr=TRUE )
  
  for ( i in 1:fitTrials )
  {
    hessPosDef <- TRUE
    localMin <- TRUE
      
    # Read in output from em
    mcPar   <- try ( read.table ( mcParFile, header = TRUE), silent=TRUE )
    mcBio   <- try ( read.table (mcBioFile), silent=TRUE )
    fitrep  <- lisread ( repFile )

    # Check if there was MCMC output, if not, go down the list of possible 
    # reasons, and rerun the fit
    if (!is.null(attr(prcOut,"status"))|!is.null(attr(evalOut,"status")))
    {
      hessPosDef <- FALSE
      # Did we run out of iterations (max grad ceiling is v. high here,
      # but this sometimes works)
      if (fitrep$iExit == 3 & abs(fitrep$maxGrad) < 1e3 )
      {
        localMin <- FALSE
        cat ( activeFileRoot, 
              " reached max iterations,\n repeat i = ",
              i, ".\n", sep = "")
        # Re-run from where the last one stopped
        prcOut <-   system (  command = procCallpar, wait =TRUE, 
                              ignore.stdout = TRUE, intern=TRUE,
                              ignore.stderr=TRUE )
        evalOut <-  system ( command = mcEvalpar, wait =TRUE, 
                             ignore.stdout = TRUE, intern=TRUE,
                             ignore.stderr=TRUE )
        next
      }
      # If the model stopped for another reason, but the max gradient is
      # too high, or the positive penalty was invoked, then Increment
      # lnBmsy and rerun.
      if ( fitrep$iExit <= 1 & (abs(fitrep$maxGrad) > 1e-4 | fitrep$fpen > 0 ) )
      {
        localMin <- FALSE
        cat ( activeFileRoot, 
              " didn't find local min,\n repeat i = ",
              i, ".\n", sep = "")
        # Increment lnBMSY and tighten sBMSY
        par$lnBmsy <- lnBmsy *(1+(i)*0.2)
        # Write to pin file 
        cat ( "## ", activeFileRoot, " initial parameter file, created ", 
              format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
              sep = "", file = pinFile, append = FALSE )
        lapply (  X = seq_along(par), FUN = writeADMB, x = par, 
                  activeFile=pinFile )
        # Now run the model and read in rep and MCMC files
        prcOut <-   system (  command = procCall, wait =TRUE, 
                              ignore.stdout = TRUE, intern=TRUE,
                              ignore.stderr=TRUE )
        evalOut <-  system ( command = mcEval, wait =TRUE, 
                             ignore.stdout = TRUE, intern=TRUE,
                             ignore.stderr=TRUE )
        next
      }      
    }
    # Break if hess Positive definite
    if (hessPosDef) break
    # Otherwise, rerun the model using the par as the pin
    prcOut <-   system (  command = procCallpar, wait =TRUE, 
                          ignore.stdout = TRUE, intern=TRUE,
                          ignore.stderr=TRUE )
    evalOut <-  system ( command = mcEvalpar, wait =TRUE, 
                         ignore.stdout = TRUE, intern=TRUE,
                         ignore.stderr=TRUE )
  }
  # Now check the MLE fit, but don't bother to refit if we've got 
  # a posterior (no matter the quality)
  if (fitrep$maxGrad > 1e-4) localMin <- FALSE

  # # Now if the hessian is NPD, get the posterior for another fit
  # if ( !hessPosDef )
  # {
  #   mcPar <- read.table (mcParBackup, header = TRUE )
  #   mcBio <- read.table (mcBioBackup)
  #   cat(  "\n", activeFileRoot, " had trouble finding PD hessian in ", 
  #         i, " trials with given pin.\n",
  #         "Copying bad pin and dat for forensics.\n", sep="" )
  #   # Save the *.pin and *.dat files for forensics.
  #   badPinFile <- file.path(  getwd(),"badfits",
  #                             paste( activeFileRoot,lab,Sys.time(),".pin", sep="" ) )
  #   file.copy( pinFile, badPinFile, overwrite=TRUE )
  #   badDatFile <- file.path(  getwd(),"badfits",
  #                             paste( activeFileRoot,lab,Sys.time(),".dat", sep="" ) )
  #   file.copy( datFile, badDatFile, overwrite=TRUE )
  # } else { 
  #   # copy Par and Bio MCMC output for use in case of bad fits
  #   file.copy( mcParFile, mcParBackup, overwrite = TRUE )
  #   file.copy( mcBioFile, mcBioBackup, overwrite = TRUE )
  # }

  # Delete output that will trick procedure in the future
  file.remove(mcParFile,mcBioFile)
  if (file.exists(psvFile)) file.remove(psvFile)

  # Return
  out <- list ( fitrep = fitrep, 
                mcPar = mcPar, 
                mcBio = mcBio,
                localMin = localMin, 
                hessPosDef = hessPosDef )
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
.seedFit <- function ( seed = 1, obj, quiet = FALSE )
{

  # Run operating model to generate data
  obj$om <- .opModel ( obj, seed = seed )

  # Create data objects for EMs
  datPar <- .makeDataLists ( obj )

  if(is.null(obj$assess$mcBurn)) obj$assess$mcBurn <- 0

  # Call EMs
  ssFit <- list()
  for (s in 1:obj$opMod$nS ) 
  {
    if (!quiet) cat ( "Fitting species ", s, " of ", obj$ctrl$nS, "\n", sep = "")
    ssFit[[s]] <- .callProcADMB ( dat = datPar$ssDat[[s]], 
                                  par = datPar$ssPar[[s]],
                                  lab=s, fitTrials = obj$assess$fitTrials, 
                                  activeFileRoot = "ssProdCV",
                                  mcTrials = obj$assess$mcTrials+obj$assess$mcBurn, 
                                  mcSave = obj$assess$mcAve, 
                                  maxfn = obj$assess$maxfn)
  }
  if (!quiet) cat ( "Fitting ", obj$opMod$nS," species hierarchically.\n", 
                    sep = "")
  msFit <- .callProcADMB (  dat = datPar$msDat, par = datPar$msPar,
                            fitTrials = obj$assess$fitTrials, 
                            activeFileRoot = "msProdCV",
                            mcTrials = obj$assess$mcTrials+obj$assess$mcBurn, 
                            mcSave = obj$assess$mcSave, 
                            maxfn = obj$assess$maxfn )

  cat ( "Completed replicate ", seed - obj$ctrl$rSeed, " of ", 
        obj$ctrl$nReps, ".\n", sep = "" )

  # Return fits
  out <- list ( om = obj$om,
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
.simEstProc <- function ( obj, quiet = TRUE )
{
  # Recover simulation controls from the control list
  nReps   <- obj$ctrl$nReps
  rSeed   <- obj$ctrl$rSeed
  nS      <- obj$opMod$nS
  nT      <- obj$opMod$nT

  # MCMC output control parameters
  mcTrials<- obj$assess$mcTrials
  mcSave  <- obj$assess$mcSave
  mcBurn  <- obj$assess$mcBurn
  nMCparSS<- obj$assess$nParSS
  nMCparMS<- obj$assess$nParMS

  # Compute total number of trials
  nMC <- (mcTrials+mcBurn)/mcSave

  # Create dimension names for arrays
  rNames <- paste("Rep",1:nReps, sep ="")
  sNames <- obj$ctrl$specNames
  mcNo <- paste("mc",1:nMC,sep="")
  mcBt <- paste("B",1:nT,sep="")

  # Create list object to store all simulated values
  om <- list  ( Bt    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                Ct    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                Ft    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                epst  = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                zetat = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                deltat= array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                It    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                ItTrue= array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                dep   = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames))
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
  am$ss <- list ( Fmsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  Bmsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  msy     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  q       = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  sigma2  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  tau2    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  dep     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  epst    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  rho     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  Bt      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  Ft      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  mlnq    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  slnq    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  mcPar   = array (NA,dim=c(nReps,nS,nMC,nMCparSS),dimnames=list(rNames,sNames,mcNo,NULL)),
                  mcBio   = array (NA,dim=c(nReps,nS,nMC,nT),dimnames=list(rNames,sNames,mcNo,mcBt)),
                  locmin  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  hesspd  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  maxGrad = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)))

  # multispecies (coastwide)
  am$ms <- list ( Fmsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  Bmsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  msy     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  q       = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  sigma2  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  tau2    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  dep     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  epst    = matrix(NA,nrow=nReps,ncol=nT,dimnames=list(rNames,1:nT)),
                  zetat   = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  rho     = vector("numeric",length=nReps),
                  chol    = array (NA,dim=c(nReps,nS,nS),dimnames=list(rNames,sNames,sNames)),
                  Bt      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  Ft      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  mlnq    = vector("numeric",length=nReps),
                  slnq    = vector("numeric",length=nReps),
                  mcPar   = array (NA,dim=c(nReps,nS,nMC,nMCparMS),dimnames=list(rNames,sNames,mcNo,NULL)),
                  mcBio   = array (NA,dim=c(nReps,nS,nMC,nT),dimnames=list(rNames,sNames,mcNo,mcBt)),
                  locmin  = vector(mode="logical",length=nReps),
                  hesspd  = vector(mode="logical",length=nReps),
                  maxGrad = vector(mode="logical",length=nReps))


  # The BLOOOOOOBBBBBB
  blob <- obj
  blob$om <- om
  blob$am <- am

  # Now create a vector of seed values, able to be lapplied over
  seeds <- rSeed + 1:nReps
  
  cat ( "Fitting ", nReps, " replicates.\n", sep = "" )

  for ( i in 1:length(seeds) )
  {
    # Run sim-est procedure
    simEst <- .seedFit ( seed = seeds[i], obj=obj, quiet = quiet)

    # Save OM values
    blob$om$Bt[i,,]         <- simEst$om$Bt
    blob$om$Ct[i,,]         <- simEst$om$Ct
    blob$om$epst[i,,]       <- simEst$om$epst
    blob$om$Ft[i,,]         <- simEst$om$Ft
    blob$om$zetat[i,,]      <- simEst$om$zetat
    blob$om$deltat[i,,]     <- simEst$om$Bt
    blob$om$It[i,,]         <- simEst$om$It
    blob$om$ItTrue[i,,]     <- simEst$om$ItTrue
    blob$om$dep[i,]         <- simEst$om$dep

    # Save AM results
    # Loop over single species
    for ( s in 1:nS )
    {
      blob$am$ss$Bmsy[i,s]        <- simEst$ssFit[[s]]$fitrep$Bmsy
      blob$am$ss$Fmsy[i,s]        <- simEst$ssFit[[s]]$fitrep$Fmsy
      blob$am$ss$msy[i,s]         <- simEst$ssFit[[s]]$fitrep$msy
      blob$am$ss$q[i,s]           <- simEst$ssFit[[s]]$fitrep$q
      blob$am$ss$sigma2[i,s]      <- simEst$ssFit[[s]]$fitrep$sigma2
      blob$am$ss$tau2[i,s]        <- simEst$ssFit[[s]]$fitrep$tau2
      blob$am$ss$dep[i,s]         <- simEst$ssFit[[s]]$fitrep$D
      blob$am$ss$epst[i,s,]       <- simEst$ssFit[[s]]$fitrep$epst
      blob$am$ss$rho[i,s]         <- simEst$ssFit[[s]]$fitrep$rho
      blob$am$ss$Bt[i,s,]         <- simEst$ssFit[[s]]$fitrep$Bt
      blob$am$ss$Ft[i,s,]         <- simEst$ssFit[[s]]$fitrep$Ut
      blob$am$ss$mlnq[i,s]        <- simEst$ssFit[[s]]$fitrep$mlnq
      blob$am$ss$slnq[i,s]        <- simEst$ssFit[[s]]$fitrep$slnq
      # Do MCMC output only if hesspd
      if (simEst$ssFit[[s]]$hessPosDef)
      {
        dimnames(blob$am$ss$mcPar)[[4]] <- colnames(simEst$ssFit[[s]]$mcPar)
        blob$am$ss$mcPar[i,s,,]     <- as.matrix(simEst$ssFit[[s]]$mcPar)
        blob$am$ss$mcBio[i,s,,]     <- as.matrix(simEst$ssFit[[s]]$mcBio)
      }
      # Estimator performance flags
      blob$am$ss$locmin[i,s]      <- simEst$ssFit[[s]]$localMin
      blob$am$ss$hesspd[i,s]      <- simEst$ssFit[[s]]$hessPosDef
      blob$am$ss$maxGrad[i,s]     <- simEst$ssFit[[s]]$fitrep$maxGrad
    }

    # Now multispecies
    blob$am$ms$Bmsy[i,]           <- simEst$msFit$fitrep$Bmsy
    blob$am$ms$Fmsy[i,]           <- simEst$msFit$fitrep$Fmsy
    blob$am$ms$msy[i,]            <- simEst$msFit$fitrep$msy
    blob$am$ms$q[i,]              <- simEst$msFit$fitrep$q
    blob$am$ms$sigma2[i,]         <- simEst$msFit$fitrep$sigma2
    blob$am$ms$tau2[i,]           <- simEst$msFit$fitrep$tau2
    blob$am$ms$dep[i,]            <- simEst$msFit$fitrep$D
    blob$am$ms$epst[i,]           <- simEst$msFit$fitrep$epst
    blob$am$ms$zetat[i,,]         <- simEst$msFit$fitrep$zetat
    blob$am$ms$rho[i]             <- simEst$msFit$fitrep$rho
    blob$am$ms$chol[i,,]          <- simEst$msFit$fitrep$chol
    blob$am$ms$Bt[i,,]            <- simEst$msFit$fitrep$Bt
    blob$am$ms$Ft[i,,]            <- simEst$msFit$fitrep$Ut
    blob$am$ms$mlnq[i]            <- simEst$msFit$fitrep$mlnq
    blob$am$ms$slnq[i]            <- simEst$msFit$fitrep$slnq
    # Split up mcPar and mcBio for the ms model
    # might have NPD hessian, so leave as NAs if so
    if (simEst$msFit$hessPosDef)
    {
      for ( s in 1:nS )
      {
        dimnames(blob$am$ms$mcPar)[[4]] <- colnames(simEst$msFit$mcPar)
        parIdx <- seq ( from = s, to = nrow(simEst$msFit$mcPar), by = nS)
        bioIdx <- seq ( from = s, to = nrow(simEst$msFit$mcBio), by = nS)
        blob$am$ms$mcPar[i,s,,]     <- as.matrix(simEst$msFit$mcPar[parIdx,])
        blob$am$ms$mcBio[i,s,,]     <- as.matrix(simEst$msFit$mcBio[bioIdx,])

      }
    }
    # Estimator Performance Flags
    blob$am$ms$locmin[i]          <- simEst$msFit$localMin
    blob$am$ms$hesspd[i]          <- simEst$msFit$hessPosDef
    blob$am$ms$maxGrad[i]         <- simEst$msFit$fitrep$maxGrad
  }

  cat ( "Completed ", nReps, " replicates.\n", sep = "" )

  blob
}


# logProdModel()
# Function to simulate a logistic based surplus production model,
# assuming that population is initially at equilibrium (2*Bmsy) modified
# by a proc error.
# inputs:     msy=maximum sustainable yield; Fmsy=fishing mortality for msy
#             nT=length of simulation; Ct=nT-vector of catch (opt)
# 						Ft=nT-vector of fishing mortality (opt)
#             epst=nT-vector of env (AR(1)) logN proc errors (bias corrected) 
#             zetat=nT-vector of logN proc errors (correlated among species)
# outputs:    Bt=nT-vector of modeled biomass; Ct=nT-vector of catch
# 						Ut=vector of exploitation rates; 
# usage:      Bexp(eps) = indices of abundance for estimation w/ logN errors eps
.logProdModel <- function ( Bmsy = 1, Fmsy = 0.1, nT = 50, Ft = NULL, Ct = NULL,
                            epst = exp(rnorm ( n = nT )-0.5), 
                            zetat = exp(rep(0,nT)) )
{
  # First, initialise a vector to hold biomass
  Bt <- numeric ( length = nT )

  if ( is.null ( Ft ) & is.null ( Ct ) ) 
  {
    cat ( "You're missing harvest, stupid.\n" )
    return ()
  }

  if ( is.null ( Ct ) ) 
  {
  	Ct <- rep(NA,nT)
  }

  # Populate the vector, with a special case for t = 1
  Bt [ 1 ] <- 2 * Bmsy * epst[1] * zetat[1]

  # Loop over remaining years
  for ( t in 2:nT )
  {
    if ( is.na(Ct[t-1]) ) Ct [ t-1 ] <- Ft[t-1] * Bt [ t-1 ]
    Bt [ t ] <-   (	Bt [ t-1 ] +                       # prev year's B
                  	2. * Fmsy * Bt [ t-1 ] * 
                  	( 1. - Bt [ t-1 ] / Bmsy / 2. ) -  # recruitment
                  	Ct [ t-1 ]                         # Catch
                  )
    ## THIS IS BOGUS ##
    Bt[t] <- max(Bt[t], 1e-1)

    Bt[t] <- Bt[t]*epst[t]*zetat[t] # Process error
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
  outList $ epst <- epst
  outList $ zetat <- zetat
  outList $ dep <- dep
  return ( outList )
}

# obsModel()
# A function which creates observations (indices of biomass) from
# given true series of biomass and observation model parameters.
# inputs:			Bt=true series of biomass; q=survey catchability parameter
#							nT=length of time series; deltat=nT-vector of standard errors
# 						tau=obs error standard deviation (uncorrelated)
# outputs:		It=nT-vector of survey observations
# usage:			creates survey observations which are supplied to the estimator
.obsModel <- function ( Bt, q = 0.05, nT = length (Bt), 
												deltat = rnorm(nT), tau = 0.4 )
{
  deltat <- deltat*tau
  tau2 <- tau*tau
	It <- q * Bt * exp ( deltat - tau2/2.0 )
	outList <- list ( It = It, deltat = deltat )
}


# fillRanWalk()
# A function that will return a bias corrected log-normal AR(1) random 
# process when given a set of deviations (assumed to be standard normal, 
# but not necessarily), uncorr. standard deviation and an autocorrelation 
# factor
# inputs:   z=vector of standard normal deviations; sd=uncorrelated sd
#           rho=autocorrelation factor
# outputs:  zcorr=vector of autocorrelated log-normal deviations
# usage:    to provide auto-correlated process errors
.fillRanWalk <- function ( z=rnorm(10), s = 1, rho = 0, c=0, normScale=FALSE )
{
  # initialise output vector
  zcorr <- numeric ( length ( z ) )

  # scale z by the uncorrelated standard dev
  z <- z * s

  # Create correlated variance term
  s2 <- s * s
  s2corr <- s2 / (1 - rho * rho )

  # loop to populate auto-correlated vector
  zcorr[1] <- c + z[1]
  for ( t in 2:length(z) )
  {
    zcorr[t] <- c + rho * zcorr[t-1] + z[t]
  }

  if (normScale) return (zcorr)

  # Make log-normal
  zcorr <-  exp ( zcorr - s2corr/2 )
  zcorr
}

# genCorrDevs()
# Function that creates sequences of deviations correlated across 
# populations (p) given a correlation matrix, vector of standard deviations 
# for each population and a matrix of uncorrelated standard normals.
# inputs:   z=nxp matrix of uncorrelated standard normals
#           Mcorr=pxp positive-definite correlation mtx
#           Sigma=p-vector of standard deviations
# outputs:  zcorr=nxp matrix of correlated standard normals
# usage:    providing correlations in error across groups/populations
.genCorrDevs <- function ( z = matrix(rnorm(30),nrow=3 ),
                          Mcorr = diag(c(1,1,1)), Sigma = c(1,1,1),
                          normScale=FALSE )
{
  # Find the Cholesky factor of the correlation matrix
  # cov <- diag(Sigma) %*% Mcorr %*% diag(Sigma)
  Mchol <- chol ( Mcorr )

  # create correlated random normals (chol returns UT matrix, so transpose)
  zcorr <- t(Mchol) %*% z
  # Scale by Sigma
  Sigma <- diag ( Sigma )
  zcorr <- Sigma %*% zcorr
  if ( normScale ) return(zcorr)
  for ( k in 1:ncol(zcorr))
  {
    zcorr[,k] <- zcorr[,k] - diag(Sigma*Sigma)/2 
  }
  zcorr <- exp(zcorr)
  zcorr
}

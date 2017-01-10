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
runSimEst <- function ( ctlFile = "simCtlFile.txt", folder=NULL, quiet=TRUE )
{ 
  # read in control file
  controlList <- .readParFile ( ctlFile )
  controlList <- .createList  ( controlList )

  # Run simEst Procedure
  blob <- .simEstProc( obj = controlList, quiet )
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
  kappa <- sqrt(obj$opMod$pars$kappa)    # longitudinal shared proc error sd
  Sigma <- sqrt(obj$opMod$pars$Sigma)    # ms proc error cov mtx
  
  # Rescale shared effects sd if desired
  if (!is.null(obj$opMod$kappaMult)) kappa <- kappa*obj$opMod$kappaMult

  # Correlation parameters
  gamma     <- obj$opMod$gamma          # longitudinal auto-correlation
  msCorr  <- obj$opMod$pars$msCorr  # ms cross correlation

  # Multiply the off-diagonal elements of msCorr if desired
  if (!is.null(obj$opMod$corrMult))
  {
    corrMult <- obj$opMod$corrMult
    msCorr <- corrMult*msCorr
    diag(msCorr) <- 1
  }

  # Obs error var
  tau <- sqrt(obj$opMod$tau2)

  # Initialise list to hold the data
  om <- list (  Bt    = matrix (NA,nrow=nS, ncol=nT ),
                Ct    = matrix (NA,nrow=nS, ncol=nT ),
                Ut    = matrix (NA,nrow=nS, ncol=nT ),
                ItTrue= matrix (NA,nrow=nS, ncol=nT ),
                epst  = numeric (length=nT ),
                zetat = matrix (NA,nrow=nS, ncol=nT ),
                deltat= matrix (NA,nrow=nS, ncol=nT ),
                It    = matrix (NA,nrow=nS, ncol=nT ),
                dep   = numeric(length = nS),
                kappa2= kappa*kappa,
                Sigma2= Sigma*Sigma,
                msCorr= msCorr,
                gamma = gamma,
                tau2  = tau*tau,
                q     = obj$opMod$q,
                nT    = nT,
                nS    = nS,
                Bmsy  = obj$opMod$pars$Bmsy,
                Umsy  = obj$opMod$pars$Umsy )

  # Now create epst and zetat vectors using the proc error components
  epst      <- .fillRanWalk(  z=rnorm(n = nT), s=kappa,
                              gamma=gamma )
  zetat     <- matrix  (rnorm ( n = nS*nT ), nrow = nS, ncol = nT)
  zetat     <- .genCorrDevs ( zetat,
                              Mcorr=msCorr,
                              Sigma=Sigma )

  # standard normal deviations for obs error (uncorrelated)
  deltat    <- matrix(rnorm ( n = nS*nT ), nrow = nS, ncol = nT) 

  # Create Ut time series
  ## REPLACE WITH A FUNCTION LATER FOR MULTIPLE correlated Ft TRAJECTORIES ##
  UtProp <- numeric ( length = nT )
  Umult <- obj$opMod$Umult
  # Then get time of peak fishing mortality
  tUpeak    <- obj$opMod$tUpeak
  tUtrough  <- obj$opMod$tUtrough
  # Compute gradients of the pieces
  Ugrad1 <- ( Umult [ 2 ] - Umult [ 1 ] ) / (tUpeak - 1 )
  Ugrad2 <- ( Umult [ 3 ] - Umult [ 2 ] ) / 
                  ( tUtrough - tUpeak )
  # Populate Ut with multipliers of Umsy
  UtProp [ 1:tUpeak ]                <- Umult [ 1 ] + 0:(tUpeak - 1) * Ugrad1
  UtProp [ (tUpeak + 1):tUtrough  ]  <- Umult [ 2 ] + 1:(tUtrough - tUpeak) * Ugrad2
  UtProp [ (tUtrough+1):nT        ]  <- Umult [ 3 ]

  # Overwrite Ft with actual F values
  Ust <- matrix ( UtProp, nrow = nS, ncol = nT, byrow = TRUE )
  for ( s in 1:nS ) Ust[s,] <- obj$opMod$pars$Umsy[s] * Ust[s,]

  # Loop over species and fill list entries with biological
  # and observational data
  for ( s in 1:nS )
  {
    bio <- .logProdModel (  Bmsy = obj$opMod$pars$Bmsy[s], 
                            Umsy = obj$opMod$pars$Umsy[s],
                            nT = nT, Ut = Ust[s,], 
                            epst = epst,
                            zetat = zetat[s,] )
    obs <- .obsModel (  Bt = bio $ Bt, q = obj$opMod$q[s], nT = nT,
                        deltat = deltat[s,], tau = tau[s] )

    om$Bt[s,]     <- bio$Bt
    om$Ct[s,]     <- bio$Ct
    om$Ut[s,]     <- bio$Ut
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
  maxU <- apply ( X = obj$om$Ut, MARGIN=1,FUN=max)
  maxU <- as.numeric(maxU)

  # loop over species
  for (s in 1:nS )
  {
    # Make dat list
    ssDat[[s]] <- list (  It = matrix(obj$om$It[s,], nrow=1),
                          Ct = matrix(obj$om$Ct[s,], nrow=1) )
    # make par list
    ssPar[[s]] <- list (  lnBmsy            = log(obj$assess$Bmsy[s]),
                          lnUmsy            = log(obj$assess$Umsy[s]),
                          lntau2            = log(obj$assess$tau2),
                          lnq               = rep( -1, 1 ) ,
                          lnqbar            = obj$assess$lnqbar,
                          lntauq2           = obj$assess$lntauq2,
                          mlnq              = obj$assess$mlnq,
                          s2lnq             = obj$assess$s2lnq,
                          lnUmsybar         = obj$assess$lnUmsybar,
                          lnsigUmsy2        = obj$assess$lnsigUmsy2,
                          mlnUmsy           = obj$assess$mlnUmsy,
                          s2lnUmsy          = obj$assess$s2lnUmsy,
                          mlnBmsy           = obj$assess$mlnBmsy[s],
                          s2lnBmsy          = obj$assess$s2lnBmsy[s],
                          eps_t             = rep( 0, nT ),
                          lnkappa2          = log( obj$assess$kappa2 ),
                          logit_gammaYr     = obj$assess$logit_gammaYr,
                          zeta_st           = matrix( 0, nrow = 1, ncol = nT ),
                          lnSigmaDiag       = 0,
                          logitSigmaOffDiag = numeric( length = 0 )
                        )
    ssMap     <-  list(   s2lnq         = factor( NA ),
                          s2lnUmsy      = factor( NA ),
                          lntauq2       = factor( NA ),
                          lnsigUmsy2    = factor( NA ),
                          mlnBmsy       = factor( rep( NA, 1 ) ),
                          s2lnBmsy      = factor( rep( NA, 1 ) ),
                          mlnUmsy       = factor( NA ),
                          mlnq          = factor( NA ),
                          zeta_st       = factor( rep( NA, nT ) ),
                          lnSigmaDiag   = factor( NA)  )
  }
  # now make dat, par and map (par masking) lists for the MS model
  msDat <- list ( It = obj$om$It,
                  Ct = obj$om$Ct
                )

  msPar <- list ( lnBmsy            = log(obj$assess$Bmsy),
                  lnUmsy            = log(obj$assess$Umsy),
                  lntau2            = log(obj$assess$tau2),
                  lnq               = rep(-1, nS),
                  lnqbar            = obj$assess$lnqbar,
                  lntauq2           = obj$assess$lntauq2,
                  mlnq              = obj$assess$mlnq,
                  s2lnq             = obj$assess$s2lnq,
                  lnUmsybar         = obj$assess$lnUmsybar,
                  lnsigUmsy2        = obj$assess$lnsigUmsy2,
                  mlnUmsy           = obj$assess$mlnUmsy,
                  s2lnUmsy          = obj$assess$s2lnUmsy,
                  mlnBmsy           = obj$assess$mlnBmsy,
                  s2lnBmsy          = obj$assess$s2lnBmsy,
                  eps_t             = rep(0, nT),
                  lnkappa2          = log(obj$assess$kappa2),              
                  logit_gammaYr     = obj$assess$logit_gammaYr,
                  zeta_st           = matrix(0, nrow = nS, ncol = nT),
                  lnSigmaDiag       = 0,
                  logitSigmaOffDiag = numeric( length = nS*(nS-1)/2 )
                )

  msMap <- list ( s2lnq       = factor( NA ),
                  s2lnUmsy    = factor( NA ),
                  lntauq2     = factor( NA ),
                  lnsigUmsy2  = factor( NA ),
                  mlnBmsy     = factor( rep( NA, nS ) ),
                  s2lnBmsy    = factor( rep( NA, nS ) ),
                  mlnUmsy     = factor( NA ),
                  mlnq        = factor( NA ) )
  # return list of dat and pin objects for running estimators
  outlist <- list ( ssDat = ssDat, 
                    ssPar = ssPar,
                    ssMap = ssMap, 
                    msDat = msDat, 
                    msPar = msPar,
                    msMap = msMap )
  return(outlist)
}


# callProcTMB()
# Function that calls an ADMB estimator for a given set of
# dat and par list objects.
# inputs:   dat=list object containing .dat file contents in order
#           par=list object containing .pin file contents in order
#           fitTrials=number of times to modify pin before giving up
# ouputs:   rep=list object containing the rep file contents
#           localMin=Logical indicating if local minimum was found
#           hessPosDef=logical indicating if MCMC was possible
# usage:    fitting TMB model to a replicate OM run's catch and CPUE
.callProcTMB <- function (  dat=ssDat[[1]], par=ssPar[[1]], map = ssMap,
                            fitTrials = 3, maxfn = 10000,
                            quiet = TRUE, TMBlib="msProd",
                            RE = c("eps_t","lnq","lnUmsy","zeta_st") )
{ 
  # First, load the dynamic library object
  dyn.load(dynlib("msProd"))
    
  # Make the AD function
  obj <- MakeADFun (  dat = dat, parameters = par, map = map,
                      random = RE, silent = quiet )

  ctrl = list ( eval.max = maxfn, iter.max = maxfn )

  # optimise the model
  fit <- try( nlminb (  start = obj$par,
                        objective = obj$fn,
                        gradient = obj$gr,
                        control = ctrl ) )

  # Run SD report on the fit if it works
  if( class ( fit ) == "try-error" ) fitrep <- NA
  else fitrep <- try( sdreport( obj ) )

  if( class( fitrep ) == "try-error" ) fitrep <- NA

  # Return
  fitrep
}

# seedFit()
# seedFit takes a given seed value and control list, and
# returns a full run of the simulation-estimation procedure using
# the supplied seed values and the control list parameters
# inputs:   seed=integer value provided to the rng
#           obj=ctlList set in controlFile.txt
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

  # Call EMs
  ssFit <- list()
  for (s in 1:obj$opMod$nS ) 
  {
    cat ( "Fitting species ", s, " of ", obj$opMod$nS, "\n", sep = "")
    ssFit[[s]] <- .callProcTMB (  dat = datPar$ssDat[[s]], 
                                  par = datPar$ssPar[[s]],
                                  map = datPar$ssMap,
                                  fitTrials = obj$assess$fitTrials, 
                                  TMBlib = obj$assess$TMBlib,
                                  maxfn = obj$assess$maxfn,
                                  quiet = obj$assess$quiet)
  }
  cat ( "Fitting ", obj$opMod$nS," species hierarchically.\n", 
                    sep = "")
  msFit <- .callProcTMB ( dat = datPar$msDat, 
                          par = datPar$msPar,
                          map = datPar$msMap,
                          fitTrials = obj$assess$fitTrials, 
                          TMBlib = "msProd",
                          maxfn = obj$assess$maxfn,
                          quiet = obj$assess$quiet )

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

  # Create dimension names for arrays
  rNames <- paste("Rep",1:nReps, sep ="")
  sNames <- obj$ctrl$specNames

  # Create list object to store all simulated values
  om <- list  ( Bt    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                Ct    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                Ut    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                epst  = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                zetat = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                deltat= array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                It    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                ItTrue= array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                dep   = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                kappa2= NA,
                Sigma2= NA
              )
                # sigma = numeric(length = nReps),
                # Sigma = numeric(length = nReps),
                # tau   = numeric(length = nReps),
                # q     = numeric(length = nReps),
                # nT    = numeric(length = nReps),
                # nS    = numeric(length = nReps),
                # msy   = numeric(length = nReps),
                # Umsy  = numeric(length = nReps) )

  # Now create a list object to ss and ms assessment model outputs
  am <- list ( ss=NULL, ms=NULL)

  # single species
  am$ss <- list ( Umsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  Bmsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  msy     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  q       = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  kappa2  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  tau2    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  dep     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  epst    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  gamma   = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  Bt      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  Ut      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  mlnq    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  slnq    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  hesspd  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  maxGrad = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)))

  # multispecies (coastwide)
  am$ms <- list ( Umsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  Bmsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  msy     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  q       = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  kappa2  = rep   (NA,length=nReps),
                  Sigma2  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  tau2    = rep   (NA,nReps),
                  dep     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  epst    = matrix(NA,nrow=nReps,ncol=nT,dimnames=list(rNames,1:nT)),
                  zetat   = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  gamma   = vector("numeric",length=nReps),
                  chol    = array (NA,dim=c(nReps,nS,nS),dimnames=list(rNames,sNames,sNames)),
                  Bt      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  Ut      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  mlnq    = vector("numeric",length=nReps),
                  slnq    = vector("numeric",length=nReps),
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
    blob$om$Ut[i,,]         <- simEst$om$Ut
    blob$om$zetat[i,,]      <- simEst$om$zetat
    blob$om$deltat[i,,]     <- simEst$om$Bt
    blob$om$It[i,,]         <- simEst$om$It
    blob$om$ItTrue[i,,]     <- simEst$om$ItTrue
    blob$om$dep[i,]         <- simEst$om$dep
    blob$om$kappa2          <- simEst$om$kappa2
    blob$om$Sigma2          <- simEst$om$Sigma2


    # Save AM results
    # Loop over single species
    for ( s in 1:nS )
    {
      if (  !is.na(simEst$ssFit[[s]]) )
      {
        # Recover report from optimisation
        fitrep    <- simEst$ssFit[[s]]
        estList   <- as.list(fitrep,"Estimate")

        # Now save estimates
        blob$am$ss$Bmsy[i,s]        <- exp(estList$lnBmsy)
        blob$am$ss$Umsy[i,s]        <- exp(estList$lnUmsy)
        blob$am$ss$msy[i,s]         <- fitrep$value["msy"]
        blob$am$ss$q[i,s]           <- exp(estList$lnq)
        blob$am$ss$kappa2[i,s]      <- exp(estList$lnkappa2)
        blob$am$ss$tau2[i,s]        <- exp(estList$lntau2)
        # blob$am$ss$dep[i,s]         <- simEst$ssFit[[s]]$fitrep$D
        blob$am$ss$epst[i,s,]       <- estList$eps_t
        blob$am$ss$gamma[i,s]       <- 2 / ( 1 + exp(-2 * estList$logit_gammaYr)) - 1
        blob$am$ss$Bt[i,s,]         <- fitrep$value[1:nT]
        # blob$am$ss$Ut[i,s,]         <- blob$om$Ct/fitrep$value[1:nT]
        blob$am$ms$qbar[i]          <- exp(estList$lnqbar)
        blob$am$ms$Umsybar[i]       <- exp(estList$lnUmsybar)
        # Estimator performance flags
        blob$am$ss$maxGrad[i,s]     <- max(fitrep$gradient.fixed)
        blob$am$ss$hesspd[i,s]      <- fitrep$pdHess
      }
     
    }

    # Now multispecies
    if (  !is.na(simEst$msFit) )
    {
      # Recover report from optimisation
      fitrep    <- simEst$msFit
      estList   <- as.list(fitrep,"Estimate")
      # Now save estimates
      blob$am$ms$Bmsy[i,]           <- exp(estList$lnBmsy)
      blob$am$ms$Umsy[i,]           <- exp(estList$lnUmsy)
      blob$am$ms$msy[i,]            <- fitrep$value["msy"]
      blob$am$ms$q[i,]              <- exp(estList$lnq)
      blob$am$ms$Sigma2[i,]         <- exp(estList$lnSigmaDiag)
      blob$am$ms$kappa2[i]          <- exp(estList$lnkappa2)
      blob$am$ms$tau2[i]            <- exp(estList$lntau2)
      # blob$am$ms$dep[i,]            <- simEst$msFit$fitrep$D
      blob$am$ms$epst[i,]           <- estList$eps_t
      blob$am$ms$zetat[i,,]         <- estList$zeta_st
      blob$am$ms$gamma[i]           <- 2 / ( 1 + exp(-2 * estList$logit_gammaYr)) - 1
      # blob$am$ms$chol[i,,]          <- simEst$msFit$fitrep$chol
      blob$am$ms$Bt[i,,]            <- matrix(fitrep$value[1:(nS*nT)],nrow=nS,byrow=FALSE)
      # blob$am$ms$Ut[i,,]            <- blob$om$Ut/matrix(fitrep$value[1:(nS*nT)],nrow=nS,byrow=FALSE)
      blob$am$ms$qbar[i]            <- exp(estList$lnqbar)
      blob$am$ms$Umsybar[i]         <- exp(estList$lnUmsybar)
    
      # Estimator Performance Flags
      blob$am$ms$hesspd[i]          <- fitrep$pdHess
      blob$am$ms$maxGrad[i]         <- max(fitrep$gradient.fixed)
    }
  }

  cat ( "Completed ", nReps, " replicates.\n", sep = "" )

  blob
}


# logProdModel()
# Function to simulate a logistic based surplus production model,
# assuming that population is initially at equilibrium (2*Bmsy) modified
# by a proc error.
# inputs:     msy=maximum sustainable yield; Umsy=fishing mortality for msy
#             nT=length of simulation; Ct=nT-vector of catch (opt)
# 						Ft=nT-vector of fishing mortality (opt)
#             epst=nT-vector of env (AR(1)) logN proc errors (bias corrected) 
#             zetat=nT-vector of logN proc errors (correlated among species)
# outputs:    Bt=nT-vector of modeled biomass; Ct=nT-vector of catch
# 						Ut=vector of exploitation rates; 
# usage:      Bexp(eps) = indices of abundance for estimation w/ logN errors eps
.logProdModel <- function ( Bmsy = 1, Umsy = 0.1, nT = 50, Ut = NULL, Ct = NULL,
                            epst = exp(rnorm ( n = nT )-0.5), 
                            zetat = exp(rep(0,nT)) )
{
  # First, initialise a vector to hold biomass
  Bt <- numeric ( length = nT )

  if ( is.null ( Ut ) & is.null ( Ct ) ) 
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
    if ( is.na(Ct[t-1]) ) Ct [ t-1 ] <- Ut[t-1] * Bt [ t-1 ]
    Bt [ t ] <-   (	Bt [ t-1 ] +                       # prev year's B
                  	2. * Umsy * Bt [ t-1 ] * 
                  	( 1. - Bt [ t-1 ] / Bmsy / 2. ) -  # recruitment
                  	Ct [ t-1 ]                         # Catch
                  )
    ## THIS IS BOGUS ##
    Bt[t] <- max(Bt[t], 1e-1)

    Bt[t] <- Bt[t]*epst[t]*zetat[t] # Process error
}  # End loop for biomass
  # Generate catch for nT
  Ct[nT] <- Ut[nT] * Bt[nT]

  # compute depletion value for comparison to model fit
  dep <- Bt [ nT ] / Bmsy / 2.

  # Return B, C, U and U series and
  # depletion value
  outList <- list ()
  outList $ Bt <- Bt
  outList $ Ct <- Ct
  outList $ Ut <- Ut
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
.fillRanWalk <- function ( z=rnorm(10), s = 1, gamma = 0, c=0, normScale=FALSE )
{
  # initialise output vector
  zcorr <- numeric ( length ( z ) )

  # scale z by the uncorrelated standard dev
  z <- z * s

  # Create correlated variance term
  s2 <- s * s
  s2corr <- s2 / (1 - gamma * gamma )

  # loop to populate auto-correlated vector
  zcorr[1] <- c + z[1]
  for ( t in 2:length(z) )
  {
    zcorr[t] <- c + gamma * zcorr[t-1] + z[t]
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

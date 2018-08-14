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
  # Copy control file to sim folder for posterity
  file.copy(from=ctlFile,to=file.path(path,"simCtlFile.txt"))
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

  # Now recover multi-species model dimensions
  fYear <- obj$opMod$fYear        # starting year of catch history
  lYear <- obj$opMod$lYear        # Final year (now) of catch history 
  nT    <- lYear - fYear + 1      # Total length of sim for each species
  sT    <- lYear - min(fYear) + 1 # initial time step of observations
  nS    <- obj$opMod$nS           # number of species

  # Create Ut time series
  ## REPLACE WITH A FUNCTION LATER FOR MULTIPLE correlated Ft TRAJECTORIES ##
  UtProp <- numeric ( length = max(nT) )
  Umult <- obj$opMod$Umult
  if (!is.null(obj$opMod$Umax)) Umult[2] <- obj$opMod$Umax
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
  UtProp [ (tUtrough+1):max(nT)        ]  <- Umult [ 3 ]

  # Overwrite Ft with actual F values
  Ust <- matrix ( UtProp, nrow = nS, ncol = max(nT), byrow = TRUE )
  for ( s in 1:nS ) Ust[s,] <- obj$opMod$Umsy[s] * Ust[s,]

  # Check to see if this is an observation error only simulation.
  if( obj$ctrl$fixProc )
  {
    # If so, check what simulation this is
    if( seed == obj$ctrl$rSeed + 1 )
    { # If this is the first one, then initialise the OM using .popInit, 
      # and save the om to a global object for future simulations
      obj$om <- .popInit( obj = obj, Ust = Ust,
                          specNames = obj$ctrl$speciesNames  )
      globalOM <<- obj$om
    } # If this is rep >= 2 then load the saved OM, and skip the .popInit
    else obj$om <- globalOM 
  } else # If proc errors are changing, then just run the OM every time
      obj$om <- .popInit( obj = obj, Ust = Ust,
                          specNames = obj$ctrl$speciesNames  )



  # Now apply the observation model, this includes the random
  # draws of species specific catchability within each survey,
  # survey specific sampling design and observation error
  obj$om <- .obsModel( obj = obj )

  return(obj)
}

# .popInit()
# Initialises population to achieve target depletion
# or biomass correlation, calls .UdevOptimNLL() to do so
# inputs:   type="corr", "dep" or "off" (character)
#           obj=opMod list object
# ouputs:   om=populated om list object
# usage:    in opModel()
.popInit <- function( obj = blob,
                      Ust = Ust, specNames )
{

  cat("MSG (popInit) Initialising population dynamics\n" )
  nS    <- obj$opMod$nS           # number of species
  fYear <- obj$opMod$fYear[1:nS]  # starting year of catch history
  lYear <- obj$opMod$lYear        # Final year (now) of catch history 
  nT    <- lYear - fYear + 1      # Total length of sim for each species
  sT    <- fYear - min(fYear) + 1 # initial time step of observations  

  # Process error component vars
  kappa <- sqrt(obj$opMod$kappa2)    # longitudinal shared proc error sd
  Sigma <- sqrt(obj$opMod$SigmaDiag[1:nS])    # ms proc error cov mtx
  
  # Rescale shared effects sd if desired
  if (!is.null(obj$opMod$kappaMult)) kappa <- kappa*obj$opMod$kappaMult

  # Correlation parameters
  gammaYr <- obj$opMod$gammaYr        # longitudinal auto-correlation
  # Create correlation matrix
  msCorr <- matrix( obj$opMod$corrOffDiag, nrow = nS, ncol = nS)
  diag(msCorr) <- 1

  # Now create epst and zetat vectors using the proc error components
  epst      <- .fillRanWalk(  z=rnorm(n = max(nT)), s=kappa,
                              gamma=gammaYr )
  zetat     <- matrix  (rnorm ( n = nS*max(nT) ), nrow = nS, ncol = max(nT))
  zetat     <- .genCorrDevs ( zetat,
                              Mcorr=msCorr,
                              Sigma=Sigma )

  # Create initial devations for each case
  if( obj$opMod$depOpt ) targYr  <- obj$opMod$targYr[1:nS] else targYr <- rep(lYear,nS)
  UdevT   <- targYr - min(fYear) + 1
  nDevs   <- sum(UdevT) 
  devs    <- rep(0,nDevs)


  # Initialise list to hold the data
  om <- list (  Nt        = matrix( NA,nrow=nS, ncol=max(nT) ),
                Bt        = matrix( NA,nrow=nS, ncol=max(nT) ),
                wbart     = matrix( NA,nrow=nS, ncol=max(nT) ),
                Ct        = matrix( NA,nrow=nS, ncol=max(nT) ),
                Ut        = matrix( NA,nrow=nS, ncol=max(nT) ),
                epst      = epst,
                zetat     = zetat,
                dep       = matrix( NA,nrow=nS, ncol=max(nT) ),
                kappa2    = kappa*kappa,
                Sigma2    = Sigma*Sigma,
                msCorr    = msCorr,
                gamma     = gamma,
                nT        = nT,
                sT        = sT,
                nS        = nS,
                Bmsy      = obj$opMod$Bmsy[1:nS],
                Umsy      = obj$opMod$Umsy[1:nS],
                specNames = specNames )


  # Transform Us to Fs for optimisation (easier for bounding)
  Fst <- log(1 / (1 - Ust) )

  if( (obj$opMod$corrOpt) | (obj$opMod$depOpt) )
  {
    cat("MSG (popInit) Optimising fishing history for desired dynamics\n" )
    # Optimise deviations
    optObj <- optim( par = devs, fn = .devOptimNLL, devT = UdevT, corr = obj$opMod$corrOpt,
                     obj = obj$opMod, specNames = specNames, baseMtx = Fst, dep = obj$opMod$depOpt,
                     method = "BFGS", om = om ) 
    devs <- optObj$par
    cat("MSG (popInit) Optimisation complete\n" )
  }
  # Populate OM biomass optimised deviations
  om <- .devOptimNLL( devs = devs, devT = UdevT, fit = FALSE, corrOpt = FALSE,
                      dep = FALSE, baseMtx = Fst, om = om,
                      obj = obj$opMod, specNames = specNames  ) 

  return(om)
}

# .UdevOptimNLL()
# Optimises deviations from a provided fishing mortality rate
# matrix (Fst)
# inputs:   devs=vectorized U deviations
#           targYr=year to achieve target, if "corr" targYr=fYear
#           
.devOptimNLL <- function( devs = rep(0,nS*max(nT)),
                          devT = UdevT,
                          baseMtx = Ust,
                          obj = obj$opMod,
                          om = om,
                          specNames,
                          corrOpt = TRUE,
                          depOpt = TRUE,
                          fit = TRUE )
{
  nS    <- om$nS            # number of species
  nT    <- om$nT            # Total length of sim for each species
  sT    <- om$sT            # initial time step of reconstruction 

  # Now
  Fst <- baseMtx
  for(s in 1:nS)
  {
    if( s == 1 ) devIdx <- 1:devT[s] 
    else devIdx <- (sum(devT[1:(s-1)])+1):sum(devT[1:s])
    Fst[s,1:devT[s]] <- baseMtx[s,1:devT[s]] * exp(devs[devIdx])
  }

  Ust <- 1 - exp( - Fst )
  
  # Compute likelihood if optimising
  if( fit )
  {
    estMtx  <- matrix(NA,nrow = nS, ncol = max(nT))
    Bst     <- estMtx
    for( s in 1:nS )
    {
      prodModel <- .logProdModel (  Bmsy = obj$Bmsy[s], 
                                    Umsy = obj$Umsy[s],
                                    nT = max(nT), Ut = Ust[s,], 
                                    epst = om$epst,
                                    zetat = om$zetat[s,],
                                    initDep = obj$initDep[s] )
      estMtx[s,] <- prodModel[[ obj$corrTargVar ]]
      Bst[s,]    <- prodModel$Bt
    }
    nll <- 0
    if( corrOpt )
    { 
      estCorr       <- cor(t(estMtx))
      targMtx       <- om$msCorr
      nll <- nll + 0.5 * sum( (estCorr - targMtx)^2 ) / 0.001
    }

    if( depOpt )
    {
      dep <- numeric(length = nS)
      for(s in 1:nS) dep[s] <- Bst[s,devT[s]] / obj$Bmsy[s] / 2
      nll <- nll + 0.5 * sum( (dep - obj$targDep)^2 ) / 0.001
    }

    # Penalise deviation size
    nll <- nll + 0.5 * sum(devs^2)/0.05
    # Penalise deviations from Umsy
    ssU <- numeric(length = nS)
    for( s in 1:nS ) ssU[s] <- sum((Ust[s,] - obj$Umsy[s])^2)
    nll <- nll + sum(ssU)

    # Penalise crashing the stock - THIS ISN'T FINISHED
    Dst <- Bst
    for( s in 1:nS ) Dst[s,] <- Dst[s,] / obj$Bmsy[s] / 2
    if( any(Dst < 0.05) ) nll <- nll + 1e6

    # Penalise first differences in Fst changes
    for( s in 1:nS ) 
      nll <- nll + 0.5*sum( (Ust[s,1:(max(nT)-1)] - Ust[s,2:max(nT)])^2 ) / 0.1

    return(nll)
  }

  # If already optimised, loop over species and fill list entries with biological
  # and observational data
  for ( s in 1:nS )
  {
    bio <- .logProdModel (  Bmsy = obj$Bmsy[s], 
                            Umsy = obj$Umsy[s],
                            nT = max(nT), Ut = Ust[s,], 
                            epst = om$epst,
                            zetat = om$zetat[s,],
                            initDep = obj$initDep[s] )

    om$Bt[s,]       <- bio$Bt
    om$Ct[s,]       <- bio$Ct
    om$Ut[s,]       <- bio$Ut
    om$dep[s,]      <- bio$dep
  }
  
  return(om)

}


.realDataOM <- function( obj )
{
  # Create dummy om object here for real data
  realDataFile <- file.path( getwd(), paste(obj$ctrl$realData,".Rdata",sep = "") )
  load( realDataFile )
  om <- list()

  # Set time frame
  fYear <- obj$assess$fYear
  lYear <- obj$assess$lYear

  minYear <- min(fYear)
  maxYear <- max(lYear)

  om$I_ost  <- data$indices[,,as.character(minYear:maxYear),1]
  om$Ct     <- data$katch[,as.character(minYear:maxYear)]

  om$sT <- fYear - minYear + 1
  om$nT <- lYear - fYear + 1

  obj$om <- om

  return(obj)
}


# makeDataLists()
# Takes an om list object produced by opModel() and creates 
# data structures that will be passed to estimators
# inputs:   obj=list object containing control and om lists
# outputs:  msDat=multispecies data list; 
#           msPar=multispecies initial parameter value list
#           ssDat=nS-list of single species data lists
#           ssPar=nS-list of single species init. par lists
.makeDataLists <- function ( obj )
{
  
  om    <- obj$om
  opMod <- obj$opMod
  
  # Needs to create nS SS dat and par lists, and 1 MS dat and par list.
  # First, SS:
  # Recover number of species
  nS    <- opMod$nS
  nSurv <- opMod$nSurv
  nT    <- om$nT
  sT    <- om$sT
  # Make dat and par lists
  ssDat <- vector ( mode = "list", length = nS )
  ssPar <- vector ( mode = "list", length = nS )
  ssMap <- vector ( mode = "list", length = nS )

  # Sum catch to get starting biomass
  sumCat <- numeric(length = nS)
  for(s in 1:nS)
  {
    sumCat[s] <- sum( om$Ct[s,sT[s]:max(nT)] )
  }
  sumCat <- as.numeric(sumCat)

  # change IG parameters for tau, kappa and Sigma if autoIG is on
  if (obj$assess$autoIG)
  {
    # recover true variances (use om, may be modified by kappaMult)
    tau2Surv    <- obj$om$tau2Surv
    kappa2      <- obj$om$kappa2
    Sigma2      <- mean(obj$om$Sigma2)

    # Now change IG parameters so that the prior mode is at the true (mean) value
    obj$assess$tau2IGb[1:nSurv] <- (obj$assess$tau2IGa[1:nSurv]+1)*tau2Surv[1:nSurv]
    obj$assess$kappa2IG[2]      <- (obj$assess$kappa2IG[1]+1)*kappa2
    obj$assess$Sigma2IG[2]      <- (obj$assess$Sigma2IG[1]+1)*Sigma2
  }

  # Create IW scale matrix
  if( obj$assess$wishType == "diag" ) wishScale <- diag( opMod$SigmaDiag[1:nS] )
  if( obj$assess$wishType == "corr" )
  {
   # Create correlation matrix - lastNegCorr makes species nS have negative correlation
    msCorr <- matrix (opMod$corrOffDiag, nrow = nS, ncol = nS)
    diag(msCorr) <- 1
    wishScale <- diag(sqrt(opMod$SigmaDiag[1:nS])) %*% msCorr %*% diag(sqrt(opMod$SigmaDiag[1:nS]))
  }

  # loop over species
  for (s in 1:nS )
  {
    # Make dat list 
    ssDat[[s]] <- list (  It              = array(om$I_ost[,s,sT[s]:max(nT)], dim = c(nSurv,1,nT[s])),
                          Ct              = matrix(om$Ct[s,sT[s]:max(nT)], nrow=1),
                          SigmaPriorCode  = obj$assess$SigmaPriorCode,
                          kappaPriorCode  = as.integer(obj$assess$estYearEff),
                          sigUPriorCode   = obj$assess$sigUPriorCode,
                          tauqPriorCode   = obj$assess$tauqPriorCode,
                          lnqPriorCode    = obj$assess$lnqPriorCode,
                          lnUPriorCode    = obj$assess$lnUPriorCode,
                          initT           = 0,
                          initBioCode     = obj$assess$initBioCode[s],
                          posPenFactor    = obj$assess$posPenFactor )
    # make par list
    ssPar[[s]] <- list (  lnBmsy            = log(sumCat[s]),
                          lnUmsy            = log(obj$assess$Umsy[s]),
                          lntau2_o          = log(obj$assess$tau2[1:nSurv]),
                          lnBinit           = log(sumCat[s]/2),
                          lnqbar_o          = rep(obj$assess$lnqbar_o, nSurv),
                          lntauq_o          = rep(obj$assess$lntauq_o, nSurv),
                          mq                = obj$assess$mq,
                          sq                = obj$assess$sq,
                          lnUmsybar         = obj$assess$lnUmsybar,
                          lnsigUmsy         = obj$assess$lnsigUmsy,
                          mUmsy             = obj$assess$mUmsy,
                          sUmsy             = obj$assess$sUmsy,
                          mBmsy             = obj$assess$mBmsy[s],
                          sBmsy             = obj$assess$sBmsy[s],
                          tau2IGa           = obj$assess$tau2IGa[1:nSurv],
                          tau2IGb           = obj$assess$tau2IGb[1:nSurv],
                          tauq2Prior        = obj$assess$tauq2Prior,
                          sigU2Prior        = obj$assess$sigU2Prior,
                          kappa2IG          = obj$assess$kappa2IG,
                          Sigma2IG          = obj$assess$Sigma2IG,
                          wishScale         = matrix(0,nrow=1,ncol=1),
                          nu                = 0,
                          eps_t             = rep( 0.01, nT[s]-1 ),
                          lnkappa2          = ifelse( obj$assess$estYearEff, log( obj$assess$kappa2 ), log( obj$assess$Sigma2 ) ),
                          zeta_st           = matrix( 0.01, nrow = 1, ncol = nT[s]-1 ),
                          lnSigmaDiag       = 0,
                          SigmaDiagMult     = 0,
                          logitSigmaOffDiag = numeric( length = 0 ),
                          logit_gammaYr     = obj$assess$logit_gammaYr
                        )
    # modify SS IG parameters to account for sum of REs
    ssPar[[s]]$kappa2IG[2] <- obj$assess$kappa2IG[2] + obj$assess$Sigma2IG[2]
    
    # Make map list to turn off parameters
    ssMap[[s]]  <-  list( mBmsy             = factor( rep( NA, 1 ) ),
                          sBmsy             = factor( rep( NA, 1 ) ),
                          lnUmsybar         = factor( NA ),
                          lnsigUmsy         = factor( NA ),
                          mUmsy             = factor( NA ),
                          sUmsy             = factor( NA ),
                          lnqbar_o          = factor( rep(NA,nSurv) ),
                          lntauq_o          = factor( rep(NA,nSurv) ),
                          mq                = factor( NA ),
                          sq                = factor( NA ),
                          zeta_st           = factor( rep( NA, nT[s]-1 ) ),
                          lnSigmaDiag       = factor( NA ),
                          SigmaDiagMult     = factor( NA ),
                          tau2IGa           = factor( rep(NA,nSurv) ),
                          tau2IGb           = factor( rep(NA,nSurv) ),
                          sigU2Prior        = factor( rep( NA, 2 ) ),
                          tauq2Prior        = factor( rep( NA, 2 ) ),
                          kappa2IG          = factor( rep( NA, 2 ) ),
                          Sigma2IG          = factor( rep( NA, 2 ) ),
                          wishScale         = factor( NA ),
                          nu                = factor( NA ) )
    if( !obj$assess$ssAR1 ) 
      ssMap[[s]]$logit_gammaYr <- factor( NA )
    if( obj$assess$fixqss )
      ssMap[[s]]$lnq_os <- factor( rep(NA,nSurv) )
    if( obj$assess$initBioCode[s] == 0 )
      ssMap[[s]]$lnBinit <- factor(NA)
  }
  
  
  # now make dat, par and map (par masking) lists for the MS model
  msDat <- list ( It              = om$I_ost,
                  Ct              = om$Ct,
                  SigmaPriorCode  = obj$assess$SigmaPriorCode,
                  kappaPriorCode  = as.integer(obj$assess$estYearEff),
                  sigUPriorCode   = obj$assess$sigUPriorCode,
                  tauqPriorCode   = obj$assess$tauqPriorCode,
                  lnqPriorCode    = obj$assess$lnqPriorCode,
                  lnUPriorCode    = obj$assess$lnUPriorCode,
                  initT           = as.integer(sT - 1)[1:nS],         # correct for indexing starting at 0 in TMB
                  initBioCode     = as.integer(obj$assess$initBioCode[1:nS]),
                  posPenFactor    = obj$assess$posPenFactor
                )
  msPar <- list ( lnBmsy            = log(sumCat[1:nS]),
                  lnUmsy            = log(obj$assess$Umsy[1:nS]),
                  lntau2_o          = log(obj$assess$tau2[1:nSurv]),
                  lnBinit           = log(sumCat[1:nS]/2),
                  lnqbar_o          = rep(obj$assess$lnqbar_o,nSurv),
                  lntauq_o          = rep(obj$assess$lntauq_o,nSurv),
                  mq                = obj$assess$mq,
                  sq                = obj$assess$sq,
                  lnUmsybar         = obj$assess$lnUmsybar,
                  lnsigUmsy         = obj$assess$lnsigUmsy,
                  mUmsy             = obj$assess$mUmsy,
                  sUmsy             = obj$assess$sUmsy,
                  mBmsy             = obj$assess$mBmsy[1:nS],
                  sBmsy             = obj$assess$sBmsy[1:nS],
                  tau2IGa           = obj$assess$tau2IGa[1:nSurv],
                  tau2IGb           = obj$assess$tau2IGb[1:nSurv],
                  tauq2Prior        = obj$assess$tauq2Prior,
                  sigU2Prior        = obj$assess$sigU2Prior,
                  kappa2IG          = obj$assess$kappa2IG,
                  Sigma2IG          = obj$assess$Sigma2IG,
                  wishScale         = wishScale,
                  nu                = nS,
                  eps_t             = rep(0.01, max(nT)-1),
                  lnkappa2          = log(obj$assess$kappa2),              
                  zeta_st           = matrix(0, nrow = nS, ncol = max(nT)-1),
                  lnSigmaDiag       = log(obj$assess$Sigma2),
                  SigmaDiagMult     = obj$assess$SigmaDiagMult[1:nS],
                  logitSigmaOffDiag = rep(0, length = nS*(nS-1)/2),
                  logit_gammaYr     = obj$assess$logit_gammaYr
                )

  lnBinitMap <- 1:nS
  lnBinitMap[ obj$assess$initBioCode[1:nS] == 0] <- NA


  msMap <- list ( mBmsy             = factor( rep( NA, nS ) ),
                  sBmsy             = factor( rep( NA, nS ) ),
                  mUmsy             = factor( NA ),
                  sUmsy             = factor( NA ),
                  mq                = factor( NA ),
                  sq                = factor( NA ),
                  SigmaDiagMult     = factor( rep( NA, nS ) ),
                  tau2IGa           = factor( rep( NA, nSurv ) ),
                  tau2IGb           = factor( rep( NA, nSurv ) ),
                  sigU2Prior        = factor( rep( NA, 2 ) ),
                  tauq2Prior        = factor( rep( NA, 2 ) ),
                  kappa2IG          = factor( rep( NA, 2 ) ),
                  Sigma2IG          = factor( rep( NA, 2 ) ),
                  wishScale         = factor( rep( NA, nS*nS ) ),
                  nu                = factor( NA ),
                  lnBinit           = factor( lnBinitMap )
                )
  # Disable autocorrelation in estimation if set.
  if( !obj$assess$msAR1 ) 
    msMap$logit_gammaYr <- factor( NA )
  if( obj$assess$msCorr == "off" ) 
    msMap$logitSigmaOffDiag <- factor( rep( NA, nS * ( nS - 1 ) / 2 ) )
  if( obj$assess$msCorr == "ident" ) 
    msMap$logitSigmaOffDiag <- factor( rep( 10, nS * ( nS - 1 ) / 2 ) )
  if( obj$assess$fixtauq2 )
    msMap$lntauq_o <- factor( rep(NA,nSurv) )
  if( obj$assess$fixsigU2 )
    msMap$lnsigUmsy <- factor( NA )
  if( obj$assess$fixUbar )
    msMap$lnUmsybar <- factor( NA )
  if( obj$assess$fixqbar )
    msMap$lnqbar_o <- factor( rep(NA, nSurv) )
  if( obj$assess$estqms == "fix" )
    msMap$lnq_os <- factor( rep( NA, nS*nSurv ) )
  if( obj$assess$estqms == "ident" )
    msMap$lnq_os <- factor( rep( 11, nS*nSurv ) )
  if( !obj$assess$estYearEff )
  {
    msMap$lnkappa2 <- factor( NA )
    msMap$eps_t <- factor( rep( NA, max(nT) - 1 ) )
  }

  # browser()
  # return list of dat and pin objects for running estimators
  outlist <- list ( 
                    ssDat = ssDat, 
                    ssPar = ssPar,
                    ssMap = ssMap,
                    msDat = msDat, 
                    msPar = msPar,
                    msMap = msMap 
                  )
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
.callProcTMB <- function (  dat=ssDat[[1]], par=ssPar[[1]], map = ssMap[[1]],
                            fitTrials = 3, maxfn = 10000, quiet = TRUE, 
                            TMBlib="msProd", UB = ssUB, LB = ssLB,
                            RE = c("eps_t"), tracePar = FALSE,
                            profiles = FALSE )
{ 
  # Make the AD function
  obj <- MakeADFun (  dat = dat, parameters = par, map = map,
                      random = NULL, silent = quiet )
  # Trace fixed effect parameter values
  obj$env$tracepar <- tracePar

  # Set max no of evaluations
  ctrl = list ( eval.max = maxfn, iter.max = maxfn )
  
  # Start counting estimation attempts
  nTries <- 1
  
  rep <- obj$report()
  zSum_os <- rep$zSum_os
  validObs <- rep$validObs
  tauq2_o <- rep$tauq2_o
  qbar_o <- rep$qbar_o
  tau2_o <- rep$tau2_o

  nO <- length(tau2_o)
  nS <- length(par$lnBmsy)

  lnqhat_os <- matrix(NA, nrow = nO, ncol = nS )

  for( oIdx in 1:nO )
    for( sIdx in 1:nS )
    {
      if( nS == 1)
        lnqhat_os[oIdx, sIdx] <- zSum_os[oIdx/sIdx] / validObs[oIdx,sIdx]
      else
      {
        lnqhat_os[oIdx,sIdx] <- ( zSum_os[oIdx,sIdx] / tau2_o[oIdx] + log(qbar_o[oIdx])/tauq2_o[oIdx] ) / ( validObs[oIdx,sIdx] / tau2_o[oIdx] + 1 / tauq2_o[oIdx] )
      }

    }


  # optimise the model with all fixed effects
  fitFE <- try( nlminb( start = obj$par + rnorm(length(obj$par),sd = 0.1),
                        objective = obj$fn,
                        gradient = obj$gr,
                        control = ctrl,
                        lower = -Inf,
                        upper = Inf ) )


  convFlag <- fitFE$convergence

  # loop to keep trying to fit if convergence doesn't happen
  while( convFlag > 0 )
  {
    # Increment counter
    nTries <- nTries + 1

    # Run optimisation from last point in par space
    fitFE <- try( nlminb( start = fitFE$par,
                          objective = obj$fn,
                          gradient = obj$gr,
                          control = ctrl,
                          lower = -Inf,
                          upper= Inf ) )

    convFlag <- fitFE$convergence
    if( nTries >= fitTrials )
      break
  }

  nTries <- 1
  bestPars <- fitFE$par

  if( length(RE) > 0 )
  {
    # Add REs
    obj <- MakeADFun (  dat = dat, parameters = par, map = map,
                        random = RE, silent = quiet )  

    # Trace fixed effect parameter values
    obj$env$tracepar <- tracePar

    # Save the best parameters from the fixed eff model
    bestPars <- bestPars[names(obj$par)]

    # optimise the model with REs
    fit <- try( nlminb( start = bestPars + rnorm(length(bestPars), sd = 0.1),
                        objective = obj$fn,
                        gradient = obj$gr,
                        control = ctrl,
                        lower = -Inf,
                        upper= Inf ) )


  } else 
    fit <- fitFE

  # Now try to improve the fit
  while( class(fit) == "try-error" )
  {
    nTries <- nTries + 1
    cat("Optimisation stopped by error, jittering initial parameter values \n" )

    fit <- try( nlminb( start = bestPars + rnorm(length(bestPars), sd = 0.2),
                        objective = obj$fn,
                        gradient = obj$gr,
                        control = ctrl,
                        lower = -Inf,
                        upper= Inf ) )
    # Break out if we reach the number of fit trials
    if( nTries >= 2*fitTrials )
      break
  }

  # Check convergence flag
  convFlag  <- fit$convergence
  # Save best version of parameters so far
  bestPars  <- fit$par
  objFunVal <- fit$objective
  
  # Now try to improve if we get false convergence by jittering the
  # best parameters
  while( convFlag > 0 )
  {
    
    nTries <- nTries + 1

    cat("False converge with given starting values, jittering previous best parameter values \n" )

    fit <- try( nlminb (  start = bestPars + rnorm(length(bestPars), sd = 0.4),
                          objective = obj$fn,
                          gradient = obj$gr,
                          control = ctrl,
                          lower = -Inf,
                          upper= Inf ) )

    if( class(fit) != "try-error" )
    {
      # Update convergence flag
      convFlag <- fit$convergence
      # Update best pars if obj fun value is lower
      if( fit$objective < objFunVal )
        bestPars <- fit$par
    }

    # Break out if we reach the number of fit trials
    if( nTries >= fitTrials )
      break
  }


  # Run SD report on the fit if it works, and get the rep file
  if( class ( fit ) == "try-error" ) 
  {
    sdrep <- NA
    CIs <- NA
    rep <- NA
  } else {
    sdrep <- try( sdreport( obj ) )
    rep   <- try( obj$report() )
  } 
  if( class( rep ) == "try-error" )
  {
    rep <- NA
    CIs <- NA
  }
  if( class( sdrep ) == "try-error" )
  {
    sdrep <- NA
    CIs <- NA
  }
  # Now we need to check the hessian - is it PD? If 
  # not, we need to rejitter and try again...


  # If sdrep exists, compute confidence intervals
  if( !any(is.na( sdrep ) ) )
  {
    CIs <-  summary( sdrep ) 
    colnames( CIs ) <- c( "val", "se" )
    CIs <-  CIs %>% 
              as.data.frame() %>%
              mutate( par = rownames( CIs ),
                      lCI = val - 1.96 * se,
                      uCI = val + 1.96 * se ) %>%
              dplyr::select( par, val, se, lCI, uCI )
    if ( profiles )
    {
      if( nS == 1 )
      {
        # Calculate SS profiles
        lnBmsyProfile   <- tmbprofile( obj, "lnBmsy", trace = FALSE )
        lnUmsyProfile   <- tmbprofile( obj, "lnUmsy", trace = FALSE )
        lnkappa2Profile <- tmbprofile( obj, "lnkappa2", trace = FALSE )
        # Save to an output list
        profileList <- list ( lnBmsy    = lnBmsyProfile,
                              lnUmsy    = lnUmsyProfile,
                              lnkappa2  = lnkappa2Profile )
      }
      if( nS > 1  )
      {
        # Calculate profiles for shared prior hyperpars
        lnUmsybarProfile    <- tmbprofile( obj, "lnUmsybar", trace = FALSE )
        lnsigUmsy2Profile   <- tmbprofile( obj, "lnsigUmsy", trace = FALSE ) 
        lnkappa2Profile     <- tmbprofile( obj, "lnkappa2", trace = FALSE )
        # Save to output list
        profileList <- list ( lnUmsybar   = lnUmsybarProfile,
                              lnsigUmsy2  = lnsigUmsy2Profile,
                              lnkappa2    = lnkappa2Profile )
      }
    }
  } else profileList <- NA

  browser()

  # Return
  out <- list( sdrep = sdrep, rep=rep, CIs = CIs, nTries = nTries )
  if( profiles ) out$profiles <- profileList

  out
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
  if( is.null(obj$ctrl$realData) ) obj <- .opModel ( obj, seed = seed )
  else obj <- .realDataOM( obj )


  # Create data objects for AMs
  datPar <- .makeDataLists ( obj )

  # Call EMs
  ssFit <- list()
  for (s in 1:obj$opMod$nS ) 
  {
    cat ( "Fitting species ", s, " of ", obj$opMod$nS, "\n", sep = "" )
    ssFit[[ s ]] <- .callProcTMB (  dat = datPar$ssDat[[ s ]], 
                                    par = datPar$ssPar[[ s ]],
                                    map = datPar$ssMap[[ s ]],
                                    fitTrials = obj$assess$fitTrials, 
                                    TMBlib = obj$assess$TMBlib,
                                    maxfn = obj$assess$maxfn,
                                    quiet = obj$assess$quiet,
                                    RE = obj$assess$ssRE,
                                    LB = datPar$ssLB[[ s ]],
                                    UB = datPar$ssUB[[ s ]],
                                    tracePar = obj$assess$tracePar,
                                    profiles = obj$assess$profiles )
  }
  cat ( "Fitting ", obj$opMod$nS," species hierarchically.\n", 
                    sep = "")
  msFit <- .callProcTMB ( dat = datPar$msDat, 
                          par = datPar$msPar,
                          map = datPar$msMap,
                          fitTrials = obj$assess$fitTrials, 
                          TMBlib = "msProd",
                          maxfn = obj$assess$maxfn,
                          quiet = obj$assess$quiet,
                          RE = obj$assess$msRE,
                          LB = datPar$msLB,
                          UB = datPar$msUB,
                          tracePar = obj$assess$tracePar,
                          profiles = obj$assess$profiles )

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
  nSurv   <- obj$opMod$nSurv
  fYear   <- obj$opMod$fYear
  lYear   <- obj$opMod$lYear
  nT      <- max(lYear - fYear + 1)
  sT      <- fYear - min(fYear) + 1

  # Create dimension names for arrays
  rNames <- paste("Rep",1:nReps, sep ="")
  sNames <- obj$ctrl$specNames[1:nS]
  survNames <- paste("Survey",1:nSurv,sep = "" )

  # Create list object to store all simulated values
  om <- list  ( Bt        = array( NA, dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                Ct        = array( NA, dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                Ut        = array( NA, dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                epst      = array( NA, dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                zetat     = array( NA, dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                delta_ost = array( NA, dim=c(nReps,nSurv,nS,nT),dimnames=list(rNames,survNames,sNames,1:nT)),
                I_ost     = array( NA, dim=c(nReps,nSurv,nS,nT),dimnames=list(rNames,survNames,sNames,1:nT)),
                dep       = array( NA, dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                q_os      = array( NA, dim=c(nReps,nSurv,nS),dimnames=list(rNames,survNames,sNames))
              )

  # Now create a list object to ss and ms assessment model outputs
  am <- list ( ss=NULL, ms=NULL)

  # single species
  am$ss <- list ( Umsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  Bmsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  msy     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  q_os    = array( NA,dim=c(nReps,nSurv,nS),dimnames=list(rNames,survNames,sNames)),
                  kappa2  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  tau2_o  = array( NA,dim=c(nReps,nSurv,nS),dimnames=list(rNames,survNames,sNames)),
                  dep     = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  epst    = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  gamma   = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  Bt      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  Ut      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  hesspd  = matrix(FALSE,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  maxGrad = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  nTries  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  sdrep   = vector( mode = "list", length = nReps ),
                  CIs     = vector( mode = "list", length = nReps ),
                  fitrep  = vector( mode = "list", length = nReps ) )

  # add a list for likelihood profiles
  if( obj$assess$profiles ) 
    am$ss$profiles <- vector( mode = "list", length = nReps )

  # multispecies (coastwide)
  am$ms <- list ( Umsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  Bmsy    = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  msy     = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  q_os    = array( NA,dim=c(nReps,nSurv,nS),dimnames=list(rNames,survNames,sNames)),                  
                  kappa2  = matrix(NA,nrow=nReps,ncol=1, dimnames=list(rNames,NA)),
                  Sigma2  = matrix(NA,nrow=nReps,ncol=nS,dimnames=list(rNames,sNames)),
                  tau2_o  = array( NA,dim=c(nReps,nSurv),dimnames=list(rNames,survNames)),
                  dep     = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  epst    = matrix(NA,nrow=nReps,ncol=nT,dimnames=list(rNames,1:nT)),
                  zetat   = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  gamma   = matrix(NA,nrow=nReps,ncol=1),
                  SigmaCor= array (NA,dim=c(nReps,nS,nS),dimnames=list(rNames,sNames,sNames)),
                  qbar_o  = array( NA,dim=c(nReps,nSurv),dimnames=list(rNames,survNames)),
                  tauq2_o = array( NA,dim=c(nReps,nSurv),dimnames=list(rNames,survNames)),
                  Umsybar = matrix(NA,nrow=nReps,ncol=1, dimnames=list(rNames,NA)),
                  sigU2   = matrix(NA,nrow=nReps,ncol=1, dimnames=list(rNames,NA)),
                  Bt      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  Ut      = array (NA,dim=c(nReps,nS,nT),dimnames=list(rNames,sNames,1:nT)),
                  hesspd  = matrix(FALSE, nrow = nReps, ncol = 1),
                  maxGrad = matrix(NA, nrow = nReps, ncol = 1),
                  nTries  = matrix(NA, nrow = nReps, ncol = 1),
                  sdrep   = vector( mode = "list", length = nReps ),
                  CIs     = vector( mode = "list", length = nReps ),
                  fitrep  = vector( mode = "list", length = nReps ) )
  
  # add a list for likelihood profiles
  if( obj$assess$profiles ) 
    am$ms$profiles <- vector( mode = "list", length = nReps )
  
  # The BLOOOOOOBBBBBB
  blob <- obj
  blob$om <- om
  blob$am <- am

  tau2mult      <- blob$assess$tau2mult
  SigmaDiagMult <- blob$assess$SigmaDiagMult

  # Now create a vector of seed values, able to be lapplied over
  seeds <- rSeed + 1:nReps
  
  cat ( "Fitting ", nReps, " replicates.\n", sep = "" )

  for ( i in 1:length(seeds) )
  {
    # Run sim-est procedure
    simEst <- .seedFit ( seed = seeds[i], obj=obj, quiet = quiet)


    # Save OM values
    # If simulating, save all biological OM values
    if( is.null( obj$ctrl$realData ) )
    {
      blob$om$Bt[i,,]         <- simEst$om$Bt
      blob$om$epst[i,,]       <- simEst$om$epst
      blob$om$Ut[i,,]         <- simEst$om$Ut
      blob$om$zetat[i,,]      <- simEst$om$zetat
      blob$om$delta_ost[i,,,] <- simEst$om$delta_ost
      blob$om$q_os[i,,]       <- simEst$om$q_os
      blob$om$dep[i,,]        <- simEst$om$dep  
    }
    # Save catch and indices if not simulating
    blob$om$Ct[i,,]         <- simEst$om$Ct
    blob$om$I_ost[i,,,]     <- simEst$om$I_ost

    # Save AM results
    # Loop over single species
    # blob$am$ss$rep <- vector(mode = "list", length = nS)
    blob$am$ss$sdrep[[i]]       <- vector( mode = "list", length = nS )
    blob$am$ss$fitrep[[i]]      <- vector( mode = "list", length = nS )
    blob$am$ss$CIs[[i]]         <- vector( mode = "list", length = nS )
    if( obj$assess$profiles ) 
      blob$am$ss$profiles[[i]] <- vector( mode = "list", length = nS )
    for( s in 1:nS )
    {
      if(  !any(is.na(simEst$ssFit[[s]]$rep)) & !any(is.null(simEst$ssFit[[s]]$rep)) )
      {
        # Recover report from optimisation
        sdrep     <- simEst$ssFit[[s]]$sdrep
        fitrep    <- simEst$ssFit[[s]]$rep
        CIs       <- simEst$ssFit[[s]]$CIs

        # Now save estimates
        blob$am$ss$Bmsy[i,s]                  <- fitrep$Bmsy
        blob$am$ss$Umsy[i,s]                  <- fitrep$Umsy
        blob$am$ss$msy[i,s]                   <- fitrep$MSY
        blob$am$ss$q_os[i,,s]                 <- fitrep$qhat_os
        blob$am$ss$kappa2[i,s]                <- fitrep$kappa2
        blob$am$ss$tau2_o[i,,s]               <- fitrep$tau2_o
        blob$am$ss$dep[i,s,sT[s]:nT]          <- fitrep$Bt/fitrep$Bmsy/2
        blob$am$ss$gamma[i,s]                 <- fitrep$gammaYr
        blob$am$ss$Bt[i,s,sT[s]:nT]           <- fitrep$Bt
        # Estimator performance flags
        if(!any(is.na(sdrep)))
        {
          blob$am$ss$maxGrad[i,s]               <- max(sdrep$gradient.fixed)
          blob$am$ss$hesspd[i,s]                <- sdrep$pdHess
          blob$am$ss$sdrep[[i]][[s]]            <- sdrep
        }
        blob$am$ss$CIs[[i]][[s]]              <- CIs
        blob$am$ss$fitrep[[i]][[s]]           <- fitrep
        blob$am$ss$nTries[i,s]                <- simEst$ssFit[[s]]$nTries              

        if( obj$assess$profiles )
          blob$am$ss$profiles[[i]][[s]] <- simEst$ssFit[[s]]$profiles
      }
     
    }

    # Now multispecies
    if (  !any(is.na(simEst$msFit$rep)) & !any(is.null(simEst$msFit$rep))  )
    {
      # Recover report from optimisation
      sdrep     <- simEst$msFit$sdrep
      fitrep    <- simEst$msFit$rep
      CIs       <- simEst$msFit$CIs
      
      # Now save estimates
      blob$am$ms$Bmsy[i,]                 <- fitrep$Bmsy
      blob$am$ms$Umsy[i,]                 <- fitrep$Umsy
      blob$am$ms$msy[i,]                  <- fitrep$MSY
      blob$am$ms$q_os[i,,]                <- fitrep$qhat_os
      blob$am$ms$Sigma2[i,]               <- fitrep$SigmaDiag
      blob$am$ms$kappa2[i,]               <- fitrep$kappa2
      blob$am$ms$Bt[i,,]                  <- fitrep$Bt
      blob$am$ms$tau2_o[i,]               <- fitrep$tau2_o
      blob$am$ms$gamma[i]                 <- fitrep$gammaYr
      blob$am$ms$qbar_o[i,]               <- fitrep$qbar_o
      blob$am$ms$tauq2_o[i,]              <- fitrep$tauq2_o
      blob$am$ms$Umsybar[i,]              <- fitrep$Umsybar
      blob$am$ms$sigU2[i,]                <- fitrep$sigUmsy2
      blob$am$ms$SigmaCor[i,,]            <- fitrep$SigmaCorr
      # Some require looping over species
      for (s in 1:nS)
      { 
        blob$am$ms$dep[i,s,sT[s]:nT]        <- blob$am$ms$Bt[i,s,sT[s]:nT]/fitrep$Bmsy[s]/2
      }
      # Estimator Performance Flags
      if(!any(is.na(sdrep)))
      {
        blob$am$ms$hesspd[i]             <- sdrep$pdHess
        blob$am$ms$maxGrad[i]            <- max(sdrep$gradient.fixed)
        blob$am$ms$sdrep[[i]]            <- sdrep
      }
      blob$am$ms$fitrep[[i]]        <- fitrep
      blob$am$ms$CIs[[i]]           <- CIs
      blob$am$ms$nTries[i]          <- simEst$msFit$nTries
      
      if( obj$assess$profiles )
        blob$am$ms$profiles[[i]] <- simEst$msFit$profiles
    }

    # First get the replicate numbers for succesful fits (MCMC runs) in BOTH models
    ssHess <- blob$am$ss$hesspd
    msHess <- blob$am$ms$hesspd
    hessPD <- ssHess
    for( s in 1:ncol(hessPD) )
      hessPD[,s] <- as.logical(ssHess[,s] * msHess)

    completed <- apply(X = hessPD, FUN = sum, MARGIN = 2, na.rm = T )
    if( all( completed >= obj$ctrl$signifReps ) )
    {
      cat("Successfuly completed ", obj$ctrl$signifReps, " replicates for each stock, ending simulation.\n" )
      break
    } 
  }
  if(any( completed < obj$ctrl$signifReps ) )
    cat (   "Completed ", nReps, " replicates without reaching ", 
            obj$ctrl$signifReps, "for each stock.\n", sep = "" )

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
                            zetat = exp(rep(0,nT)), initDep = 1 )
{
  # First, initialise a vector to hold biomass
  Bt <- numeric ( length = nT )

  if ( is.null ( Ut ) & is.null ( Ct ) ) 
  {
    cat ( "ERR (logProdModel) You're missing harvest, stupid.\n" )
    return ()
  }

  if ( is.null ( Ct ) ) 
  {
  	Ct <- rep(NA,nT)
  }

  # Populate the vector, assuming originally at B0
  Bt [ 1 ] <- 2 * Bmsy * initDep

  # Loop over remaining years
  for ( t in 2:nT )
  {
    if ( is.na(Ct[t-1]) ) Ct [ t-1 ] <- Ut[t-1] * Bt [ t-1 ]
    Bt [ t ] <-   (	Bt [ t-1 ] +                       # prev year's B
                  	2. * Umsy * Bt [ t-1 ] * 
                  	( 1. - Bt [ t-1 ] / Bmsy / 2. ) -  # recruitment
                  	Ct [ t-1 ]                         # Catch
                  )

    Bt[t] <- Bt[t]*epst[t]*zetat[t] # Process error
    ## THIS IS BOGUS ##
    Bt[t] <- max(Bt[t], Ut[t]*Bt[t] + 1)
  }  # End loop for biomass
  # Generate catch for nT
  Ct[nT] <- Ut[nT] * Bt[nT]

  # compute depletion value for comparison to model fit
  dep <- Bt / Bmsy / 2

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
# Creates observations (indices of biomass) for the survey(s)
# given true series of biomass and observation model parameters.
# inputs:     Bt=true series of biomass; q=survey catchability parameter
#             nT=length of time series; deltat=nT-vector of standard errors
#             sT=time step of initial observations
#             tau=obs error standard deviation (uncorrelated)
# outputs:    It=nT-vector of survey observations
# usage:      creates survey observations which are supplied to the estimator
.obsModel <- function (  obj = obj )
{
  cat("MSG (obsModel) Generating random observations\n" )
  om        <- obj$om
  opMod     <- obj$opMod

  # Recover model dimensions
  nT        <- om$nT            # Number of time steps (for each species)
  nS        <- om$nS            # No. of species
  sT        <- om$sT            # time step of start of reconstruction (nS-numeric)
  nSurv     <- opMod$nSurv      # No. of surveys
  fYearSurv <- opMod$fYearSurv  # first year of each survey
  lYearSurv <- opMod$lYearSurv  # last year of each survey
  fYear     <- opMod$fYear      # first year of reconstruction (5-numeric)
  lYear     <- opMod$lYear      # last year of reconstruction (1-numeric)
  survFreq  <- opMod$surveyFreq # Frequency of surveys

  # Derive other dimensions from above
  sTSurv    <- fYearSurv - min(fYear) + 1
  fTSurv    <- lYearSurv - min(fYear) + 1 

  # Recover biomass, obs error variance
  Bst       <- om$Bt
  tauSurv   <- opMod$tauSurv
  tau2Surv  <- tauSurv^2

  # Now create a matrix of q values for each species and survey
  q_os      <- matrix( 0, nrow = nSurv, ncol = nS )
  qSurvM    <- opMod$qSurvM
  qSurvSD   <- opMod$qSurvSD
  for(survIdx in 1:nSurv)
  {
    q_os[ survIdx, ] <- exp( rnorm( n = nS, 
                                    mean = log( qSurvM[ survIdx ] ), 
                                    sd = qSurvSD[ survIdx ] ) )
  }

  # standard normal deviations for obs error (uncorrelated)
  delta_ost    <- array( rnorm (  n = nS*max(nT)*nSurv ), 
                                  dim = c( nSurv, nS, max(nT) ),
                                  dimnames = list(  paste("surv",1:nSurv,sep = "" ),
                                                    obj$ctrl$specNames[1:nS],
                                                    paste("t",1:max(nT), sep = "" ) ) )
  # array to hold observations
  I_ost <- array( -1, dim = c( nSurv, nS, max(nT) ),
                      dimnames = list(  paste("surv",1:nSurv,sep = "" ),
                                        obj$ctrl$specNames[1:nS],
                                        paste("t",1:max(nT), sep = "" ) ) )

  for( survIdx in 1:nSurv )
  {
    tIdxSurv <- seq( from = sTSurv[survIdx], to = fTSurv[survIdx], by = survFreq[survIdx] )
    # Scale standard normal deviations by their SD
    delta_ost[ survIdx, , ] <- delta_ost[ survIdx, , ] * tauSurv[ survIdx ]
    # Now start applying the observation model
    for( s in 1:nS )
    {
      tIdx <- sT[s]:nT[s]
      I_ost[ survIdx, s, tIdxSurv ] <-  q_os[ survIdx, s ] * 
                                        Bst[ s, tIdxSurv ] * 
                                        exp ( delta_ost[ survIdx, s, tIdxSurv ] - 
                                              .5*tau2Surv[ survIdx ] )
      delta_ost[ survIdx, s, - tIdxSurv ] <- NA
    }  
  }
  om$I_ost      <- I_ost
  om$delta_ost  <- delta_ost
  om$q_os       <- q_os
  om$nSurv      <- nSurv
  om$tau2Surv   <- tau2Surv

  om
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

# --------------------------------------------------------------------------
# plots.R
# 
# Plots functions script for the coastwide (single stock) hierarchical
# multispecies assessment model simulation-estimation experiments.
# 
# Author: Samuel D N Johnson
# Date: 8 September, 2016
#
# --------------------------------------------------------------------------

# plotMCMCpar()
# Function that will plot MCMC output from ADMB models for a nominated
# parameter, simulation and replicate. Uses the coda package.
# inputs:   rep=replicate number; 
#           sim=number indicating blob to load from project dir
#           par=character name of parameter
# output:   NULL
# usage:    post-sim, showing MCMC performance
plotMCMCpar <- function ( rep=1, sim=1, par="BMSY" )
{
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path="./project",full.names = FALSE,
                        recursive=FALSE)
  simList <- dirList[grep(pattern="sim",x=dirList)]
  folder <- simList[sim]

  # Load the nominated blob
  blobFileName <- paste(folder,".RData",sep="")
  blobPath <- file.path(getwd(),"project",folder,blobFileName)
  load ( file = blobPath )

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctl$scenario
  mp        <- blob$ctl$mp
  stamp     <- paste(scenario,":",mp,sep="")
  repCount  <- paste("Replicate ",rep,"/",blob$ctl$nReps,sep="")

  # Recover blob elements for plotting
  nS    <- blob$ctl$nS
  nT    <- blob$ctl$nT

  # Species names
  specNames <- blob$ctl$specNames

  # Recover MCMC output for each model and species
  parMCoutSS <- t(blob$am$ss$mcPar[rep,,,as.character(par)])
  parMCoutMS <- t(blob$am$ms$mcPar[rep,,,as.character(par)])
  # There's an issue with having fewer MC replicates than expected
  # so reduce the length of mc output 
  nonNAss <- numeric(nS)
  for (s in 1:nS)
  { 
    nonNAss[s] <- length ( parMCoutSS[,s][!is.na(parMCoutSS[,s])])
  }
  nonNAms <- length(parMCoutMS[!is.na(parMCoutMS)])/nS
  # Set up mc output thinning
  parMCoutSS <- as.mcmc(parMCoutSS[1:min(nonNAss),])
  parMCoutMS <- as.mcmc(parMCoutMS[1:nonNAms,])

  # Now to plot each model's posterior distributions
  titleSS <- paste ( "SS MCMC for ", par, sep = "")
  titleMS <- paste ( "MS MCMC for ", par, sep = "")

  # plot SS trace and posteriors
  plot(parMCoutSS)
  mtext(side=3,text=titleSS,outer=TRUE,padj=2)

  # create a new device and plot MS trace and posteriors
  dev.new()
  plot(parMCoutMS)
  mtext(side=3,text=titleMS,outer=TRUE,padj=2)  

  return()
}

# plotMCMCspecies()
# Function that will plot MCMC output from ADMB models for a nominated
# species, simulation and replicate. Uses the coda package.
# inputs:   rep=replicate number; 
#           sim=number indicating blob to load from project dir
#           spec=number indicating the species
# output:   NULL
# usage:    post-sim, showing MCMC performance
plotMCMCspecies <- function ( rep=1, sim=1, spec=1 )
{
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path="./project",full.names = FALSE,
                        recursive=FALSE)
  simList <- dirList[grep(pattern="sim",x=dirList)]
  folder <- simList[sim]

  # Load the nominated blob
  blobFileName <- paste(folder,".RData",sep="")
  blobPath <- file.path(getwd(),"project",folder,blobFileName)
  load ( file = blobPath )

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctl$scenario
  mp        <- blob$ctl$mp
  stamp     <- paste(scenario,":",mp,sep="")
  repCount  <- paste("Replicate ",rep,"/",blob$ctl$nReps,sep="")

  # Recover blob elements for plotting
  nS    <- blob$ctl$nS
  nT    <- blob$ctl$nT

  # Species names
  if (!is.null(blob$ctl$specNames))
  {
    specNames <- blob$ctl$specNames
    spec <- specNames[spec]
  }
  
  # Recover MCMC output for each model and species
  parMCoutSS <- blob$am$ss$mcPar[rep,spec,,]
  parMCoutMS <- blob$am$ms$mcPar[rep,spec,,]
  # There's an issue with having fewer MC replicates than expected
  # so reduce the length of mc output 
  nonNAss <- length(parMCoutSS[!is.na(parMCoutSS)])/ncol(parMCoutSS)
  nonNAms <- length(parMCoutMS[!is.na(parMCoutMS)])/ncol(parMCoutMS)
  # Set up mc output thinning
  parMCoutSS <- as.mcmc(parMCoutSS[nonNAss,])
  parMCoutMS <- as.mcmc(parMCoutMS[nonNAms,])

  # Now to plot each model's posterior distributions
  titleSS <- paste ( "SS MCMC for ", spec, sep = "")
  titleMS <- paste ( "MS MCMC for ", spec, sep = "")

  # plot SS trace and posteriors
  plot(parMCoutSS)
  mtext(side=3,text=titleSS,outer=TRUE,padj=2)

  # create a new device and plot MS trace and posteriors
  dev.new()
  plot(parMCoutMS)
  mtext(side=3,text=titleMS,outer=TRUE,padj=2)  

  return()
}


# plotMCMCbio()
# Function that plots true and posterior biomass for each species for a
# given simulation and replicate.
# inputs:   rep=replicate number
#           sim=number indicating blob to load from project dir
#           quant=numeric of percentiles to be calculated
plotMCMCbio <- function ( sim = 1, rep = 1, quant=c(0.025,0.5,0.975))
{
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path="./project",full.names = FALSE,
                        recursive=FALSE)
  simList <- dirList[grep(pattern="sim",x=dirList)]
  folder <- simList[sim]

  # Load the nominated blob
  blobFileName <- paste(folder,".RData",sep="")
  blobPath <- file.path(getwd(),"project",folder,blobFileName)
  load ( file = blobPath )

  # Recover control pars
  nS <- blob$ctl$nS
  nT <- blob$ctl$nT


  # load biomass trajectories
  omBio <- blob$om$Bt[rep,,]
  ssBio <- blob$am$ss$mcBio[rep,,,]
  msBio <- blob$am$ms$mcBio[rep,,,]

  # Compute # of MCMC samples
  nSamp <- length(ssBio[1,,1])

  # reduce the biomass to provided quantiles
  ssBio <- apply(X=ssBio,FUN=quantile,MARGIN=c(1,3),na.rm=TRUE,probs=quant)
  msBio <- apply(X=msBio,FUN=quantile,MARGIN=c(1,3),na.rm=TRUE,probs=quant)

  # Now plot!
  par(mfrow = c(3,1), mar = c(1,1,1,1), oma=c(1,1,1,1))
  for (s in 1:nS)
  {
    # Compute max B value
    yLim <- 1.1*max (ssBio[,s,])
    # Blank plot
    plot ( x=c(1,nT), y =c(1,yLim), type = "n", xlab = "", ylab = "Biomass (t)")
    # Create polygon vertices
    xPoly <- c(1:nT,nT:1)
    yPolySS <- c(ssBio[1,s,1:nT],ssBio[3,s,nT:1])
    yPolyMS <- c(msBio[1,s,1:nT],msBio[3,s,nT:1])
    # plot CI bands
    # SS
    polygon (x = xPoly, y = yPolySS, col = "steelblue", density=1)
    # MS
    polygon (x = xPoly, y = yPolyMS, col = "salmon", density=1)
    # Plot Medians
    lines ( x = 1:nT, y = ssBio[2,s,], col = "steelblue",lty=2,lwd=2 )
    lines ( x = 1:nT, y = msBio[2,s,], col = "salmon",lty=2,lwd=2 )
    # Plot true biomass
    lines ( x=1:nT, y = omBio[s,],col="black", lwd = 0.8)
  }
  return()
}

# plotRepBCF()
# Function that will open a supplied saved blob object and plot a 
# given replicate's true and estimated time series of biomass, 
# catch, fishing mortality and IoA data. 
# inputs:   rep=replicate number for plotting; 
#           sim=number indicating blob to load (alpha order in project folder)
#           folder=name of folder/blob file (supercedes sim number)
# output:   NULL
# usage:    post-simulation run, plotting performance
plotBCF <- function ( rep = 1, sim = 1, folder = NULL, est="MLE" )
{
  # First, if a folder name isn't nominated, use sim number
  # to find the folder name
  if (is.null(folder))
  {
    # List directories in project folder, remove "." from list
    dirList <- list.dirs (path="./project",full.names = FALSE,
                          recursive=FALSE)
    simList <- dirList[grep(pattern="sim",x=dirList)]
    folder <- simList[sim]
  }
  # Load the nominated blob
  blobFileName <- paste(folder,".RData",sep="")
  blobPath <- file.path(getwd(),"project",folder,blobFileName)
  load ( file = blobPath )

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctl$scenario
  mp        <- blob$ctl$mp
  stamp     <- paste(scenario,":",mp,sep="")
  repCount  <- paste("Replicate ",rep,"/",blob$ctl$nReps,sep="")

  # Recover blob elements for plotting
  nS <- blob$ctl$nS
  nT <- blob$ctl$nT

  # Species names
  specNames <- blob$ctl$specNames

  # True OM quantities
  omBt  <- blob$om$Bt[rep,,]
  Ct    <- blob$om$Ct[rep,,]
  It    <- blob$om$It[rep,,]
  Ft    <- blob$om$Ft[rep,,]

  
  if ( est == "MLE" )
  { # Single species model
    ssBt  <- blob$am$ss$Bt[rep,,]
    ssq   <- blob$am$ss$q[rep,]

    # Multispecies model
    msBt  <- blob$am$ms$Bt[rep,,]
    msq   <- blob$am$ms$q[rep,]  
  }
  if ( est == "MCMC" )
  {
    # single species
    mcBioSSMed <- apply(X=blob$am$ss$mcBio,FUN=quantile,MARGIN=c(1,2,4),na.rm=TRUE,probs=0.5)
    mcParSSMed <- apply(X=blob$am$ss$mcPar,FUN=quantile,MARGIN=c(1,2,4),na.rm=TRUE,probs=0.5)
    ssBt       <- mcBioSSMed[rep,,]
    ssq        <- mcParSSMed[rep,,"q"]

    # multispeces
    mcBioMSMed <- apply(X=blob$am$ms$mcBio,FUN=quantile,MARGIN=c(1,2,4),na.rm=TRUE,probs=0.5)
    mcParMSMed <- apply(X=blob$am$ms$mcPar,FUN=quantile,MARGIN=c(1,2,4),na.rm=TRUE,probs=0.5)
    msBt       <- mcBioMSMed[rep,,]
    msq        <- mcParMSMed[rep,,"q"]
  }
  

  # Set up plot window
  par ( mfrow = c(3,nS), mar = c(1,4,1,0), oma = c(3,0,1,0.5) )
  # Plot biomass, actual and estimated, including 2 index series,
  # scaled by model estimated q
  for ( s in 1:nS )
  {
    if ( s == 1 ) yLab <- "Biomass (t)" else yLab <- ""
    maxBt <- 1.05*max ( omBt[s,], ssBt[s,], msBt[s,],It[s,]/ssq[s],It[s,]/msq[s])
    plot    ( x = c(1,nT), y = c(0,maxBt), type = "n", 
              ylim = c(0,maxBt), ylab = yLab, las = 1, xlab = "" ,
              main = specNames[s] )
    lines   ( x = 1:nT, y = omBt[s,], col = "black", lwd = 2)
    lines   ( x = 1:nT, y = ssBt[s,], col = "steelblue", lwd = 2, lty = 2 )
    lines   ( x = 1:nT, y = msBt[s,], col = "salmon", lwd = 2, lty = 2 )
    points  ( x = 1:nT, y = It[s,]/ssq[s], pch = 2, cex = 0.8, col="grey50" )
    points  ( x = 1:nT, y = It[s,]/msq[s], pch = 5, cex = 0.8, col="grey50" )
  }
  # Now plot catch
  for ( s in 1:nS )
  {
    maxCt <- max ( Ct[s,])
    if ( s == 1 ) yLab <- "Catch (t)" else yLab <- ""
    plot  ( x = 1:nT, y = Ct[s,], col = "blue", lwd = 2, type = "l",
            ylim = c(0,maxCt), ylab = yLab, las = 1, xlab = "" )
  }
  # Now F
  for ( s in 1:nS )
  {
    if ( s == 1 ) yLab <- "Fishing Mortality" else yLab <- ""
    maxFt <- max ( Ft[s,])
    plot  ( x = 1:nT, y = Ft[s,], col = "black", lwd = 2, type = "l",
            ylim = c(0,maxFt), ylab = yLab, las = 1, xlab = "" )
  }
  mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5)
  mtext ( text = c(stamp,repCount),side=1, outer = TRUE, 
          at = c(0.9,0.1),padj=2,col="grey50",cex=0.8 )
}


# plotSimPerf()
# A function to plot simulation-estimation performance for a whole
# collection of replicates. Resulting plot is a dot and bar graph
# showing relative error distributions for nominated parameters
# inputs:   sim = numeric indicator of simulation in project folder
#           folder = optional character indicating sim folder name
#           pars = nominated estimated leading and derived parameters 
# outputs:  NULL
# usage:    
plotSimPerf <- function ( sim = 1, folder = NULL, 
                          pars = c("Bmsy","Fmsy","q","dep","BnT") )
{
  # First, if a folder name isn't nominated, use sim number
  # to find the folder name
  if (is.null(folder))
  {
    # List directories in project folder, remove "." from list
    dirList <- list.dirs (path="./project",full.names = FALSE,
                          recursive=FALSE)
    simList <- dirList[grep(pattern="sim",x=dirList)]
    folder <- simList[sim]
  }
  # Load the nominated blob
  blobFileName <- paste(folder,".RData",sep="")
  blobPath <- file.path(getwd(),"project",folder,blobFileName)
  load ( file = blobPath )

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctl$scenario
  mp        <- blob$ctl$mp
  stamp     <- paste(scenario,":",mp,sep="")

  # Recover blob elements for plotting
  nS        <- blob$ctl$nS
  specNames <- blob$ctl$specNames

  # Recover relative error distributions
  ssRE <- blob$am$ss$err.mle
  msRE <- blob$am$ms$err.mle

  # Create a wrapper function for generating quantiles
  quantWrap <- function ( entry = 1, x = ssRE, ... )
  {
    quants <- apply ( X = x[[entry]], FUN = quantile, ... )
    return(t(quants))
  }

  # generate quantiles
  ssQuant <- lapply ( X = seq_along(ssRE), FUN = quantWrap, x=ssRE, 
                      MARGIN = 2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  msQuant <- lapply ( X = seq_along(msRE), FUN = quantWrap, x=msRE, 
                      MARGIN = 2, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)


  # get names of parameters
  ssPars <- names ( ssRE )
  msPars <- names ( msRE )
  names ( ssQuant) <- ssPars
  names ( msQuant) <- msPars

  # Create an array of quantile values 
  # (dims = species,percentiles,parameters,models)
  quantiles <- array ( NA, dim = c(nS,3,length(pars),2))
  dimnames (quantiles) <- list (  1:nS,
                                  c("2.5%","50%","97.5%"), 
                                  pars,
                                  c("ss","ms") )
  # populate table
  for ( s in 1:nS )
  {
    for (par in 1:length(pars))
    {
      quantiles[s,,par,1] <- ssQuant[pars[par]][[1]][s,]
      quantiles[s,,par,2] <- msQuant[pars[par]][[1]][s,]
    }
  }
  
  # Set plotting window
  par (mfrow = c(1,3), mar = c(3,0,1,1), oma = c(3,0,2,0) )

  for ( s in 1:nS )
  {
    med <- quantiles[s,"50%",,]
    q975 <- quantiles[s,"97.5%",,]
    q025 <- quantiles[s,"2.5%",,]

    # Plot main dotchart
    plotTitle <- specNames[s]
    if ( s == 1) dotchart ( x = t(med), xlim = c(-1.5,1.5
      ),main=plotTitle,
                            pch = 16)
    else dotchart ( x = t(med), xlim = c(-1.5,1.5), main = plotTitle,
                    pch = 16)
    
    # Now add segments
    for ( par in length(pars):1)
    {
      parIdx <- length(pars) - par + 1
      plotY <- 2 * par + 2 * (par-1)
      segments( q025[pars[parIdx],"ss"],plotY-1,q975[pars[parIdx],"ss"],plotY-1,lty=1,col="grey60", lwd=2)  
      segments( q025[pars[parIdx],"ms"],plotY,q975[pars[parIdx],"ms"],plotY,lty=1,col="grey60", lwd=2)  
    }

    abline ( v = 0, lty = 3, lwd = 0.5 )

  } 
  title <- paste("Relative error distributions, nReps = ", length(blob$goodReps))
  mtext ( text = title, outer = TRUE, side = 3 )
  mtext ( text = "Relative Error", side = 1, outer = TRUE)
  mtext ( text = c(stamp),side=1, outer = TRUE, 
          at = c(0.75),padj=1,col="grey50" )

}
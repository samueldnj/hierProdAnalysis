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

# plotRep()
# Function that will open a supplied saved blob object and plot a 
# given replicate's true and estimated time series of biomass, 
# catch, fishing mortality and IoA data. 
# inputs:   rep=replicate number for plotting; 
#           sim=number indicating blob to load (alpha order in project folder)
#           folder=name of folder/blob file (supercedes sim number)
# output:   NULL
# usage:    post-simulation run, plotting performance
plotReps <- function ( rep = 1, sim = 1, folder = NULL )
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

  # Recover blob elements for plotting
  nS <- blob$ctl$nS
  nT <- blob$ctl$nT

  # True OM quantities
  omBt  <- blob$om$Bt[rep,,]
  Ct    <- blob$om$Ct[rep,,]
  It    <- blob$om$It[rep,,]
  Ft    <- blob$om$Ft[rep,,]

  # Single species model
  ssBt  <- blob$am$ss$Bt[rep,,]
  ssq   <- blob$am$ss$q[rep,]

  # Multispecies model
  msBt  <- blob$am$ms$Bt[rep,,]
  msq   <- blob$am$ms$q[rep,]

  # Set up plot window
  par ( mfrow = c(3,nS), mar = c(2,4,1,1), oma = c(3,1,1,1) )
  # Plot biomass, actual and estimated, including 2 index series,
  # scaled by model estimated q
  for ( s in 1:nS )
  {
    if ( s == 1 ) yLab <- "Biomass (t)" else yLab <- ""
    maxBt <- 1.05*max ( omBt[s,], ssBt[s,], msBt[s,],It[s,]/ssq[s],It[s,]/msq[s])
    plot    ( x = 1:nT, y = omBt[s,], col = "red", lwd = 2, type = "l", 
              ylim = c(0,maxBt), ylab = yLab, las = 1, xlab = "" ,
              main = paste("Species ", s, sep = "") )
    lines   ( x = 1:nT, y = ssBt[s,], col = "red", lwd = 2, lty = 2 )
    lines   ( x = 1:nT, y = msBt[s,], col = "red", lwd = 2, lty = 3 )
    points  ( x = 1:nT, y = It[s,]/ssq[s], pch = 2, cex = 1 )
    points  ( x = 1:nT, y = It[s,]/msq[s], pch = 5, cex = 1 )
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
  mtext ( text = "Year", outer = TRUE, side = 1)
}

plotSimPerf <- function ( sim = 1, folder = NULL, 
                          pars = c("Bmsy","Fmsy","q","dep","BnT"))
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

  # Recover blob elements for plotting
  nS <- blob$ctl$nS
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
    plotTitle <- paste ( "Species ", s, sep = "")
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

}
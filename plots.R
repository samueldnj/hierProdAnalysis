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

plotPriorSens <- function ( tableName = "qSens_MRE",
                            par = "q" )
{
  # plotPriorSens()
  # Reads in a stats table for prior sensitivity analyses
  # and plots the performance
  # inputs:     tableName = character name of stats table
  #             par=paremeter that is affected by sens analysis
  # outputs:    NULL
  # side-effs:  plots to quartz
}


plotObsRasters <- function  ( tableName = "obsErr_MARE",
                              pars = c("BnT","Umsy","q","Dep"),
                              wdh = 18, hgt = 18,
                              breakRange = c(-1,1),
                              compMean = FALSE,
                              ... 
                            )
{
  # plotObsRasters()
  # Reads in a statistics table produced by .statTableXXX()
  # and produces plots of performance contours.
  # inputs:     tableName=charactre vector of file name root
  #             pars=character vector of parameters to plot contours for
  #             wdh=numeric width (inches) of plot
  #             hgt=numeric height (inches) of plot
  #             breakRange=range of colour breaks for raster plots
  #             compMean=Switch for computing the mean (TRUE) or 
  #                       sum (FALSE) for the complex contour
  #             ...=args to pass to plotObsCompContour()
  # outputs:    NA
  # side eff:   Prints pdf plots of each MP and nS combination
  #             to getwd()/project/figs/<tableName>/
  
  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )

  tabFolder <- file.path(getwd(),"project","figs",tableName )
  dir.create(tabFolder)
  # Now get the list of MPs
  MPs <- unique( table$mp )
  # browser()
  # Loop over MPs, and plot one set of raster fields for
  # each of the other levels
  for (mp in MPs)
  {
    mpFolder  <- file.path(tabFolder,mp)
    dir.create(mpFolder)
    mpTab     <- table %>% filter( mp == mp )
    nS        <- unique(mpTab$nS)
    levels    <- list (nS = nS )
    grid      <- expand.grid(levels)
    # browser()
    for( g in 1:nrow(grid) )
    {
      # Set up output file
      outFile <- paste( grid$nS[g], 
                        mp,
                        "obsErrComp.pdf", sep = "")
      outFile <- file.path(mpFolder,outFile)
      pdf ( file = outFile, width = wdh, height = hgt )
      plotObsCompContour( tableName = fileName,
                          mpLabel = mp,
                          pars = pars,
                          nSp = grid[g,"nS"],
                          tcex = 2,
                          breakRange = breakRange,
                          compMean = compMean,
                          ... )
      dev.off()
      outFile <- paste( grid$nS[g], 
                        mp,
                        "obsErr_MS.pdf", sep = "")
      outFile <- file.path(mpFolder,outFile)
      pdf ( file = outFile, width = wdh, height = hgt )
      plotObsModContour(  tableName = fileName,
                          mpLabel = mp,
                          pars = pars,
                          nSp = grid[g,"nS"],
                          model = "ms", tcex = 2,
                          breakRange = breakRange,
                          compMean = compMean )
      dev.off()
      outFile <- paste( grid$nS[g], 
                        "obsErr_SS.pdf", sep = "")
      outFile <- file.path( mpFolder, outFile )
      pdf ( file = outFile, width = wdh, height = hgt )
      plotObsModContour(  tableName = fileName,
                          mpLabel = mp,
                          pars = pars,
                          nSp = grid[g,"nS"],
                          model = "ss", tcex =2,
                          breakRange = breakRange,
                          compMean = compMean )
      dev.off()
    }
  }
  # No return
}


plotObsModContour <- function (   tableName   = "obsErr_MARE.csv", 
                                  pars        = c("BnT", "Umsy", "q", "Dep"),
                                  nSp         = 2,
                                  method      = "deviance",
                                  compMean    = FALSE,                                
                                  tcex        = 1,
                                  breakRange  = c(-1,1),
                                  mpLabel     = "corrTRUE_SigPriorIG",
                                  model       = "ss" 
                              )
{
  # plotObsModContour()
  # Plots the observation model experiment MSE contours for
  # estimation performance of given pars by species
  # and for the complex.
  # inputs:     tableName=charactre vector of file name root
  #             pars=character vector of parameters to plot contours for
  #             nSp=number of species (scenario selector)
  #             method=method for comparing ms and ns models,
  #                     deviance => log2(ss/ms),
  #                     difference => abs(ss) - abs(ms)
  #             compMean=compute mean for complex (TRUE), or sum (FALSE)?
  #             tcex=cex for text calls (prints numbers over raster cells)
  #             breakRange=range of colour breaks for raster plots
  #             mpLabel=character vector selecting MP
  #             model=single species "ss" or multi-species "ms"
  # outputs:    NA
  # side eff:   Prints plots of model estimation performance contours for
  #             selected pars and MP

  # Load stat table
  tablePath <- file.path ( getwd(), "project/Statistics", tableName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )
  # reduce to the correct MP
  table <-  table  %>% 
            filter( mp == mpLabel, nS == nSp ) %>%
            mutate( CV = round( sqrt( exp( tau2OM ) -1 ), digits = 2 ) )

  # browser()
  
  species   <- as.character(unique( table$species))
  nS <- length(species)
  specList  <- vector( mode="list", length=nS)

  # Axis ticks:
  CVs   <- unique( table[ , "CV" ] )
  CVs   <- CVs[ order( CVs ) ]
  # scenario name values
  tau2  <- unique( table[ , "tau2OM" ] ) 
  tau2  <- tau2[ order( tau2 ) ]
 
  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    z <- array( NA, dim=c(length(CVs),length(CVs),length(pars)+1),
                dimnames=list(CVs[length(CVs):1],CVs,c(pars,"HessPD")))

    for (xIdx in 1:length(CVs))
    {
      for (yIdx in 1:length(CVs))
      {
        # Get x and y values
        xVal <- tau2[xIdx]
        yVal <- tau2[yIdx]
        # Get rows
        scenName <- paste( "nS", nS, "_CVc(", xVal, ",", yVal, ")", sep = "" )
        # Get rows
        rows <- which ( specTab[,"scenario"] == scenName )
        # browser()
        if( length(rows) == 0 ) next
        colNames <- paste(model,c(pars,"HessPD"),sep="")
        z[length(CVs)-yIdx+1,xIdx,] <- as.numeric(specTab[rows,colNames])
      }
    } 
    # browser()
    specList[[s]] <- brick( z, xmn = min(CVs), xmx = max(CVs), ymn = min(CVs), ymx=max(CVs) )
  }

  names(specList) <- species

  par(mfrow = c((nS+1),length(pars)+1),  mar = c(1.5,1.5,1.5,1.5), oma = c(5,5,5,5), las=1)
  for (s in 1:nS)
  {
    # Create a raster from the stat table info
    plotBrick <- specList[[s]]
    # Make titles for the plots
    specName <- species[s]

    titles <- names(plotBrick)
    l <- length(titles)

    colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
    colBreaksHess <- seq(0,200,length = 8)
    colScale  <- brewer.pal(8, "Reds")
    colScale      <- c(colScale[8:1],colScale)
    
    # plot the rasters, if first species plot main titles
    if ( s == 1 )
    {
      for(k in 1:(l-1))
      {
        plot( plotBrick[[ k ]],main = titles[ k ], legend = FALSE, 
              asp=NA, ylab=specName, col = colScale, breaks = colBreaks )
        text( plotBrick[[ k ]], digits = 2, halo = T, cex = tcex ) 
      }
      plot( plotBrick[[ l ]],main = titles[ l ], legend = FALSE, 
              asp=NA, ylab=specName, col = colScale[9:16], breaks = colBreaksHess )
        text( plotBrick[[ l ]], digits = 2, halo = T, cex = tcex ) 
    }
    else {
      # browser()
      for(k in 1:( l - 1 ) )
      {
        plot( plotBrick[[k]],main = "", legend = FALSE, 
              asp=NA, ylab=specName, col = colScale, breaks = colBreaks )
        text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )  
      }
      plot( plotBrick[[l]],main = "", legend = FALSE, 
              asp=NA, ylab=specName, col = colScale[9:16], breaks = colBreaksHess )
      text( plotBrick[[l]], digits = 2, halo = T, cex = tcex ) 
    }   
  }

  # browser()
  # Now create an average/sum raster brick for the complex
  for(k in 1:l)
  {
    # use leftover plotBrick in memory to compute complex
    plotBrick[[k]] <- 0
    for( s in 1:length(species) )
    {
      plotBrick[[k]] <- plotBrick[[k]] + specList[[s]][[k]]
    }
    if( compMean ) plotBrick[[k]] <- plotBrick[[k]]/length(species)
  }

  titles <- names(plotBrick)
  
  # browser()
  
  for(k in 1:( l - 1 ) )
  {
    plot( plotBrick[[k]],main = "", legend = FALSE, 
          asp=NA, ylab="Complex", col = colScale, breaks = colBreaks )
    text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )   
  }
  plot( plotBrick[[l]],main = "", legend = FALSE, 
        asp=NA, ylab="Complex", col = colScale[9:16], breaks = nSp*colBreaksHess )
  text( plotBrick[[l]], digits = 2, halo = T, cex = tcex )    

  grid.text(  "Dover Sole Survey CV", 
              x=unit(0.5, "npc"), y=unit(0.02, "npc"), rot=0)
  grid.text(  "English Sole Survey CV", 
              x=unit(0.02, "npc"), y=unit(0.5, "npc"), rot=90)
  grid.text(  paste("Raw error ", model, " model.", sep = ""), 
              x=unit(0.5, "npc"), y=unit(0.98, "npc"), rot=0)
  grid.text(  tableName,
              x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
  grid.text(  mpLabel,
              x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)
  # grid.text(  paste( "lastNegCorr = ", negCorr, sep = "" ),
              # x=unit(0.15, "npc"), y=unit(0.02, "npc"), rot=0)
  # grid.text(  paste( "tUtrough = ", tUtr, sep = "" ),
  #             x=unit(0.15, "npc"), y=unit(0.05, "npc"), rot=0)
}


# plotObsCompContour()
# Plots the correlation comparative experiment performance contours for
# estimation of BnT and Umsy (abundance and productivity) by species
# and for the complex. The plots are panel plots of Rasters, so need
# a little love as far as sizing goes.
# inputs: table=character name of table generated by .statTableXXX
#         mpLabel=character name of MP label
#         axes=character 2-vector for axes of raster plots
#         nSp=integer number of species in scenario
#         method="deviance" for log2(SS/MS), difference for SS-MS
#         compMean=compute the mean for the complex?
#         tUtr=integer value for tUtrough in Fhist scenarios
#         tcex=cex of text in raster cells
plotObsCompContour <- function (  tableName   = "obsErr_MARE.csv", 
                                  axes        = c("corr", "kappaMult"),
                                  pars        = c("BnT", "Umsy", "q", "Dep"),
                                  nSp         = 2,
                                  method      = "deviance",
                                  compMean    = FALSE,                                
                                  tcex        = 2,
                                  breakRange  = c(-2,2),
                                  mpLabel     = "corrTRUE_SigPriorIG" 
                                  )
{
  # Load stat table
  tablePath <- file.path ( getwd(),"project/Statistics",tableName)
  table <- read.csv (tablePath, header=TRUE, stringsAsFactors=FALSE)
  # reduce to the correct MP
  table <- table %>% filter(mp == mpLabel )

  # Filter by other 
  table <- table %>% 
           filter( nS == nSp ) %>%
           mutate ( CV = round( sqrt( exp(tau2OM) - 1 ), digits = 2 ) )


  # browser()
  # Now make a pallette for the colours
  colScaleHi    <- brewer.pal( 8, "Greens" )
  colScaleLo    <- brewer.pal( 8, "Reds" )
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colScale      <- c(colScaleLo[8:1],colScaleHi)

  # calculate response
  # Deviance is log2(ss/ms)
  if (method == "deviance" )
  {
    table <- table %>% mutate(  BnT     = log2(abs(ssBnT/msBnT)),
                                Umsy    = log2(abs(ssUmsy/msUmsy)),
                                q       = log2(abs(ssq/msq)),
                                Dep     = log2(abs(ssDep/msDep)) )
  }   
  # Difference is just ss - ms
  if ( method == "difference")
  {
    table <- table %>% mutate(  BnT     = abs(ssBnT) - abs(msBnT),
                                Umsy    = abs(ssUmsy) - abs(msUmsy),
                                q       = abs(ssq) - abs(msq),
                                Dep     = abs(ssDep) - abs(msDep) )
  }

  species   <- as.character(unique( table$species))
  nS <- length(species)
  specList  <- vector( mode="list", length=length(species))

  CVs   <- unique( table[ , "CV" ] )
  CVs   <- CVs[ order( CVs ) ]

  tau2  <- unique( table[ , "tau2OM" ] ) 
  tau2  <- tau2[ order( tau2 ) ]
  # browser()
  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    z <- array( NA, dim = c( length( CVs ), length( CVs ), length(pars) ),
                dimnames=list(  CVs[ length( CVs ):1 ], CVs,
                                pars ) )

    for (xIdx in 1:length(CVs))
    {
      for (yIdx in 1:length(CVs))
      {
        # Get x and y values
        xVal <- tau2[xIdx]
        yVal <- tau2[yIdx]
        scenName <- paste( "nS", nS, "_CVc(", xVal, ",", yVal, ")", sep = "" )
        # Get rows
        rows <- which ( specTab[,"scenario"] == scenName )
        # browser()
        if (length(rows) == 0 ) next
        # Recover MSE comparison
        z[ length( CVs ) - yIdx + 1, xIdx, ] <- as.numeric(specTab[rows,pars])
      }
    } 
    specList[[s]] <- brick( z, xmn = min(CVs), xmx = max(CVs), ymn = min(CVs), ymx=max(CVs) )
  }
  names(specList) <- species

  # Start plotting

  par(mfrow = c(nS+1,length(pars)), mar = c(1.5,1.5,1.5,1.5), oma = c(5,5,5,5), las=1)
  for (s in 1:nS)
  {
    # Create a raster from the stat table info 
    plotBrick <- specList[[s]]
    # Make titles for the plots
    titles <- names(plotBrick)
    # browser()
    # plot the rasters, without legends
    for ( k in 1:length(titles) )
    {
      if(s == 1) plot( plotBrick[[k]],main = titles[k], legend = FALSE, 
                        breaks = colBreaks, col=colScale, asp=NA )
      else plot( plotBrick[[k]],main = "", legend = FALSE, 
                  breaks = colBreaks, col=colScale, asp=NA )
      if( k == 1 ) mtext ( text = species[s], side = 2 )
      text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )  
    }
  }

  # browser()
  # Now create the complex plot brick - idea: group by scenario
  # and MP, there should be a unique set of species for each, then
  # apply the complex rule
  compTable <- table %>%  group_by( scenario, mp ) %>%
                          summarise(  ssBnT = sum( ssBnT ),
                                      msBnT = sum( msBnT ),
                                      ssUmsy = sum( ssUmsy ),
                                      msUmsy = sum( msUmsy ),
                                      ssq = sum( ssq ),
                                      msq = sum( msq ),
                                      ssHessPD = sum( ssHessPD ),
                                      msHessPD = sum( msHessPD ),
                                      # ssMSY = sum( ssMSY ),
                                      # msMSY = sum( msMSY ),
                                      ssDep = sum( ssDep ),
                                      msDep = sum( msDep ),
                                      corr = mean(corr),
                                      kappaMult = mean(kappaMult),
                                      tUpeak = mean(tUpeak),
                                      tUtrough = mean(tUtrough),
                                      Umax = mean(Umax),
                                      tau2OM = mean (tau2OM),
                                      CV = mean ( CV ) )

  # Now make the comparison based on the compMean logical flag
  if( method == "deviance" )
  {
    compTable <-  compTable %>% 
                  mutate( BnT = log2( abs(ssBnT) / abs(msBnT) ),
                          Umsy = log2( abs(ssUmsy) / abs(msUmsy) ),
                          q = log2( abs(ssq) / abs(msq) ),
                          # MSY = log2( ssMSY / msMSY ),
                          Dep = log2 ( abs(ssDep) / abs(msDep) ) 
                        )
  }

  if( method == "difference" )
  {
    compTable <-  compTable %>% 
                  mutate( BnT = abs(ssBnT) - abs(msBnT) ,
                          Umsy = abs(ssUmsy) - abs(msUmsy) ,
                          q = abs(ssq) - abs(msq) ,
                          # MSY = abs(ssMSY) - abs(msMSY) ,
                          Dep = abs(ssDep) - abs(msDep) 
                        )
    if( compMean )
      compTable <-  compTable %>% 
                    mutate( BnT = BnT/nS,
                            Umsy = Umsy/nS,
                            q = q/nS,
                            # MSY = MSY/nS,
                            Dep = Dep/nS )
  }

  # Make a raster brick with the selected pars
  z <- array( NA, dim = c( length( CVs ), length( CVs ), length(pars) ),
              dimnames=list(  CVs[ length(CVs):1], CVs,
                              pars ) )

  # browser()
  for (xIdx in 1:length(CVs))
  {
    for (yIdx in 1:length(CVs))
    {
      # Get x and y values
      xVal <- tau2[xIdx]
      yVal <- tau2[yIdx]
      # Get rows
      scenName <- paste( "nS", nS, "_CVc(", xVal, ",", yVal, ")", sep = "" )
      # Get rows
      rows <- which ( specTab[,"scenario"] == scenName )
      if( length(rows) == 0 ) next
      # Recover MSE comparison
      z[length(CVs) - yIdx + 1,xIdx,] <- as.numeric(compTable[rows,pars])
    }
  } 
  compBrick <- brick( z, xmn = min(CVs), xmx = max(CVs), ymn = min(CVs), ymx=max(CVs) )

  for(k in 1:length(titles))
  {
    plot( compBrick[[k]],main = "", legend = FALSE, 
          asp=NA, col = colScale, breaks = colBreaks )
      if( k == 1 ) mtext ( text = "Complex", side = 2 )
    text( compBrick[[k]], digits = 2, halo = T, cex = tcex )  
  }

  # Plot the legends
  plot( plotBrick[[1]], legend.only=TRUE, col=colScale, fill=colScale,
        legend.width=1, legend.shrink=0.75,
        axis.args = list( at = colBreaks,
                          labels = round(colBreaks,2),
                          cex.axis = 0.6),
        zlim = c(-3,3),
        legend.args = list( text="log2(ssMSE/msMSE)", 
                            side=4, font=2, 
                            line=3, cex=0.75, las=0) )

  # Labels
  grid.text(  "Dover Sole Survey CV", 
              x=unit(0.5, "npc"), y=unit(0.05, "npc"), rot=0)
  grid.text(  "English Sole Survey CV", 
              x=unit(0.02, "npc"), y=unit(0.5, "npc"), rot=90)
  grid.text(  "Relative errors ss/ms", 
              x=unit(0.5, "npc"), y=unit(0.95, "npc"), rot=0)
  grid.text(  tableName,
              x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
  grid.text(  mpLabel,
              x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)
  # grid.text(  paste( "lastNegCorr = ", negCorr, sep = "" ),
  #             x=unit(0.15, "npc"), y=unit(0.02, "npc"), rot=0)
  # grid.text(  paste( "tUtrough = ", tUtr, sep = "" ),
  #             x=unit(0.15, "npc"), y=unit(0.05, "npc"), rot=0)

}

# plotRasters()
# Reads in a statistics table produced by .statTableXXX()
# and produces plots of performance contours.
# inputs:     tableName=charactre vector of file name root
#             axes=
plotRasters <- function ( tableName = "RE_coarse_MARE",
                          axes = c("corr","kappaMult"),
                          pars = c("BnT","Umsy","q","Dep"),
                          wdh = 14, hgt = 18,
                          breakRange = c(-1,1),
                          compMean = FALSE,
                           ... )
{
  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )

  # browser()
  tabFolder <- file.path(getwd(),"project","figs",tableName )
  dir.create(tabFolder)
  # Now get the list of MPs
  MPs <- unique( table$mp )
  # browser()
  # Loop over MPs, and plot one set of raster fields for
  # each of the other levels
  for (mp in MPs)
  {
    mpFolder  <- file.path(tabFolder,mp)
    dir.create(mpFolder)
    mpTab     <- table %>% filter( mp == mp )
    nCorr     <- unique(mpTab$lastNegCorr)
    tUtrough  <- unique(mpTab$tUtrough)
    nS        <- unique(mpTab$nS)
    levels    <- list (nCorr = nCorr, tUtrough = tUtrough, nS = nS )
    grid      <- expand.grid(levels)
    for( g in 1:nrow(grid) )
    {
      # Set up output file
      outFile <- paste( grid[g,]$nS, 
                        "S_tUtr", grid[g,]$tUtrough,
                        "_nCorr", grid[g,]$nCorr,
                        ".pdf", sep = "")
      outFile <- file.path(mpFolder,outFile)
      pdf ( file = outFile, width = wdh, height = hgt )
      plotCompContour(  tableName = fileName,
                        mpLabel = mp,
                        axes = axes,
                        pars = pars,
                        nSp = grid[g,"nS"],
                        tUtr = grid[g,"tUtrough"],
                        negCorr = grid[g,"nCorr"],
                        tcex = 2,
                        breakRange = breakRange,
                        compMean = compMean,
                        ...)
      dev.off()
      outFile <- paste( grid[g,]$nS, 
                        "S_tUtr", grid[g,]$tUtrough,
                        "_nCorr", grid[g,]$nCorr,
                        "_MS.pdf", sep = "")
      outFile <- file.path(mpFolder,outFile)
      pdf ( file = outFile, width = wdh, height = hgt )
      plotModContour(  tableName = fileName,
                        mpLabel = mp,
                        axes = axes,
                        pars = pars,
                        nSp = grid[g,"nS"],
                        tUtr = grid[g,"tUtrough"],
                        negCorr = grid[g,"nCorr"],
                        model = "ms", tcex = 2,
                        breakRange = breakRange,
                        compMean = compMean )
      dev.off()
      outFile <- paste( grid[g,]$nS, 
                        "S_tUtr", grid[g,]$tUtrough,
                        "_nCorr", grid[g,]$nCorr,
                        "_SS.pdf", sep = "")
      outFile <- file.path( mpFolder, outFile )
      pdf ( file = outFile, width = wdh, height = hgt )
      plotModContour(  tableName = fileName,
                        mpLabel = mp,
                        axes = axes,
                        pars = pars,
                        nSp = grid[g,"nS"],
                        tUtr = grid[g,"tUtrough"],
                        negCorr = grid[g,"nCorr"],
                        model = "ss", tcex = 2,
                        breakRange = breakRange,
                        compMean = compMean )
      dev.off()
    }
  }
  # No return
}


# plotModContour()
# Plots the correlation  experiment MSE contours for
# estimation of BnT and Umsy (abundance and productivity) by species
# and for the complex.
plotModContour <- function (  tableName = "2sRE.csv", 
                              mpLabel   = "bothEff",
                              axes      = c("corr","kappaMult"),
                              pars      = c("BnT","Umsy","q"),
                              nSp       = 2,
                              model     = "ss",
                              compMean  = FALSE,
                              tUtr      = 15,
                              negCorr   = FALSE,
                              tcex      = 2,
                              breakRange= c(-1,1) )
{
  # Load stat table
  tablePath <- file.path ( getwd(), "project/Statistics", tableName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )
  # reduce to the correct MP
  table <- table  %>% filter( mp == mpLabel )
  # Filter by other 
  if( "nS" %in% names(table) ) 
    table <- table %>% 
             filter( nS == nSp )
  if( "tUtrough" %in% names(table) ) 
    table <- table %>% 
             filter( tUtrough == tUtr )
  if( "lastNegCorr" %in% names(table) ) 
    table <- table %>% 
             filter( lastNegCorr == negCorr )

  species   <- as.character(unique( table$species))
  nS <- length(species)
  specList  <- vector( mode="list", length=nS)

  if ( length(axes) != 2) return ("Error! length(axes) != 2.\n")
  # Pick out columns for axes
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }

  x   <- unique(table[,axesCols[1]])
  x   <- x[order(x)]
  y   <- unique(table[,axesCols[2]])
  y   <- y[order(y)]

  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    z <- array( NA, dim=c(length(y),length(x),length(pars)+1),
                dimnames=list(y[length(y):1],x,c(pars,"HessPD")))

    for (xIdx in 1:length(x))
    {
      for (yIdx in 1:length(y))
      {
        # Get x and y values
        xVal <- x[xIdx]
        yVal <- y[yIdx]
        # Get rows
        xRows <- which ( specTab[,axesCols[1]] == xVal)
        yRows <- which ( specTab[,axesCols[2]] == yVal)
        rows <- intersect(xRows,yRows)
        if( length(rows) == 0 ) next
        colNames <- paste(model,c(pars,"HessPD"),sep="")
        z[length(y)-yIdx+1,xIdx,] <- as.numeric(specTab[rows,colNames])
      }
    } 
    # browser()
    specList[[s]] <- brick( z, xmn = min(x), xmx = max(x), ymn = min(y), ymx=max(y) )
  }

  names(specList) <- species

  par(mfrow = c((nS+1),length(pars)+1),  mar = c(1.5,1.5,1.5,1.5), oma = c(5,5,5,5), las=1)
  for (s in 1:nS)
  {
    # Create a raster from the stat table info
    plotBrick <- specList[[s]]
    # Make titles for the plots
    specName <- species[s]

    titles <- names(plotBrick)
    l <- length(titles)

    colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
    colBreaksHess <- seq(0,200,length = 8)
    colScale  <- brewer.pal(8, "Reds")
    colScale      <- c(colScale[8:1],colScale)

    # plot the rasters, if first species plot main titles
    if ( s == 1 )
    {
      for(k in 1:(l-1))
      {
        plot( plotBrick[[k]],main = titles[k], legend = FALSE, 
              asp=NA, ylab=specName, col = colScale, breaks = colBreaks )
        text( plotBrick[[k]], digits = 2, halo = T, cex = tcex ) 
      }
      plot( plotBrick[[l]],main = titles[l], legend = FALSE, 
            asp=NA, ylab=specName, col = colScale[9:16], breaks = colBreaksHess )
      text( plotBrick[[l]], digits = 2, halo = T, cex = tcex ) 
    }
    else {
      # browser()
      for(k in 1:(l-1))
      {
        plot( plotBrick[[k]],main = "", legend = FALSE, 
              asp=NA, ylab=specName, col = colScale, breaks = colBreaks )
        text( plotBrick[[k]], digits = 2, halo = T, cex = tcex ) 
      }
      plot( plotBrick[[l]],main = "", legend = FALSE, 
            asp=NA, ylab=specName, col = colScale[9:16], breaks = colBreaksHess )
      text( plotBrick[[l]], digits = 2, halo = T, cex = tcex ) 
    }   
  }

  # browser()
  # Now create an average/sum raster brick for the complex
  for(k in 1:l)
  {
    # use leftover plotBrick in memory to compute complex
    plotBrick[[k]] <- 0
    for( s in 1:length(species) )
    {
      plotBrick[[k]] <- plotBrick[[k]] + specList[[s]][[k]]
    }
    if( compMean ) plotBrick[[k]] <- plotBrick[[k]]/length(species)
  }

  titles <- names(plotBrick)

  # Get colour breaks, based on range of data
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colBreaksHess <- seq(0,200,length = 8)
  colScale  <- brewer.pal(8, "Reds")
  colScale      <- c(colScale[8:1],colScale)

  for( k in 1:( l - 1 ) )
  {
    plot( plotBrick[[k]],main = "", legend = FALSE, 
          asp=NA, ylab="Complex", col = colScale, breaks = colBreaks )
    text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )  
  }
  plot( plotBrick[[l]],main = "", legend = FALSE, 
        asp=NA, ylab="Complex", col = colScale, breaks = nSp * colBreaksHess )
  text( plotBrick[[l]], digits = 2, halo = T, cex = tcex )  

  grid.text(  axes[1], 
              x=unit(0.5, "npc"), y=unit(0.02, "npc"), rot=0)
  grid.text(  axes[2], 
              x=unit(0.02, "npc"), y=unit(0.5, "npc"), rot=90)
  grid.text(  paste("Raw error ", model, " model.", sep = ""), 
              x=unit(0.5, "npc"), y=unit(0.98, "npc"), rot=0)
  grid.text(  tableName,
              x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
  grid.text(  mpLabel,
              x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)
  grid.text(  paste( "lastNegCorr = ", negCorr, sep = "" ),
              x=unit(0.15, "npc"), y=unit(0.02, "npc"), rot=0)
  grid.text(  paste( "tUtrough = ", tUtr, sep = "" ),
              x=unit(0.15, "npc"), y=unit(0.05, "npc"), rot=0)
}



# plotCorrCompContour()
# Plots the correlation comparative experiment performance contours for
# estimation of BnT and Umsy (abundance and productivity) by species
# and for the complex. The plots are panel plots of Rasters, so need
# a little love as far as sizing goes.
# inputs: table=character name of table generated by .statTableXXX
#         mpLable=character name of MP label
#         axes=character 2-vector for axes of raster plots
#         nSp=integer number of species in scenario
#         method="deviance" for log2(SS/MS), difference for SS-MS
#         compMean=compute the mean for the complex?
#         tUtr=integer value for tUtrough in Fhist scenarios
#         negCorr=logical value of simCtl par lastNegCorr in corr scenarios
#         tcex=cex of text in taster cells
plotCompContour <- function ( tableName   = "2sRE.csv", 
                              mpLabel     = "corrTRUE_SigPriorIG",
                              axes        = c("corr","kappaMult"),
                              pars        = c("BnT","Umsy","q"),
                              nSp         = 2,
                              method      = "deviance",
                              compMean    = FALSE,
                              tUtr        = 15,
                              negCorr     = FALSE,
                              tcex        = 2,
                              breakRange  = c(-2,2) 
                              )
{
  # browser()
  # Load stat table
  tablePath <- file.path ( getwd(),"project/Statistics",tableName)
  table <- read.csv (tablePath, header=TRUE, stringsAsFactors=FALSE)
  # reduce to the correct MP
  table <- table %>% filter(mp == mpLabel )

  # Filter by other 
  if( "nS" %in% names(table) ) 
    table <- table %>% 
             filter( nS == nSp )
  if( "tUtrough" %in% names(table) ) 
    table <- table %>% 
             filter( tUtrough == tUtr )
  if( "lastNegCorr" %in% names(table) ) 
    table <- table %>% 
             filter( lastNegCorr == negCorr )

  # browser()
  # Now make a pallette for the colours
  colScaleHi    <- brewer.pal( 8, "Greens" )
  colScaleLo    <- brewer.pal( 8, "Reds" )
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colScale      <- c(colScaleLo[8:1],colScaleHi)

  # calculate response
  # Deviance is log2(ss/ms)
  if (method == "deviance" )
  {
    table <- table %>% mutate(  BnT     = log2(abs(ssBnT/msBnT)),
                                Umsy    = log2(abs(ssUmsy/msUmsy)),
                                q       = log2(abs(ssq/msq)),
                                Dep     = log2(abs(ssDep/msDep)) )
  }   
  # Difference is just ss - ms
  if ( method == "difference")
  {
    table <- table %>% mutate(  BnT     = abs(ssBnT) - abs(msBnT),
                                Umsy    = abs(ssUmsy) - abs(msUmsy),
                                q       = abs(ssq) - abs(msq),
                                Dep     = abs(ssDep) - abs(msDep) )
  }

  species   <- as.character(unique( table$species))
  nS <- length(species)
  specList  <- vector( mode="list", length=length(species))

  # Pick out columns for axes
  # browser()
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }

  x   <- unique(table[,axesCols[1]])
  x   <- x[order(x)]
  y   <- unique(table[,axesCols[2]])
  y   <- y[order(y)]
  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    z <- array( NA, dim = c( length( y ), length( x ), length(pars) ),
                dimnames=list(  y[ length( y ):1 ], x,
                                pars ) )

    for (xIdx in 1:length(x))
    {
      for (yIdx in 1:length(y))
      {
        # Get x and y values
        xVal <- x[xIdx]
        yVal <- y[yIdx]
        # Get rows
        xRows <- which ( specTab[,axesCols[1]] == xVal)
        yRows <- which ( specTab[,axesCols[2]] == yVal)
        rows <- intersect(xRows,yRows)
        # Recover MSE comparison
        # browser()
        z[ length( y ) - yIdx + 1, xIdx, ] <- as.numeric(specTab[rows,pars])
      }
    } 
    specList[[s]] <- brick( z, xmn = min(x), xmx = max(x), ymn = min(y), ymx=max(y) )
  }
  names(specList) <- species

  # Start plotting

  par(mfrow = c(nS+1,length(pars)), mar = c(1.5,1.5,1.5,1.5), oma = c(5,5,5,5), las=1)
  for (s in 1:nS)
  {
    # Create a raster from the stat table info 
    plotBrick <- specList[[s]]
    # Make titles for the plots
    titles <- names(plotBrick)
    # browser()
    # plot the rasters, without legends
    for ( k in 1:length(titles) )
    {
      if(s == 1) plot( plotBrick[[k]],main = titles[k], legend = FALSE, 
                        breaks = colBreaks, col=colScale, asp=NA )
      else plot( plotBrick[[k]],main = "", legend = FALSE, 
                  breaks = colBreaks, col=colScale, asp=NA )
      if( k == 1 ) mtext ( text = species[s], side = 2 )
      text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )  
    }
  }

  # browser()
  # Now create the complex plot brick - idea: group by scenario
  # and MP, there should be a unique set of species for each, then
  # apply the complex rule
  compTable <- table %>%  group_by( scenario, mp ) %>%
                          summarise(  ssBnT = sum( ssBnT ),
                                      msBnT = sum( msBnT ),
                                      ssUmsy = sum( ssUmsy ),
                                      msUmsy = sum( msUmsy ),
                                      ssq = sum( ssq ),
                                      msq = sum( msq ),
                                      ssHessPD = sum( ssHessPD ),
                                      msHessPD = sum( msHessPD ),
                                      # ssMSY = sum( ssMSY ),
                                      # msMSY = sum( msMSY ),
                                      ssDep = sum( ssDep ),
                                      msDep = sum( msDep ),
                                      corr = mean(corr),
                                      kappaMult = mean(kappaMult),
                                      tUpeak = mean(tUpeak),
                                      tUtrough = mean(tUtrough),
                                      Umax = mean(Umax) )

  # Now make the comparison based on the compMean logical flag
  if( method == "deviance" )
  {
    compTable <-  compTable %>% 
                  mutate( BnT = log2( abs(ssBnT) / abs(msBnT) ),
                          Umsy = log2( abs(ssUmsy) / abs(msUmsy) ),
                          q = log2( abs(ssq) / abs(msq) ),
                          # MSY = log2( ssMSY / msMSY ),
                          Dep = log2 ( abs(ssDep) / abs(msDep) ) 
                        )
  }

  if( method == "difference" )
  {
    compTable <-  compTable %>% 
                  mutate( BnT = abs(ssBnT) - abs(msBnT) ,
                          Umsy = abs(ssUmsy) - abs(msUmsy) ,
                          q = abs(ssq) - abs(msq) ,
                          # MSY = abs(ssMSY) - abs(msMSY) ,
                          Dep = abs(ssDep) - abs(msDep) 
                        )
    if( compMean )
      compTable <-  compTable %>% 
                    mutate( BnT = BnT/nS,
                            Umsy = Umsy/nS,
                            q = q/nS,
                            # MSY = MSY/nS,
                            Dep = Dep/nS )
  }

  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(compTable) == axes[i])
  }

  # Make a raster brick with the selected pars
  z <- array( NA, dim = c( length( y ), length( x ), length(pars) ),
              dimnames=list(  y[ length(y):1], x,
                              pars ) )
  # browser()
  for (xIdx in 1:length(x))
  {
    for (yIdx in 1:length(y))
    {
      # Get x and y values
      xVal <- x[xIdx]
      yVal <- y[yIdx]
      # Get rows
      xRows <- which ( compTable[,axesCols[1]] == xVal)
      yRows <- which ( compTable[,axesCols[2]] == yVal)
      rows <- intersect(xRows,yRows)
      # Recover MSE comparison
      z[length(y) - yIdx + 1,xIdx,] <- as.numeric(compTable[rows,pars])
    }
  } 
  compBrick <- brick( z, xmn = min(x), xmx = max(x), ymn = min(y), ymx=max(y) )

  for(k in 1:length(titles))
  {
    plot( compBrick[[k]],main = "", legend = FALSE, 
          asp=NA, col = colScale, breaks = colBreaks )
      if( k == 1 ) mtext ( text = "Complex", side = 2 )
    text( compBrick[[k]], digits = 2, halo = T, cex = tcex )  
  }

  # Plot the legends
  plot( plotBrick[[1]], legend.only=TRUE, col=colScale, fill=colScale,
        legend.width=1, legend.shrink=0.75,
        axis.args = list( at = colBreaks,
                          labels = round(colBreaks,2),
                          cex.axis = 0.6),
        zlim = c(-3,3),
        legend.args = list( text="log2(ssMSE/msMSE)", 
                            side=4, font=2, 
                            line=3, cex=0.75, las=0) )

  # Labels
  grid.text(  axes[1], 
              x=unit(0.5, "npc"), y=unit(0.05, "npc"), rot=0)
  grid.text(  axes[2], 
              x=unit(0.02, "npc"), y=unit(0.5, "npc"), rot=90)
  grid.text(  "Relative errors ss/ms", 
              x=unit(0.5, "npc"), y=unit(0.95, "npc"), rot=0)
  grid.text(  tableName,
              x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
  grid.text(  mpLabel,
              x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)
  grid.text(  paste( "lastNegCorr = ", negCorr, sep = "" ),
              x=unit(0.15, "npc"), y=unit(0.02, "npc"), rot=0)
  grid.text(  paste( "tUtrough = ", tUtr, sep = "" ),
              x=unit(0.15, "npc"), y=unit(0.05, "npc"), rot=0)

}

# plotRepScan()
# Scans through BCU plots for every replicate in a simulation
# blob. Good for checking the single reps 
# inputs:   sim=number of the simulation in the project folder (alphabetical)
#           rep=optional number of replicate to start scanning from
#           est=character indicating which esitimates for plotting (MLE or MCMC)
plotRepScan <- function ( sim =1, rep = 1, est = "MLE")
{
  # First load the simulation
  .loadSim(sim)

  # Recover number of reps
  nReps <- blob$ctrl$nReps

  if (rep > nReps ) 
  {
    cat ( "(Error MSG) Starting rep too high, this sim has ",
          nReps, " replicates.", sep = "")
    return()
  }

  for (r in rep:nReps)
  {
    plotBCU(r,est)
    invisible(readline(prompt="Press [enter] for the next plot."))
  }
  return(invisible())
}

# plotMCMCpar()
# Function that will plot MCMC output from ADMB models for a nominated
# parameter, simulation and replicate. Uses the coda package.
# inputs:   rep=replicate number; 
#           par=character name of parameter
#           sim=number indicating blob to load from project dir
# output:   NULL
# usage:    post-sim, showing MCMC performance
plotMCMCpar <- function ( rep=1, par="Bmsy", sim=1 )
{
  # Check if blob is loaded, if not, load the
  if (!exists(x="blob",where=1)) .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")
  repCount  <- paste("Replicate ",rep,"/",blob$ctrl$nReps,sep="")

  # Recover blob elements for plotting
  nS    <- blob$opMod$nS
  nT    <- blob$opMod$nT

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
  plot(parMCoutSS,las=1)
  mtext(side=3,text=titleSS,outer=TRUE,padj=2)

  # create a new device and plot MS trace and posteriors
  dev.new()
  plot(parMCoutMS,las=1)
  mtext(side=3,text=titleMS,outer=TRUE,padj=2)  

  return()
}

# plotMCMCspecies()
# Function that will plot MCMC output from ADMB models for a nominated
# species and replicate from the loaded blob. Uses the coda package.
# inputs:   rep=replicate number to start from; 
#           spec=number indicating the species
#           sim=number indicating simulation to load if no blob present
# output:   NULL
# usage:    post-sim, showing MCMC performance
plotMCMCspecies <- function ( rep=1, spec=1, sim = 1 )
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  if (!exists(x="blob",where=1)) .loadSim(sim)
 
  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")
  repCount  <- paste("Replicate ",rep,"/",blob$ctl$nReps,sep="")

  # Recover blob elements for plotting
  nS    <- blob$opMod$nS
  nT    <- blob$opMod$nT

  # Species names
  if (!is.null(blob$ctrl$specNames))
  {
    specNames <- blob$ctrl$specNames
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
  parMCoutSS <- as.mcmc(parMCoutSS[1:nonNAss,])
  parMCoutMS <- as.mcmc(parMCoutMS[1:nonNAms,])

  # Now to plot each model's posterior distributions
  titleSS <- paste ( "SS MCMC for ", spec, sep = "")
  titleMS <- paste ( "MS MCMC for ", spec, sep = "")

  # plot SS trace and posteriors
  plot(parMCoutSS,las=1)
  mtext(side=3,text=titleSS,outer=TRUE,padj=2)

  # create a new device and plot MS trace and posteriors
  dev.new()
  plot(parMCoutMS,las=1)
  mtext(side=3,text=titleMS,outer=TRUE,padj=2)  

  return()
}

# plotCIbioScan()
# First load the simulation
plotCIbioScan <- function ( sim =1, rep = 1 )
{
  .loadSim(sim)

  # Recover number of reps
  nReps <- blob$ctrl$nReps

  if (rep > nReps ) 
  {
    cat ( "(Error MSG) Starting rep too high, this sim has ",
          nReps, " replicates.", sep = "")
    return()
  }

  for (r in rep:nReps)
  {
    plotCIbio( r )
    invisible( readline( prompt = "Press [enter] for the next plot." ) )
  }
  return(invisible())
}

# plotCICbio()
# Function that plots true and posterior biomass for each species for a
# given simulation and replicate.
# inputs:   rep=replicate number
#           sim=number indicating blob to load from project dir
#           quant=numeric of percentiles to be calculated
plotCIbio <- function ( rep = 1, sim=1)
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  if (!exists(x="blob",where=1)) .loadSim(sim)

  # Recover control pars
  nS <- blob$opMod$nS
  nT <- blob$opMod$nT
  species <- blob$ctrl$speciesNames[1:nS]

  # load biomass trajectories
  omBio <- blob$om$Bt[rep,,]
  msBio <- array( NA, dim = c( nS, nT, 3 ), 
                  dimnames = list( species, 1:nT, c("lCI", "val", "uCI") ) )
  ssBio <- array( NA, dim = c( nS, nT, 3 ), 
                  dimnames = list( species, 1:nT, c("lCI", "val", "uCI") ) )

  # Fill MS model biomass
  if( !is.null( blob$am$ms$CIs[[ rep ]] ) )
  {
    msCIs <- blob$am$ms$CIs[[ rep ]]
  } 
  else {
    msStdErr <- blob$am$ms$sdrep[[ rep ]]
    if( is.null( msStdErr ) ) 
      msCIs <- NA
    else {  
      msStdErr <- summary( msStdErr )
      colnames( msStdErr ) <- c( "val", "se" )
      msCIs <-  msStdErr %>%
                as.data.frame() %>%
                mutate( par = rownames(msStdErr),
                        lCI = val - 1.645*se,
                        uCI = val + 1.645*se ) %>%
                dplyr::select( par, val, se, lCI, uCI )
    }
  }
  # browser()
  if( !is.na(msCIs) )
  {
    Btrows <- which( msCIs$par == "Bt" )
    msBio[ , , 1 ] <- matrix( msCIs[ Btrows , "uCI" ], nrow = nS, ncol = nT, byrow = FALSE )
    msBio[ , , 2 ] <- matrix( msCIs[ Btrows, "val" ], nrow = nS, ncol = nT, byrow = FALSE )
    msBio[ , , 3 ] <- matrix( msCIs[ Btrows, "lCI" ], nrow = nS, ncol = nT, byrow = FALSE ) 
  }

  # Now fill SS model biomas
  if( !is.null(blob$am$ss$CIs) ) 
  {
    ssCIs <- blob$am$ss$CIs[[rep]]
  }
  else {
    # Create CIs from sd report objects
    ssCIs <- vector( mode = "list", length = 3 )
    ssSDreps <- blob$am$ss$sdrep[[rep]]
    # The following makes up for some coding errors in earlier sims
    if( length( ssSDreps ) == 1 ) k <- nS else k <- 1
    # loop over valid sd objects
    for( s in k:nS )
    {
      ssStdErr <- ssSDreps[[ s ]]
      # Save single species sd rep obj
      if( is.null(ssStdErr) ) 
        ssCIs[[ s ]] <- NA
      else {
      ssStdErr <- summary( ssStdErr )
      colnames( ssStdErr ) <- c("val","se")
      ssCIs[[ s ]] <-   msStdErr %>%
                        as.data.frame() %>%
                        mutate( par = rownames(msStdErr),
                                lCI = val - 1.645*se,
                                uCI = val + 1.645*se ) %>%
                                dplyr::select( par, val, se, lCI, uCI )    
      }
    }
  }
  for( s in 1:nS )
  {
    if( is.null( ssCIs[[ s ]] ) ) next
    if( is.na( ssCIs[[ s ]] ) ) next
    Btrows <- which( ssCIs[[ s ]]$par == "Bt" )
    ssBio[ s, , 1 ] <- matrix( ssCIs[[ s ]][ Btrows, "uCI" ], nrow = 1, ncol = nT )
    ssBio[ s, , 2 ] <- matrix( ssCIs[[ s ]][ Btrows, "val" ], nrow = 1, ncol = nT )
    ssBio[ s, , 3 ] <- matrix( ssCIs[[ s ]][ Btrows, "lCI" ], nrow = 1, ncol = nT )  
  }
  

  # Fit diagnostics (hessian)
  hpdSS <- blob$am$ss$hesspd[rep,]
  hpdMS <- blob$am$ms$hesspd[rep]

  # Set colours
  ssCol <- "steelblue"
  msCol <- "salmon"
  # Create transparent fill colours for polys
  ssFill <- col2rgb(ssCol)/255
  ssFill <- rgb(ssFill[1],ssFill[2],ssFill[3],alpha=0.5)
  msFill <- col2rgb(msCol)/255
  msFill <- rgb(msFill[1],msFill[2],msFill[3],alpha=0.5)

  # Now plot!
  par(mfrow = c(nS,1), mar = c(2,4,1,1), oma=c(4,4,1,1))
  for (s in 1:nS)
  {
    # Compute max B value
    yLim <- 2*max (omBio[s,])
    # Blank plot
    plot (  x=c(1,nT), y =c(1,yLim), type = "n", xlab = "", ylab = "Biomass (t)",
            las=1 )
    # Create polygon vertices
    xPoly <- c(1:nT,nT:1)
    yPolySS <- c(ssBio[s,1:nT,1],ssBio[s,nT:1,3])
    yPolyMS <- c(msBio[s,1:nT,1],msBio[s,nT:1,3])
    if (s == 1) panLegend ( x=0.2,y=1,legTxt=c("ss","ms"),
                            fill=c(ssFill,msFill), border=c(NA,NA), 
                            cex=c(0.7), bty="n" )
    if ( hpdSS[s] ) panLab (x=0.9,y=0.9,txt="h",col=ssCol,cex=0.7)
    if ( hpdMS ) panLab (x=0.9,y=0.85,txt="h",col=msCol,cex=0.7)
    # plot CI bands
    # SS
    polygon (x = xPoly, y = yPolySS,border=NA,col = ssFill,density=NA)
    # MS
    polygon (x = xPoly, y = yPolyMS,border=NA,col = msFill,density=NA)
    # Plot Medians
    lines ( x = 1:nT, y = ssBio[s,,2], col = ssCol,lty=2,lwd=2 )
    lines ( x = 1:nT, y = msBio[s,,2], col = msCol,lty=2,lwd=2 )
    # Plot true biomass
    lines ( x=1:nT, y = omBio[s,],col="black", lwd = 2)
  }
}

# plotBCU()
# Function that will open a supplied saved blob object and plot a 
# given replicate's true and estimated time series of biomass, 
# catch, exploitation and IoA data. 
# inputs:   rep=replicate number for plotting; 
#           sim=number indicating blob to load (alpha order in project folder)
#           folder=name of folder/blob file (supercedes sim number)
# output:   NULL
# usage:    post-simulation run, plotting performance
plotBCU <- function ( rep = 1, est="MLE", sim=1, legend=TRUE,
                      data = FALSE, labSize = 2, tickSize = 1.2 )
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  if (!exists(x="blob",where=1)) .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")
  repCount  <- paste("Replicate ",rep,"/",blob$ctrl$nReps,sep="")

  # Recover blob elements for plotting
  nS <- blob$opMod$nS
  nT <- blob$opMod$nT

  # Species names
  specNames <- blob$ctrl$speciesNames

  # True OM quantities
  omBt  <- blob$om$Bt[rep,,]
  Ct    <- blob$om$Ct[rep,,]
  It    <- blob$om$It[rep,,]
  Ut    <- blob$om$Ut[rep,,]

  
  
  if ( est == "MLE" )
  { # Single species model
    ssBt  <- blob$am$ss$Bt[rep,,]
    ssq   <- blob$am$ss$q[rep,]

    # Multispecies model
    msBt  <- blob$am$ms$Bt[rep,,]
    msq   <- blob$am$ms$q[rep,]  
  }

  # Estimated Ut
  ssUt <- Ct / ssBt
  msUt <- Ct / msBt

  
  # Recover diagnostics for the fits
  hpdSS <- blob$am$ss$hesspd[rep,]
  grdSS <- blob$am$ss$maxGrad[rep,]
  hpdMS <- blob$am$ms$hesspd[rep]
  grdMS <- blob$am$ms$maxGrad[rep]

  # Set colours for each model
  ssCol <- "steelblue"
  msCol <- "salmon"

  # Set up plot window
  par ( mfrow = c(3,nS), mar = c(1,4.5,2,0), oma = c(3,1,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )
  # Plot biomass, actual and estimated, including 2 index series,
  # scaled by model estimated q
  for ( s in 1:nS )
  {
    if ( s == 1 ) yLab <- "Biomass (t)" else yLab <- ""
    maxBt <- 1.2*max ( omBt[s,],na.rm=TRUE)
    plot    ( x = c(1,nT), y = c(0,maxBt), type = "n",
              ylim = c(0,maxBt), ylab = yLab, las = 1, xlab = "" ,
              main = specNames[s] )
    if (s == 1) panLegend ( x=0.2,y=1,legTxt=c("ss","ms"),
                            col=c(ssCol,msCol), lty = c(2,2), 
                            lwd = c(2,2), cex=c(1), bty="n" )
    # if ( minSS[s] ) panLab (x=0.85,y=0.9,txt="c",col=ssCol,cex=1.1)
    if ( hpdSS[s] & !is.na(hpdSS[s]) ) panLab (x=0.9,y=0.9,txt="h",col=ssCol,cex=1.1)
    # if ( minMS ) panLab (x=0.85,y=0.85,txt="c",col=msCol,cex=1.1)
    if ( hpdMS & !is.na(hpdMS) ) panLab (x=0.9,y=0.85,txt="h",col=msCol,cex=1.1)
    if ( data ) points  ( x = 1:nT, y = It[s,]/ssq[s], pch = 2, cex = 0.6, col="grey70" )
    if ( data ) points  ( x = 1:nT, y = It[s,]/msq[s], pch = 5, cex = 0.6, col="grey70" )
    lines   ( x = 1:nT, y = omBt[s,], col = "black", lwd = 2)
    lines   ( x = 1:nT, y = ssBt[s,], col = ssCol, lwd = 2, lty = 2 )
    lines   ( x = 1:nT, y = msBt[s,], col = msCol, lwd = 2, lty = 2 )
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
    if ( s == 1 ) yLab <- "Exploitation Rate" else yLab <- ""
    maxUt <- max ( Ut[s,])
    plot  ( x = 1:nT, y = Ut[s,], col = "black", lwd = 2, type = "l",
            ylim = c(0,maxUt), ylab = yLab, las = 1, xlab = "" )
    lines (x = 1:nT, y = ssUt[s,], col= ssCol, lwd = 2, lty = 2 )
    lines (x = 1:nT, y = msUt[s,], col= msCol, lwd = 2, lty = 2 )
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
plotSimPerf <- function ( pars = c("Bmsy","Umsy","q","dep","BnT","tau2","totRE"), sim=1, 
                          est="MLE" )
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  if (!exists(x="blob",where=1)) .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")

  # Recover blob elements for plotting
  nS        <- blob$opMod$nS
  specNames <- blob$ctrl$specNames[1:nS]

  # Recover relative error distributions
  if (est == "MLE")
  {
    ssRE <- blob$am$ss$err.mle[pars]
    msRE <- blob$am$ms$err.mle[pars]
  }
  
  # Create a wrapper function for generating quantiles
  quantWrap <- function ( entry = 1, x = ssRE, ... )
  {
    if (is.null(dim(x[[entry]]))) x[[entry]] <- matrix(x[[entry]],ncol=1)
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
    for (p in 1:length(pars))
    {
      quantiles[s,,p,1] <- ssQuant[pars[p]][[1]][s,]
      quantiles[s,,p,2] <- msQuant[pars[p]][[1]][s,]
    }
  }
  # Set plotting window
  par (mfrow = c(1,nS), mar = c(3,0,1,1), oma = c(3,0,3,0) )

  for ( s in 1:nS )
  {
    med <- quantiles[s,"50%",,]
    q975 <- quantiles[s,"97.5%",,]
    q025 <- quantiles[s,"2.5%",,]

    # Plot main dotchart
    plotTitle <- specNames[s]
    if ( s == 1) dotchart ( x = t(med), xlim = c(-1.5,1.5), main=plotTitle,
                            pch = 16)
    else dotchart ( x = t(med), xlim = c(-1.5,1.5), main = plotTitle,
                    pch = 16)
    
    # Now add segments
    for ( p in length(pars):1)
    {
      parIdx <- length(pars) - p + 1
      plotY <- 2 * p + 2 * (p-1)
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


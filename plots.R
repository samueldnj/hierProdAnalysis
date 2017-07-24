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

plotStatTableGraphs <- function(  tableRoot = "allSame_infoScenarios",
                                  resp = c("BnT","Umsy","Bmsy","Dep","HessPD","q_1","q_2","tau2_1","tau2_2"),
                                  axes = c("fYear","nDiff"),
                                  MARE = TRUE,
                                  RE = TRUE,
                                  groupPars = FALSE )
{
  # plotStatTableGraphs()
  # Plots the statistics table information and saves into
  # the figs directory.
  # inputs:   tableRoot = RE tableName root in statistics folder
  #           resp = response variables for plotting
  #           axes = 2-character giving factors to plot over
  # outputs:  NULL
  # side-eff: plots to quartz or figs directory

  # Create save directory
  saveDir <- file.path(getwd(),"project","figs",tableRoot)
  dir.create(saveDir)

  # First, plot MARE values
  # a place to save them
  MAREdir <- file.path(saveDir,"MARE")
  dir.create(MAREdir)
  # Now load the RE table
  MAREtabRoot  <- paste(tableRoot,"_MARE", sep = "")
  MAREtabFile  <- paste(tableRoot,"_MARE.csv", sep = "")
  MAREtabPath <- file.path(getwd(),"project","Statistics",MAREtabFile)
  MAREtab     <- read.csv( MAREtabPath, header=TRUE, stringsAsFactors=FALSE )
  nSpp        <- unique(MAREtab$nS)
  MPs         <- unique(MAREtab$mp)



  # Loop over complex sizes
  if( MARE )
  {
    for( sIdx in 1:length(nSpp) )
    {
      # Make folder for this complex size
      nSp           <- nSpp[sIdx]
      nSpFolderName <- paste(nSp,"stocks",sep = "_")
      nSpFolder     <- file.path(MAREdir,nSpFolderName)
      dir.create(nSpFolder)
      nSppMAREtab   <- MAREtab %>% filter( nS == nSp )
      spp           <- unique(nSppMAREtab$species) 

      for( rIdx in 1:length(resp) )
      {
        response = resp[rIdx]
        plotFile      <- paste(response,"_",nSp,"stocks_",MAREtabRoot,"_byMP.pdf", sep = "") 
        plotPath      <- file.path(nSpFolder,plotFile)
        pdf( file = plotPath, width = 11, height = 17 )
        plotMAREs(  tableName = MAREtabRoot, table = nSppMAREtab,
                    axes = c(axes,"mp"), nSp = nSp, spec = "Stock1", resp = response )
        dev.off()
      }
    }
  }


  

  # Now plot relative error distributions
  REdir <- file.path(saveDir,"RE")
  dir.create(REdir)
  # Now load the RE table
  REtabRoot   <- paste(tableRoot,"_RE", sep = "")
  REtabFile   <- paste(tableRoot,"_RE.csv", sep = "")
  REtabPath   <- file.path(getwd(),"project","Statistics",REtabFile)
  if( RE )
  {
    REtab       <- read.csv( REtabPath, header=TRUE, stringsAsFactors=FALSE )
    nSpp        <- unique(REtab$nS)

    # Loop over complex sizes
    for( sIdx in 1:length(nSpp) )
    {
      # Make folder for this complex size
      nSp           <- nSpp[sIdx]
      nSpFolderName <- paste(nSp,"stocks",sep = "_")
      nSpFolder     <- file.path(REdir,nSpFolderName)
      dir.create(nSpFolder)
      nSppREtab   <- REtab %>% filter( nS == nSp )
      spp           <- unique(nSppREtab$species) 

      for( rIdx in 1:length(resp) )
      {
        response = resp[rIdx]
        plotFile      <- paste(response,"_",nSp,"stocks_",REtabRoot,"_byMP.pdf", sep = "") 
        plotPath      <- file.path(nSpFolder,plotFile)
        pdf( file = plotPath, width = 11, height = 17 )
        plotREdists(  tableName = REtabRoot, 
                      table = nSppREtab,
                      axes = c(axes,"mp"), 
                      nSp = nSp, 
                      spec = "Stock1", 
                      resp = response )
        dev.off()
      }
    }
  }
  
  if( groupPars )
  {
    # Now plot group parameters
    MREtabRoot <- paste(tableRoot,"_MRE", sep = "" )
    for( nSp in nSpp )
    {
      for( mIdx in 1:length(MPs))
      {
        MP <- MPs[mIdx]
        plotFile <- paste( nSp, "nS_", MP,"_groupPars_", tableRoot,".pdf", sep = "" ) 
        plotPath <- file.path( saveDir, plotFile )
        pdf( file = plotPath, width = 11, height = 17 )
        plotGroupPars(  tableName = MREtabRoot,
                        axes = axes,
                        nSp = nSp,
                        MP = MP,
                        spec = "Stock1"
                     )  
        dev.off()
      } 
    }  
  }
  cat("Table plotting complete for ", tableRoot, ".\n", sep = "")
}



plotGroupPars <- function(  tableName = "rKq_msInc__MRE",
                            axes = c("qOM","UmsyOM"),
                            nSp = 10,
                            MP = "allJointPriors",
                            spec = "Stock1" )
{
  # plotGroupPars()
  # Plots group prior parameter estimates 
  # as a function of axes[1], grouped by axes[2].
  # inputs:   tableName = RE table in statistics folder
  #           nSp = number of species to restrict to
  #           axes = 2-numeric giving factor levels to plot over
  #           par = character name of group level parameter to plot
  # outputs:  NULL
  # side-eff: plots to quartz

  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE ) 

  if(!is.null(spec)) table <- table %>% filter(species == spec)

  # restrict to correct number of species, summarise group pars since we're
  # only using 2 axes
  table  <-   table %>%
              filter( mp == MP, nS == nSp, species == spec ) %>%
              group_by( kappaMult, corr, Umax, tUpeak ) %>%
              summarise_if(.predicate = "is.numeric",.funs="median") %>%
              ungroup()


  # Get levels for axes
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }
  # order levels
  x   <- unique(table %>% pull(axesCols[1]))
  x   <- x[order(x)]
  y   <- unique(table %>% pull(axesCols[2]))
  y   <- y[order(y)]

  groupPars <- c("Umsybar","tauq2","sigU2")
  
  # Create an array to hold info
  z <- array( NA, dim = c( length(x), length(y), length(groupPars), 3 ),
              dimnames = list(  paste(axes[1],x,sep = ""),
                                paste(axes[2],y,sep = ""),
                                groupPars,
                                c(0.025, 0.5, 0.975) ) )

  # Fill the array, looping over species and axes
  for( xIdx in 1:length(x) )
  {
    for( yIdx in 1:length(y) )
    {
      # Get x and y values
      xVal <- x[ xIdx ]
      yVal <- y[ yIdx ]
      # Get rows
      xRows <- which ( table[ ,axesCols[ 1 ] ] == xVal )
      yRows <- which ( table[ ,axesCols[ 2 ] ] == yVal )
      rows <- intersect( xRows, yRows )
      # browser()
      if( length(rows) == 0 ) next
      filtered <- table[rows,]
      if(nrow(filtered) > 1) filtered <- filtered[1,]
      # Fill array
      z[xIdx,yIdx,,2] <- as.numeric(filtered[1,groupPars])
      z[xIdx,yIdx,,1] <- as.numeric(filtered[1,paste(groupPars,"025",sep = "")])
      z[xIdx,yIdx,,3] <- as.numeric(filtered[1,paste(groupPars,"975",sep = "")])
    }
  }
  # Set up plotting colours and a jitter for different axes
  cols      <- brewer.pal( n = length(y), name = "Dark2" )
  xShift    <- 1/(length(y)+2)
  leftShift <- xShift*length(y)/2
  # Set up plotting environment
  par(mfrow = c(2,2), mar = c(2,2,1,1), oma = c(3,3,2,2) )
  # Loop over the group parameters
  for( p in 1:length(groupPars) )
  {
    plot( x = c(1-leftShift,length(x) + leftShift), 
          y = range(z[,,groupPars[p],], na.rm=T), axes = F, type = "n",
          xlab = "", ylab = "", main = groupPars[p] )
    axis( side = 1, at = 1:length(x), labels = x )
    axis( side = 2, las = 1 )
    for( yIdx in 1:length(y) )
    {
      lines(  x = 1:length(x) - leftShift + xShift*(yIdx-1), 
              y = z[,yIdx,groupPars[p],2], col = cols[yIdx], lwd = 0.5 )
      points( x = 1:length(x) - leftShift + xShift*(yIdx-1), 
              y = z[,yIdx,groupPars[p],2], col = cols[yIdx], pch = 16 )
      segments( x0 = 1:length(x) - leftShift + xShift*(yIdx-1),
                x1 = 1:length(x) - leftShift + xShift*(yIdx-1),
                y0 = z[,yIdx,groupPars[p],1], 
                y1 = z[,yIdx,groupPars[p],3], col = cols[yIdx] )
    }
  }
  
  panLegend(  legTxt = paste( axes[2], " = ", y ),
              lty = 1, pch = 16, cex = 1,
              lwd = 0.8, col = cols,
              x = 0.1, y = 0.9, bty = "n" )
  mtext( side = 1, outer = TRUE, text = axes[1], line = 1.5 )  
}

plotMAREs <- function(  tableName = "allSame_RE_msIncr__MARE",
                        axes = c("corr","kappaMult","mp"),
                        resp = "Umsy",
                        spec = "Stock1",
                        MP = NULL,
                        nSp = 6,
                        table = NULL,
                        ssModel = NULL,
                        msModel = NULL
                      )
{
  # plotMAREs()
  # Plots MARE/MRE (response variable) values as a function
  # of given axes (depends on experimental factors). 
  # Will group points using lines corresponding to the first axis (w) and 
  # plot across values of the second axes (x), with multi-panel rows 
  # corresponding to the third (y). Panel columns are 
  # single and Joint AMs.
  # inputs:   tableName = RE table in statistics folder
  #           axes = 3-char vector naming RE table columns
  #           resp = response variable of interest for RE values 
  #                  (usually a leading or derived AM parameter)
  #           table = optional df to replace loading file, useful for 
  #                   automated plotting to directory
  # outputs:  NULL

  # Read in table if none supplied
  if( is.null(table) )
  {
    # Load table
    fileName <- paste( tableName, ".csv", sep = "" )
    tablePath <- file.path ( getwd(), "project/Statistics", fileName )
    table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE ) 
  }
  tab <- table
  # Restrict to subtables if requested
  if(!is.null(MP))    tab  <-   tab %>% filter( mp == MP ) 
  if(!is.null(spec))  tab  <-   tab %>% filter( species == spec ) 
  if(!is.null(nSp))   tab  <-   tab %>% filter( nS == nSp )

  # # Run predictions if explanatory model supplied
  # predTab <- tab %>% mutate( ssPred = NA, msPred = NA )
  # if( !is.null(ssModel) )
  # {
  #   predTab <-  predTab %>% 
  #               mutate( ssPred = predict( ssModel, newdata = predTab ) )
  # }

  # if( !is.null(msModel) )
  # {
  #   predTab <-  predTab %>% 
  #               mutate( msPred = predict( msModel, newdata = predTab ) )
  # }

  # response colnames
  respCols <- paste(c("ss","ms"),resp,sep = "")


  # Get levels for axes
  axesCols <- integer(length=length(axes))
  # predAxesCols <- integer(length=length(axes))
  for ( i in 1:3)
  {
    axesCols[i]     <- which(names(tab) == axes[i])
    # predAxesCols[i] <- which(names(predTab) == axes[i])
  }
  # order levels
  w   <- unique(tab %>% pull(axesCols[1]))
  w   <- w[order(w)]
  x   <- unique(tab %>% pull(axesCols[2]))
  x   <- x[order(x)]
  y   <- unique(tab %>% pull(axesCols[3]))
  y   <- y[order(y)]

  # Create an array to hold info
  z <- array( NA, dim = c( length(w), length(x), length(y), 2 ),
              dimnames = list(  paste(axes[1],w,sep = ""),
                                paste(axes[2],x,sep = ""),
                                paste(axes[3],y,sep = ""),
                                respCols ) )

  # Create an array to hold info
  zPred <- array( NA, dim = c( length(w), length(x), length(y), 2 ),
                  dimnames = list(  paste(axes[1],w,sep = ""),
                                    paste(axes[2],x,sep = ""),
                                    paste(axes[3],y,sep = ""),
                                    respCols ) )


  # Fill the array, looping over species and axes
  for( wIdx in 1:length(w) )
  {
    for( xIdx in 1:length(x) )
    {
      for( yIdx in 1:length(y) )
      {
        # Get x and y values
        wVal <- w[ wIdx ]
        xVal <- x[ xIdx ]
        yVal <- y[ yIdx ]
        # Get rows
        wRows <- which ( tab[,axesCols[1] ] == wVal)
        xRows <- which ( tab[,axesCols[2] ] == xVal)
        yRows <- which ( tab[,axesCols[3] ] == yVal)
        # wRowsPred <- which ( predTab[,predAxesCols[1] ] == wVal)
        # xRowsPred <- which ( predTab[,predAxesCols[2] ] == xVal)
        # yRowsPred <- which ( predTab[,predAxesCols[3] ] == yVal)
        rows      <- intersect(intersect(wRows,xRows),yRows)
        # predRows  <- intersect(intersect(wRowsPred,xRowsPred),yRowsPred)
        if( length(rows) == 0 ) next
        filtered <- tab[rows,respCols]
        if( nrow(filtered) > 1 ) 
          filtered <-  apply(X = filtered, FUN = mean, MARGIN = 2, na.rm = T)
        # Fill array
        z[wIdx,xIdx,yIdx,respCols] <- as.numeric( filtered )
        # Now fill zPred array  
        # predTab <- predTab[predRows,c("ssPred","msPred")]
        # if(any(!is.na(predTab)))
        # {
        #   if( nrow(predTab > 1) ) predTab <- apply(X = predTab, FUN = mean, MARGIN = 2, na.rm = T)
        #   zPred[wIdx,xIdx,yIdx,respCols] <- as.numeric( predTab[ , c( "ssPred", "msPred" ) ] )
        # }  
      }
    }
  }
  if(any(grepl(pattern="msHessPD", x = respCols)) ) z[,,,"msHessPD"] <- z[,,,"msHessPD"] + 100

  zDelta <- log2(z[,,,1] / (z[,,,2]))

  cols <- brewer.pal( n = length(w), name = "Dark2" )
  xShift <- 1 / (length(w) + 2)
  leftShift <- length(w)*xShift / 2

  
  # Plot contents. Use w values as line groupings, x values as
  # horizontal axis points, and y values as panel row groups
  par( mfrow = c(length(y), 3), mar = c(2,2,2,1), oma = c(3,3,2,2) )
  for( yIdx in 1:length(y) )
  {
    # Plot SS model
    plot( x = c(1 - leftShift,length(x)+leftShift), 
          y = c(0,max(z, na.rm = T )), type = "n",
          axes = FALSE )
      if( yIdx == 1 ) mtext( side = 3, text = "Single Species Model", cex = 0.8)
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      for( wIdx in 1:length(w) )
      {
        rect( xleft = 1:length(x) - leftShift + xShift*(wIdx-1), ybottom = 0, 
              xright = 1:length(x) - leftShift + xShift*(wIdx), ytop = z[ wIdx, , yIdx, respCols[1] ],
              col = cols[wIdx] )
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[3], " = ", y[yIdx], sep = "" ) )
    # plot MS model
    plot( x = c(1 - leftShift,length(x)+leftShift), 
          y = c(0,max(z,na.rm = T)), type = "n",
          axes = FALSE )
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      if( yIdx == 1 ) mtext( side = 3, text = "Joint Model", cex = 0.8)
      for( wIdx in 1:length(w) )
      {
        rect( xleft = 1:length(x) + xShift*(wIdx-1) - leftShift, ybottom = 0, 
              xright = 1:length(x) + xShift*(wIdx) - leftShift , ytop = z[ wIdx, , yIdx, respCols[2] ],
              col = cols[wIdx] )
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[3], " = ", y[yIdx], sep = "" ) )

    # plot Delta statistic
    plot( x = c(1 - leftShift,length(x)+leftShift), 
          y = range(zDelta,na.rm = T), type = "n",
          axes = FALSE )
      if( yIdx == 1 ) panLegend(  legTxt = paste( axes[1], " = ", w ),
                                  lty = 1, pch = 16, cex = 1,
                                  lwd = 0.8, col = cols,
                                  x = 0.1, y = 0.7, bty = "n" )
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      if( yIdx == 1 ) mtext( side = 3, text = "Delta Values", cex = 0.8)
      for( wIdx in 1:length(w) )
      {
        xShift <- 1/(length(w) + 2)
        rect( xleft = 1:length(x) + xShift*(wIdx-1) - leftShift, ybottom = 0, 
              xright = 1:length(x) + xShift*(wIdx) - leftShift , ytop = zDelta[ wIdx, , yIdx ],
              col = cols[wIdx] )
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[3], " = ", y[yIdx], sep = "" ) )
  }
  mtext( side = 3, outer = T, text = tableName, cex = 1.2 )
  mtext( side = 1, outer =T, text = axes[2], cex = 1.2 )
  mtext( side = 2, outer = T, text = "Relative Error (%)", cex = 1.2, line = 1 )
  mtext( side = 1, outer = T, text = resp, cex = 1, col = "grey50",
          adj = 0.9, line = 1.3 )
  mtext( side = 1, outer = T, text = spec, cex = 1, col = "grey50",
          adj = 0.1, line = 1.3 )
}

plotREdists <- function(  tableName = "allSame_RE_msIncr__RE",
                          axes = c("corr","kappaMult","mp"),
                          resp = "Umsy",
                          spec = "Stock1",
                          MP = NULL,
                          nSp = 6,
                          table = NULL,
                          ssModel = NULL,
                          msModel = NULL
                        )
{
  # plotREdists()
  # Plots relative error distributions as  a function
  # of given axes (depends on experimental factors). 
  # Will group points using lines corresponding to the first axis (w) and 
  # plot across values of the second axes (x), with multi-panel rows 
  # corresponding to the third (y). Panel columns are 
  # single and Joint AMs.
  # inputs:   tableName = RE table in statistics folder
  #           axes = 3-char vector naming RE table columns
  #           resp = response variable of interest for RE values 
  #                  (usually a leading or derived AM parameter)
  #           table = optional df to replace loading file, useful for 
  #                   automated plotting to directory
  # outputs:  NULL

  # Read in table if none supplied
  if( is.null(table) )
  {
    # Load table
    fileName <- paste( tableName, ".csv", sep = "" )
    tablePath <- file.path ( getwd(), "project/Statistics", fileName )
    table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE ) 
  }
  tab <- table
  # Restrict to subtables if requested
  if(!is.null(MP))    tab  <-   tab %>% filter( mp == MP ) 
  if(!is.null(spec))  tab  <-   tab %>% filter( species == spec ) 
  if(!is.null(nSp))   tab  <-   tab %>% filter( nS == nSp )

  # Run predictions if explanatory model supplied
  predTab <- tab %>% mutate( ssPred = NA, msPred = NA )
  if( !is.null(ssModel) )
  {
    predTab <-  predTab %>% 
                mutate( ssPred = predict( ssModel, newdata = predTab ) )
  }

  if( !is.null(msModel) )
  {
    predTab <-  predTab %>% 
                mutate( msPred = predict( msModel, newdata = predTab ) )
  }

  # response colnames
  respCols <- paste(c("ss","ms"),resp,sep = "")


  # Get levels for axes
  axesCols <- integer(length=length(axes))
  predAxesCols <- integer(length=length(axes))
  for ( i in 1:3)
  {
    axesCols[i]     <- which(names(tab) == axes[i])
    predAxesCols[i] <- which(names(predTab) == axes[i])
  }
  # order levels
  w   <- unique(tab %>% pull(axesCols[1]))
  w   <- w[order(w)]
  x   <- unique(tab %>% pull(axesCols[2]))
  x   <- x[order(x)]
  y   <- unique(tab %>% pull(axesCols[3]))
  y   <- y[order(y)]

  # Create an array to hold info
  z <- array( NA, dim = c( length(w), length(x), length(y), 2, 3 ),
              dimnames = list(  paste(axes[1],w,sep = ""),
                                paste(axes[2],x,sep = ""),
                                paste(axes[3],y,sep = ""),
                                respCols, c("0.025","0.5","0.975") ) )

  # Create an array to hold info
  zPred <- array( NA, dim = c( length(w), length(x), length(y), 2 ),
                  dimnames = list(  paste(axes[1],w,sep = ""),
                                    paste(axes[2],x,sep = ""),
                                    paste(axes[3],y,sep = ""),
                                    respCols ) )

  # Fill the array, looping over species and axes
  for( wIdx in 1:length(w) )
  {
    for( xIdx in 1:length(x) )
    {
      for( yIdx in 1:length(y) )
      {
        # Get x and y values
        wVal <- w[ wIdx ]
        xVal <- x[ xIdx ]
        yVal <- y[ yIdx ]
        # Get rows
        wRows <- which ( tab[,axesCols[1] ] == wVal)
        xRows <- which ( tab[,axesCols[2] ] == xVal)
        yRows <- which ( tab[,axesCols[3] ] == yVal)
        wRowsPred <- which ( predTab[ , predAxesCols[ 1 ] ] == wVal)
        xRowsPred <- which ( predTab[ , predAxesCols[ 2 ] ] == xVal)
        yRowsPred <- which ( predTab[ , predAxesCols[ 3 ] ] == yVal)
        rows      <- intersect(intersect(wRows,xRows),yRows)
        predRows  <- intersect(intersect(wRowsPred,xRowsPred),yRowsPred)
        if( length(rows) == 0 ) next
        filtered <- tab[rows,respCols]
        quantiles <- apply( X = filtered, FUN = quantile, 
                            probs = c(0.025,0.5, 0.975), MARGIN = 2,
                            na.rm = T)
        medPred <- apply( X = predTab[predRows,c("ssPred","msPred")],
                          FUN = median, na.rm = T, MARGIN = 2 )
        # Fill array
        # browser()
        z[wIdx,xIdx,yIdx,,] <- t(quantiles)
        zPred[wIdx,xIdx,yIdx,] <- as.numeric(medPred)
      }
    }
  }

  cols <- brewer.pal( n = length(w), name = "Dark2" )
  xShift <- 1 / (length(w) + 2)
  leftShift <- length(w)*xShift / 2
  
  # Plot contents. Use w values as line groupings, x values as
  # horizontal axis points, and y values as panel row groups
  par( mfrow = c(length(y), 2), mar = c(2,2,2,1), oma = c(3,3,2,2) )
  for( yIdx in 1:length(y) )
  {
    # Plot SS model
    plot( x = c(1-leftShift,length(x)+leftShift), 
          y = range(z, na.rm = T ), type = "n",
          axes = FALSE )
      if( yIdx == 1 ) panLegend(  legTxt = paste( axes[1], " = ", w ),
                                  lty = 1, pch = 16, cex = 1,
                                  lwd = 0.8, col = cols,
                                  x = 0.1, y = 0.7, bty = "n" )
      if( yIdx == 1 ) mtext( side = 3, text = "Single Species Model", cex = 0.8)
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      # axis( side = 3, at = 1:length(x), labels = FALSE )
      abline(h = 0, lty = 2, lwd = 0.8)
      for( wIdx in 1:length(w) )
      {
        lines(  x = 1:length(x) - leftShift + xShift*(wIdx-1), y = z[ wIdx, , yIdx, respCols[1], 2], col = cols[wIdx],
                lwd = 0.6, lty = 3 )
        points( x = 1:length(x) - leftShift + xShift*(wIdx-1), y = z[ wIdx, , yIdx, respCols[1], 2 ], col = cols[wIdx],
                pch = 16, cex = 1 )
        segments( x0 = 1:length(x) - leftShift + xShift*(wIdx-1), y0 = z[ wIdx, , yIdx, respCols[1], 1],
                  x1 = 1:length(x) - leftShift + xShift*(wIdx-1), y1 = z[ wIdx, , yIdx, respCols[1], 3],
                  col = cols[wIdx], lwd = 2 )
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[3], " = ", y[yIdx], sep = "" ) )
    # plot MS model
    plot( x = c(1-leftShift,length(x)+leftShift), 
          y = range(z,na.rm = T), type = "n",
          axes = FALSE )
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      # axis( side = 3, at = 1:length(x), labels = FALSE)
      abline(h = 0, lty = 2, lwd = 0.8)
      if( yIdx == 1 ) mtext( side = 3, text = "Joint Model", cex = 0.8)
      for( wIdx in 1:length(w) )
      {
        lines(  x = 1:length(x) - leftShift + xShift*(wIdx-1), y = z[ wIdx, , yIdx, respCols[2], 2], col = cols[wIdx],
                lwd = 0.6, lty = 3 )
        points( x = 1:length(x) - leftShift + xShift*(wIdx-1), y = z[ wIdx, , yIdx, respCols[2], 2 ], col = cols[wIdx],
                pch = 16, cex = 1 )
        segments( x0 = 1:length(x) - leftShift + xShift*(wIdx-1), y0 = z[ wIdx, , yIdx, respCols[2], 1],
                  x1 = 1:length(x) - leftShift + xShift*(wIdx-1), y1 = z[ wIdx, , yIdx, respCols[2], 3],
                  col = cols[wIdx], lwd = 2 )
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[3], " = ", y[yIdx], sep = "" ) )
  }
  mtext( side = 3, outer = T, text = tableName, cex = 1.2 )
  mtext( side = 1, outer =T, text = axes[2], cex = 1.2 )
  mtext( side = 2, outer = T, text = "Relative Error (%)", cex = 1.2, line = 1 )
  mtext( side = 1, outer = T, text = resp, cex = 1, col = "grey50",
          adj = 0.9, line = 1.3 )
  mtext( side = 1, outer = T, text = MP, cex = 1, col = "grey50",
          adj = 0.1, line = 1.3 )
}

plotREtableSlice <- function( tableName = "Fhist_pub_MARE",
                              nSp = 5,
                              axes = c("tUpeak","Umax"),
                              specVals = c( Dover = 0.30,
                                            English = 0.40,
                                            Rock = 0.35,
                                            Petrale = 0.25,
                                            Arrowtooth = 0.15 ),
                              par = "Umsy"
                            )
{
  # plotREtableSlice()
  # Plots MARE/MRE (response variable) values as a function
  # of given specVals (OM values of q, Umsy, etc). 
  # Will group lines corresponding to the first axis and create panel rows 
  # corresponding to the second. Panel columns are single and Joint
  # AMs.
  # inputs:   tableName = RE table in statistics folder
  #           nSp = number of species to restrict to
  #           axes = 2-char vector naming RE table columns
  #           specVals = named numeric of OM quantities to use as plot x axis
  #           par = estimated parameter of interest for MRE/MARE values
  # outputs:  NULL
  # side-eff: plots to quartz

  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE ) 

  # restrict to correct number of species
  table  <-   table %>%
              filter( nS == nSp )

  specVals <- specVals[1:nSp]
  specVals <- specVals[order(specVals)]


  # Get levels for axes
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }
  # order levels
  x   <- unique(table %>% pull(axesCols[1]))
  x   <- x[order(x)]
  y   <- unique(table %>% pull(axesCols[2]))
  y   <- y[order(y)]

  mps <- unique( table %>% pull(mp) )

  # Create an array to hold info
  z <- array( NA, dim = c( length(x), length(y), length(specVals), 2, length(mps) ),
              dimnames = list( x, y, names(specVals), c("ss","ms"), mps ) )

  # Fill the array, looping over species and axes
  for(sIdx in 1:length(specVals) )
  {
    spec <- names( specVals )[ sIdx ]
    specTab <- table %>% filter( species == spec )
    for( xIdx in 1:length(x) )
    {
      for( yIdx in 1:length(y) )
      {
        # Get x and y values
        xVal <- x[ xIdx ]
        yVal <- y[ yIdx ]
        # Get rows
        xRows <- which ( specTab[,axesCols[1] ] == xVal)
        yRows <- which ( specTab[,axesCols[2] ] == yVal)
        rows <- intersect(xRows,yRows)
        if( length(rows) == 0 ) next
        # Fill array
        parColNames <- paste(c("ss","ms"),par, sep = "")
        z[xIdx,yIdx,sIdx,,] <- as.numeric(specTab[rows,parColNames])
      }
    }
  }

  cols <- brewer.pal( n = length(x), name = "Dark2" )
  
  # Plot contents. Use x values as line groupings, y values as
  # panel groupings, and species productivities as x values
  par( mfrow = c(length(y), 2), mar = c(2,2,2,1), oma = c(3,3,2,2) )
  for( yIdx in 1:length(y) )
  {
    # Plot SS model
    plot( x = c(1,length(specVals)), 
          y = c( min( z ), max( z ) ), type = "n",
          axes = FALSE )
      if( yIdx == 1 ) panLegend(  legTxt = paste( axes[1], " = ", x ),
                                  lty = 1, pch = 16, cex = 1,
                                  lwd = 0.8, col = cols,
                                  x = 0.1, y = 0.7, bty = "n" )
      if( yIdx == 1 ) mtext( side = 3, text = "Single Species Model", cex = 0.8)
      axis( side = 1, at = 1:length(specVals),labels = specVals )
      axis( side = 2, las = 1 )
      axis( side = 3, at = 1:length(specVals), labels = names(specVals))
      abline(h = 0, lty = 2, lwd = 0.8)
      for( m in 1:length(mps) )
      {
        for( xIdx in 1:length(x) )
              {
                lines(  x = 1:length(specVals), y = z[ xIdx, yIdx, , "ss", m ], col = cols[xIdx],
                        lwd = 0.8, lty = m )
                points( x = 1:length(specVals), y = z[ xIdx, yIdx, , "ss", m ], col = cols[xIdx],
                        pch = 16, cex = 0.5 )
              }
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[2], " = ", y[yIdx], sep = "" ) )

    plot( x = c(1,length(specVals)), 
          y = c( min( z ), max( z ) ), type = "n",
          axes = FALSE )
      axis( side = 1, at = 1:length(specVals),labels = specVals )
      axis( side = 2, las = 1 )
      axis( side = 3, at = 1:length(specVals), labels = names(specVals))
      abline(h = 0, lty = 2, lwd = 0.8)
      if( yIdx == 1 ) mtext( side = 3, text = "Joint Model", cex = 0.8)
      for( m in 1:length(mps) )
      {
        for( xIdx in 1:length(x) )
        {
          lines(  x = 1:length(specVals), y = z[ xIdx, yIdx, , "ms", m ], col = cols[xIdx],
                  lwd = 0.8, lty = m )
          points( x = 1:length(specVals), y = z[ xIdx, yIdx, , "ms", m ], col = cols[xIdx],
                  pch = 16, cex = 0.5 )
        }  
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[2], " = ", y[yIdx], sep = "" ) )
  }
  mtext( side = 3, outer = T, text = tableName, cex = 1.2 )
  mtext( side = 1, outer = T, text = par, cex = 1, col = "grey80",
          adj = 0.9, line = 1.5 )
}

plotFhist <- function ( sims = 1:2,
                        nTraj = 40,
                        seed = 2,
                        save = FALSE,
                        fileName = "depPlot.pdf",
                        ... )
{
  # plotFhist()
  # This function plots traces of biomass trajectories
  # for different fishing histories. 
  # This one's a little path specific, and only needs
  # to be pointed to the correct sims
  # inputs      sims=numeric vector indicating sims in project dir
  # I chose two levels for t_d = 6,14, and 2 levels for U_d = 0.4,2
  
  .loadSim(sims[1])

  # Pull parameters
  Bmsy  <- blob$opMod$Bmsy
  nT    <- max( blob$opMod$lYear - blob$opMod$fYear + 1 )
  nReps <- blob$ctrl$signifReps
  nS    <- blob$opMod$nS


  # create an array for holding biomass
  depBlob <- array (  NA, dim = c(length(sims),nTraj,nS,nT),
                      dimnames = list(sims,1:nTraj,blob$ctrl$speciesNames[1:nS],1:nT) )

  Fblob   <- array (  NA, dim = c(length(sims),nS,nT),
                      dimnames = list(sims,blob$ctrl$speciesNames[1:nS],1:nT) ) 

  # Also collect Fhist info
  initDep   <- blob$opMod$initDep
  Umax      <- numeric( length(sims) )

  # set random seed, select traces
  set.seed(seed)

  traces <- sample( 1:nReps, nTraj )

  # Save first sim's biomass and Fhist info
  for( sIdx in 1:nS )
  {
    depBlob[ 1, , sIdx , ]   <- blob$om$Bt[ traces, sIdx, ] / Bmsy[sIdx] / 2
    Fblob[ 1, sIdx,  ]       <- blob$om$Ut[ 1, sIdx, ]
  }
  Umax[ 1 ]         <- blob$opMod$Umult[2]
  

  # Now loop over the other three and collect the info
  for( simIdx in 2:length(sims) )
  {
    .loadSim( sims[simIdx] )
    for( sIdx in 1:nS )
    {
      depBlob[ simIdx, , sIdx , ]   <- blob$om$Bt[ traces, sIdx, ] / Bmsy[sIdx] / 2
      Fblob[ simIdx, sIdx,  ]       <- blob$om$Ut[ 1, sIdx, ]
    }
    Umax[ simIdx ]       <- blob$opMod$Umult[2]
  }

  ordUmax <- unique( Umax ) [ order( unique( Umax  ) ) ]


  nrows <- length( ordUmax )
  ncols <- nS
  # calculate median depletion trajectory
  medDep <- apply ( X = depBlob, FUN = median, MARGIN = c(1,3,4), na.rm = T )

  # If saving to figs folder
  if( save )
  {
    outFile <- file.path(getwd(),"project","figs", fileName )
    pdf( file = outFile, ... )
  } 
  # Make plot labels
  # labs <- c("(a)", "(b)", "(c)", "(d)" )
  counter <- 1

  # now make x axis year ticks
  years <- seq(1984,2017,by=5)
  years <- c(years,2017)


  par( mfrow = c( nrows, ncols ), mar = c(1.5,1.5,1,1.5), oma = c(3,3,1,3) )
  # Loop over U and t values
  for( U in ordUmax )
  {
    for ( sIdx in 1:nS )
    { 
      sim <- which ( Umax == U )
      plot( x = c(1984,2017), y = c(0,1.2), type = "n", axes = FALSE,
            xlab = "", ylab = "" )
      if( length(sim) == 0 ) next
        for( traj in 1:nTraj ) 
          lines( x = 1984:2017, 
                 y = depBlob[ sim, traj, sIdx, ],
                 col = "grey80", lwd = 0.5 )
        lines( x = 1984:2017, y = medDep[ sim, sIdx, ], lwd = 2, lty = 2 )
        lines( x = 1984:2017, y = Fblob[ sim, sIdx, ], lty = 1, lwd = 2)
        axis( side = 1, at = years )
        axis( side = 2, las = 1 )
        if( sIdx < nS ) axis( side = 4, las = 1, labels = FALSE )
        else axis( side = 4, las = 1 )
        # panLab( txt = labs[counter], x = 0.8, y = 0.9, cex = 1.2 )
      counter <- counter + 1
    }
  }
  mtext( side = 2, outer = TRUE, text = "Depletion", cex = 1.5, line = 1.5 )
  mtext( side = 1, outer = TRUE, text = "Year", cex = 1.5, line = 1.5 )
  mtext( side = 4, outer = TRUE, text = "Harvest Rate", cex = 1.5, line = 1.5 )
  if( save ) dev.off()
}


plotUPriorSens <- function (  tableName = "priorSens_pub_MRE",
                              nSp = 2,
                              betaU = c(.01,.3),
                              showLegend = FALSE )
{
  # plotUPriorSens()
  # Reads in a stats table for prior sensitivity analyses
  # and plots the performance
  # inputs:     tableName = character name of stats table
  #             par=paremeter that is affected by sens analysis
  # outputs:    NULL
  # side-effs:  plots to quartz

  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )
  
  # restrict to the table with the 
  table     <-  table %>%
                # filter( fixqbar == TRUE ) %>%
                filter( nS == nSp ) %>%
                group_by( mp, species) %>%
                summarise(  msUmsy = mean(msUmsy),
                            msUmsy025 = mean(msUmsy025),
                            msUmsy975 = mean(msUmsy975),
                            UmsyOM = mean(UmsyOM),
                            sigU2 = mean(sigU2),
                            sigU2P2 = mean(sigU2P2),
                            Umsybar = mean(Umsybar),
                            Umsybar025 = mean(Ubar025),
                            Umsybar975 = mean(Ubar975),
                            sigU2025 = mean(sigU2025),
                            sigU2975 = mean(sigU2975),
                            ssUmsy = mean(ssUmsy),
                            ssUmsy025 = mean(ssUmsy025),
                            ssUmsy975 = mean(ssUmsy975),
                            s2U       = mean(s2Umsy) ) %>%
                mutate( msUmsy      = UmsyOM * msUmsy + UmsyOM,
                        msUmsy025   = UmsyOM * msUmsy025 + UmsyOM,
                        msUmsy975   = UmsyOM * msUmsy975 + UmsyOM,
                        ssUmsy      = UmsyOM * ssUmsy + UmsyOM,
                        ssUmsy025   = UmsyOM * ssUmsy025 + UmsyOM,
                        ssUmsy975   = UmsyOM * ssUmsy975 + UmsyOM ) %>%
                dplyr::select(  msUmsy, msUmsy025, msUmsy975, 
                                ssUmsy, ssUmsy025, ssUmsy975, 
                                Umsybar, Umsybar025, Umsybar975,
                                mp, species, UmsyOM, 
                                sigU2, sigU2025, sigU2975, s2U,
                                sigU2P2 )

  # plotting quantities
  s2U <- unique(table$s2U)
  s2U <- s2U[ order(s2U) ]

  species     <- unique( table$species )
  nS          <- length( species )

  # Species plotting locations
  specX       <- seq (-0.3,0.3, length=nS )
  fullX       <- c()
  for( k in 1:length(s2U) ) fullX <- c( fullX, k + specX ) 

  # pch for each species
  pchSpec <- c(0,1,2,5,6)
  
  textLab <- c("(a)","(b)","(c)","(d)")

  par( mfrow = c(length(betaU),1), mar = c(1,1,2,1), oma = c(3,3,0,0) )
  for( bIdx in 1:length(betaU) )
  {
    beta <- betaU[bIdx]
    # plotMat <- matrix(NA, nrow=nS*length(s2U), ncol = 7 )
    plotTable <-  table %>% filter( sigU2P2 == beta ) 
    # browser()
    Umsybar <- plotTable %>% group_by(s2U) %>% 
            dplyr::select(s2U,Umsybar,Umsybar025,Umsybar975)

    sigU2 <- plotTable %>% group_by(s2U) %>% 
            dplyr::select(s2U,sigU2,sigU2025,sigU2975)

    Uspec <- plotTable %>% group_by(s2U) %>% 
            dplyr::select(s2U,msUmsy,msUmsy025,msUmsy975,UmsyOM,ssUmsy,ssUmsy025,ssUmsy975 )
    
    plot( x = c( 0, length(s2U) + 1 ), y = c( 0.9*min(plotTable$msUmsy025, plotTable$UmsyOM), 1.1*max(plotTable$msUmsy975, plotTable$UmsyOM) ), 
          type = "n",
          axes = FALSE, xlab = "", ylab = "" )
      polygon(  x = c( fullX, rev(fullX) ),
                y = c(Umsybar$Umsybar025, rev(Umsybar$Umsybar975) ),
                col = "grey90", border = NA )
      lines( x = fullX, y = Umsybar$Umsybar, lty = 2, lwd = 2 )
      abline( h = exp(-1.23), lty = 3, lwd = 2 )
      points( x = fullX, y = Uspec$msUmsy, pch = pchSpec[1:nS], cex = 1.5 )
      points( x = fullX+0.05, y = Uspec$ssUmsy, pch = pchSpec[1:nS], cex = 1.5 )
      points( x = fullX, y = Uspec$UmsyOM, pch = 20, cex = 1.5 )
      segments( x0 = fullX, y0 = as.numeric( Uspec$msUmsy025 ),
                y1 = as.numeric( Uspec$msUmsy975 ), 
                col = "grey40", lwd = 2 )
      segments( x0 = fullX+0.05, y0 = as.numeric( Uspec$ssUmsy025 ),
                y1 = as.numeric( Uspec$ssUmsy975 ), 
                col = "grey40", lwd = 2, lty = 4 )
      panLab( x = 0.05, y = 0.8, txt = textLab[bIdx] )
      if( showLegend & bIdx == 1)
        panLegend(  x = 0.05, y = 0.1, legTxt = c( "q OM", "q Est", "qbar"),
                    bty = "n", pch = c( 1, 16, NA ), 
                    lty = c( NA, NA, 2 ),
                    lwd = c( NA, NA, 2 ),
                    cex = c( 1.5, 1.5, NA ) )
      axis( side = 1, at = 1:length(s2U), labels = s2U )
      axis( side = 2, las = 1 )
      if( bIdx > 1 ) axis( side = 3, at = 1:length(s2U), labels = FALSE )
  }

  mtext ( side = 2, text = "Optimal exploitation rate (Umsy)", las = 0 , outer = TRUE,
          line = 2 )

  mtext ( side = 1, text = "Mean Hyperprior Variance (v_U)", 
          las = 0 , outer = TRUE,
          line = 2 )
}


plotqPriorSens <- function (  tableName = "priorSens_pub_MRE",
                              nSp = 2,
                              betaq = c(.01,.3),
                              showLegend = FALSE )
{
  # plotqPriorSens()
  # Reads in a stats table for prior sensitivity analyses
  # and plots the performance
  # inputs:     tableName = character name of stats table
  #             par=paremeter that is affected by sens analysis
  # outputs:    NULL
  # side-effs:  plots to quartz

  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )

  DERPAcols <- brewer.pal(n = 5, name = "Set1" )

  # restrict to the table with the 
  table <-      table %>%
                # filter( fixqbar == TRUE ) %>%
                filter( nS == nSp, tauq2P2 %in% betaq ) %>%
                group_by( mp, species) %>%
                summarise(  msq = mean(msq),
                            qOM = mean(qOM),
                            tauq2 = mean(tauq2),
                            qbar = mean(qbar),
                            msq025 = mean(msq025),
                            msq975 = mean(msq975),
                            qbar025 = mean(qbar025),
                            qbar975 = mean(qbar975),
                            tauq2025 = mean(tauq2025),
                            tauq2975 = mean(tauq2975),
                            tauq2P2 = mean(tauq2P2),
                            s2q     = mean(s2q),
                            ssq     = mean(ssq),
                            ssq025  = mean(ssq025),
                            ssq975  = mean(ssq975) ) %>%
                mutate( msq     = qOM * msq + qOM,
                        msq025  = qOM * msq025 + qOM,
                        msq975  = qOM * msq975 + qOM,
                        ssq     = qOM * ssq + qOM,
                        ssq025  = qOM * ssq025 + qOM,
                        ssq975  = qOM * ssq975 + qOM ) %>%
                dplyr::select(  msq, msq025, msq975, 
                                qbar, qbar025, qbar975,
                                mp, species, qOM, s2q, tauq2P2,
                                tauq2, tauq2025, tauq2975, species,
                                ssq025, ssq975, ssq )
  # plotting quantities
  s2q <- unique(table$s2q)
  s2q <- s2q[ order(s2q) ]

  species     <- unique( table$species )
  nS          <- length( species )

  # Species plotting locations
  specX       <- seq (-0.3,0.3, length=nS )
  fullX       <- c()
  for( k in 1:length(s2q) ) fullX <- c( fullX, k + specX ) 

  # pch for each species
  pchSpec <- 1:5
  
  textLab <- c("(a)","(b)","(c)","(d)")

  par( mfrow = c(length(betaq),1), mar = c(1,1,2,1), oma = c(3,3,0,0) )
  for( bIdx in 1:length(betaq) )
  {
    beta <- betaq[bIdx]
    # plotMat <- matrix(NA, nrow=nS*length(s2q), ncol = 7 )
    plotTable <-  table %>% filter( tauq2P2 == beta ) 
    # browser()
    qbar <- plotTable %>% group_by(s2q) %>% 
            dplyr::select(s2q,qbar,qbar025,qbar975)

    tauq2 <- plotTable %>% group_by(s2q) %>% 
            dplyr::select(s2q,tauq2,tauq2025,tauq2975)

    qSpec <- plotTable %>% group_by(s2q) %>% 
            dplyr::select(s2q,msq,msq025,msq975,qOM,ssq,ssq025,ssq975 )
    
    plot( x = c( 0, length(s2q) + 1 ), y = c( 0.9*min(plotTable$msq025, plotTable$qOM), 1.1*max(plotTable$msq975, plotTable$qOM) ), 
          type = "n",
          axes = FALSE, xlab = "", ylab = "" )
      polygon(  x = c( fullX, rev(fullX) ),
                y = c(qbar$qbar025, rev(qbar$qbar975) ),
                col = "grey90", border = NA )
      lines( x = fullX, y = qbar$qbar, lty = 2, lwd = 2 )
      abline( h = exp(-0.51), lty = 3, lwd = 2 )
      points( x = fullX, y = qSpec$msq, pch = pchSpec[1:nS], cex = 1.5,
              col = DERPAcols )
      points( x = fullX+0.05, y = qSpec$ssq, pch = pchSpec[1:nS], cex = 1.5,
              col = DERPAcols )
      points( x = fullX, y = qSpec$qOM, pch = 20, cex = 1.5, col = DERPAcols )
      segments( x0 = fullX, y0 = as.numeric( qSpec$msq025 ),
                y1 = as.numeric( qSpec$msq975 ), 
                col = DERPAcols, lwd = 2 )
      segments( x0 = fullX+0.05, y0 = as.numeric( qSpec$ssq025 ),
                y1 = as.numeric( qSpec$ssq975 ), 
                col = DERPAcols, lwd = 2, lty = 4 )
      panLab( x = 0.05, y = 0.8, txt = textLab[bIdx] )
      if( showLegend & bIdx == 1)
        panLegend(  x = 0.05, y = 0.1, legTxt = c( "q OM", "q Est", "qbar"),
                    bty = "n", pch = c( 1, 16, NA ), 
                    lty = c( NA, NA, 2 ),
                    lwd = c( NA, NA, 2 ),
                    cex = c( 1.5, 1.5, NA ) )
      axis( side = 1, at = 1:length(s2q), labels = s2q )
      axis( side = 2, las = 1 )
      if( bIdx > 1 ) axis( side = 3, at = 1:length(s2q), labels = FALSE )
  }

  mtext ( side = 2, text = "Catchability", las = 0 , outer = TRUE,
          line = 2 )

  mtext ( side = 1, text = "Mean Hyperprior Variance (v_q)", 
          las = 0 , outer = TRUE,
          line = 2 )
}


plotObsRasters <- function  ( tableName = "obsErr_pub_MARE",
                              pars = c("BnT","Umsy","q","Dep"),
                              wdh = 18, hgt = 18,
                              breakRange = c(-1,1),
                              compMean = TRUE,
                              hess = FALSE,
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
                        "obsErrDelta.pdf", sep = "")
      outFile <- file.path(mpFolder,outFile)
      pdf ( file = outFile, width = wdh, height = hgt )
      plotObsCompContour( tableName = fileName,
                          mpLabel = mp,
                          pars = pars,
                          nSp = grid[g,"nS"],
                          tcex = 1.5,
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
                          model = "ms", tcex = 1.5,
                          breakRange = breakRange,
                          compMean = compMean,
                          hess = hess )
      dev.off()
      outFile <- paste( grid$nS[g], mp,
                        "obsErr_SS.pdf", sep = "")
      outFile <- file.path( mpFolder, outFile )
      pdf ( file = outFile, width = wdh, height = hgt )
      plotObsModContour(  tableName = fileName,
                          mpLabel = mp,
                          pars = pars,
                          nSp = grid[g,"nS"],
                          model = "ss", tcex =1.5,
                          breakRange = breakRange,
                          compMean = compMean,
                          hess = hess )
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
                                  model       = "ss",
                                  devLabels   = FALSE,
                                  hess        = FALSE 
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

  if( hess ) pars <- c( pars, "HessPD" )
 
  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    z <- array( NA, dim=c(length(CVs),length(CVs),length(pars)),
                dimnames=list(CVs[length(CVs):1],CVs,pars) )

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
        colNames <- paste(model,pars,sep="")
        z[length(CVs)-yIdx+1,xIdx,] <- as.numeric(specTab[rows,colNames])
      }
    } 
    # browser()
    specList[[s]] <- brick( z, xmn = 0, xmx = length(CVs), ymn = 0, ymx=length(CVs) )
  }

  names(specList) <- species

  par(mfrow = c((nS+1),length(pars) ),  mar = c(1,1,1,1), oma = c(5,5,2,2) )
  for (s in 1:nS)
  {
    # Create a raster from the stat table info
    plotBrick <- specList[[s]]
    # Make titles for the plots
    specName <- species[s]

    titles <- names(plotBrick)
    l <- length(titles)
    if( !hess ) l <- l + 1

    colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
    colBreaksHess <- seq(0,200,length = 8)
    colScale  <- brewer.pal(8, "Reds")
    colScale      <- c(colScale[8:1],colScale)
    
    # plot the rasters, if first species plot main titles
    if ( s == 1 )
    {
      for(k in 1:(l-1))
      {
        plotBrick[[ k ]][plotBrick[[k]] > 10 ] <- 99
        plot( plotBrick[[ k ]],main = titles[ k ], legend = FALSE, 
              asp=NA, ylab = "", xlab = "", axes = FALSE, 
              col = colScale, breaks = colBreaks )
          if( k == 1 ) mtext ( side = 2, text = specName, line = 2, las = 0)
          text( plotBrick[[ k ]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
          axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
      }
      if( hess )
      {
        plot( plotBrick[[ l ]],main = titles[ l ], legend = FALSE, 
              asp=NA, ylab = "", xlab = "", axes = FALSE, 
              col = colScale[9:16], breaks = colBreaksHess )
          text( plotBrick[[ l ]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
          axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )  
      }
    }
    else {
      # browser()
      for(k in 1:( l - 1 ) )
      {
        plotBrick[[ k ]][plotBrick[[k]] > 10 ] <- 99
        plot( plotBrick[[k]],main = "", legend = FALSE, 
              asp=NA, ylab = "", xlab = "", axes = FALSE, 
              col = colScale, breaks = colBreaks )
          text( plotBrick[[k]], digits = 2, halo = T, cex = tcex ) 
          if( k == 1 ) mtext ( side = 2, text = specName, line = 2, las = 0)
          axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
          axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
      }
      if( hess )
      {
        plot( plotBrick[[l]],main = "", legend = FALSE, 
                    asp=NA, ylab = "", xlab = "", axes = FALSE, 
                    col = colScale[9:16], breaks = colBreaksHess )
          text( plotBrick[[l]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
          axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
      }
    }   
  }

  # browser()
  # Now create an average/sum raster brick for the complex
  # fudge factor for hess logical flag
  m <- l
  if( !hess ) m <- l-1
  for(k in 1:m)
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
  
  for(k in 1:m )
  {
    plotBrick[[ k ]][plotBrick[[k]] > 10 ] <- 99
    plot( plotBrick[[k]],main = "", legend = FALSE, 
          asp=NA, ylab="", xlab = "", axes = FALSE, 
          col = colScale, breaks = colBreaks )
      if( k == 1 ) mtext ( side = 2, text = "Complex", line = 2, las = 0)
      text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )
      axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
      axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )   
  }
  if( hess )
  {
    plot( plotBrick[[l]],main = "", legend = FALSE, 
          asp=NA, ylab="", xlab = "", axes = FALSE,
          col = colScale[9:16], breaks = nSp*colBreaksHess )
      text( plotBrick[[l]], digits = 2, halo = T, cex = tcex )   
      axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
      axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )   
  }
  

  mtext( side = 1, text = "Dover Sole Survey CV", outer = TRUE, las = 0,
          line = 3 )
  mtext( side = 2, text = "English Sole Survey CV", outer = TRUE, las = 0,
          line = 3 )
  if( devLabels )
  {
    grid.text(  paste("Raw error ", model, " model.", sep = ""), 
                x=unit(0.5, "npc"), y=unit(0.98, "npc"), rot=0)
    grid.text(  tableName,
                x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
    grid.text(  mpLabel,
                x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)  
  }
  
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
                                  pars        = c("BnT", "Umsy", "q", "Dep"),
                                  nSp         = 2,
                                  method      = "deviance",
                                  compMean    = FALSE,                                
                                  tcex        = 2,
                                  breakRange  = c(-2,2),
                                  mpLabel     = "corrTRUE_SigPriorIG",
                                  legend      = FALSE,
                                  devLabels   = FALSE 
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
                                Dep     = log2(abs(ssDep/msDep)),
                                Bmsy    = log2(abs(ssBmsy/msBmsy)) )
  }   
  # Difference is just ss - ms
  if ( method == "difference")
  {
    table <- table %>% mutate(  BnT     = abs(ssBnT) - abs(msBnT),
                                Umsy    = abs(ssUmsy) - abs(msUmsy),
                                q       = abs(ssq) - abs(msq),
                                Dep     = abs(ssDep) - abs(msDep),
                                Bmsy    = log2(abs(ssBmsy/msBmsy)) )
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
    specList[[s]] <- brick( z, xmn = 0, xmx = length(CVs), ymn = 0, ymx=length(CVs) )
  }
  names(specList) <- species

  # Start plotting

  par(mfrow = c(nS+1,length(pars)), mar = c(1.5,1.5,1.5,1.5), oma = c(5,5,2,2) )
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
      plotBrick[[ k ]][plotBrick[[k]] > 100 ] <- 99
      plotBrick[[ k ]][plotBrick[[k]] < -100 ] <- -99 
      if(s == 1) plot( plotBrick[[k]],main = titles[k], legend = FALSE, 
                        breaks = colBreaks, col=colScale, asp=NA, axes = FALSE )
      else plot( plotBrick[[k]],main = "", legend = FALSE, axes = FALSE, 
                  breaks = colBreaks, col=colScale, asp=NA )
      if( k == 1 ) mtext ( text = species[s], side = 2, las = 0, line = 2 )
      axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
      axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
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
                                      ssBmsy = sum( ssBmsy ),
                                      msBmsy = sum( msBmsy ),
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
                          Bmsy = log2( abs(ssBmsy) / abs(msBmsy) ),
                          Dep = log2 ( abs(ssDep) / abs(msDep) ) 
                        )
  }

  if( method == "difference" )
  {
    compTable <-  compTable %>% 
                  mutate( BnT = abs(ssBnT) - abs(msBnT) ,
                          Umsy = abs(ssUmsy) - abs(msUmsy) ,
                          q = abs(ssq) - abs(msq) ,
                          Bmsy = abs(ssBmsy) - abs(msBmsy) ,
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
      rows <- which ( compTable[,"scenario"] == scenName )
      if( length(rows) == 0 ) next
      # Recover MSE comparison
      z[length(CVs) - yIdx + 1,xIdx,] <- as.numeric(compTable[rows,pars])
    }
  } 
  compBrick <- brick( z, xmn = 0, xmx = length(CVs), ymn = 0, ymx=length(CVs) )

  for(k in 1:length(titles))
  {
    compBrick[[ k ]][compBrick[[k]] > 10 ] <- 99
    compBrick[[ k ]][compBrick[[k]] < -10 ] <- -99
    plot( compBrick[[k]],main = "", legend = FALSE, 
          asp=NA, col = colScale, breaks = colBreaks,
          axes = FALSE )
      if( k == 1 ) mtext ( text = "Complex", side = 2, las = 0, line = 2 )
    text( compBrick[[k]], digits = 2, halo = T, cex = tcex )
    axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
    axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
  }

  if( legend )
  {
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
  
  }
  
  # Labels
  mtext( side = 1, text = "Dover Sole Survey CV", outer = TRUE, las = 0,
          line = 3 )
  mtext( side = 2, text = "English Sole Survey CV", outer = TRUE, las = 0,
          line = 3 )
  if( devLabels )
  {
    grid.text(  "Relative errors ss/ms", 
                x=unit(0.5, "npc"), y=unit(0.95, "npc"), rot=0)
    grid.text(  tableName,
                x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
    grid.text(  mpLabel,
                x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)
  }
  # grid.text(  paste( "lastNegCorr = ", negCorr, sep = "" ),
  #             x=unit(0.15, "npc"), y=unit(0.02, "npc"), rot=0)
  # grid.text(  paste( "tUtrough = ", tUtr, sep = "" ),
  #             x=unit(0.15, "npc"), y=unit(0.05, "npc"), rot=0)

}

plotComplexPubRasters <- function(  tableRoot = "allSame_RE_msIncr_",
                                    axes = c("corr","kappaMult"),
                                    pars = c("BnT","Umsy","q","Dep"),
                                    wdh = 14, hgt = 18,
                                    breakRange = c(-1,1),
                                    tcex = 1.5,
                                    legend = FALSE,
                                    mpLabel   = "allJointPriors",
                                    hess = FALSE,
                                    nSp       = 5,
                                    tUtr      = 15,
                                    negCorr   = FALSE,
                                    devLabels = FALSE,
                                    save      = FALSE
                                  )
{
  # Load in MARE and MRE tables
  # browser()
  MAREfile  <- paste( tableRoot, "_MARE.csv", sep = "" )
  MREfile   <- paste( tableRoot, "_MRE.csv", sep = "" )

  MAREpath  <- file.path ( getwd(), "project/Statistics", MAREfile )
  MREpath   <- file.path ( getwd(), "project/Statistics", MREfile )

  MREtable  <- read.csv( MREpath, header = TRUE, stringsAsFactors = FALSE )
  MAREtable <- read.csv( MAREpath, header = TRUE, stringsAsFactors = FALSE )

  # reduce to the correct subset of sims
  MREtable <- MREtable  %>% 
              filter( mp == mpLabel,
                      nS == nSp,
                      tUtrough == tUtr,
                      lastNegCorr == negCorr ) %>%
              group_by( scenario, mp ) %>%
              summarise( ssBnT = mean( ssBnT ),
                         msBnT = mean( msBnT ),
                         ssUmsy = mean( ssUmsy ),
                         msUmsy = mean( msUmsy ),
                         ssq = mean( ssq ),
                         msq = mean( msq ),
                         ssDep = mean( ssDep ),
                         msDep = mean( msDep ),
                         corr = mean(corr),
                         kappaMult = mean(kappaMult) )
              
  
  MAREtable <-  MAREtable  %>% 
                filter( mp == mpLabel,
                        nS == nSp,
                        tUtrough == tUtr,
                        lastNegCorr == negCorr ) %>%
                group_by( scenario, mp ) %>%
                summarise(  ssBnT = mean( ssBnT ),
                            msBnT = mean( msBnT ),
                            ssUmsy = mean( ssUmsy ),
                            msUmsy = mean( msUmsy ),
                            ssq = mean( ssq ),
                            msq = mean( msq ),
                            ssDep = mean( ssDep ),
                            msDep = mean( msDep ),
                            corr = mean(corr),
                            kappaMult = mean(kappaMult) )
  # Compute Delta values              
  MAREtable <-  MAREtable %>%
                mutate( BnT     = log2(abs(ssBnT/msBnT)),
                        Umsy    = log2(abs(ssUmsy/msUmsy)),
                        q       = log2(abs(ssq/msq)),
                        Dep     = log2(abs(ssDep/msDep)) )

  # Now start pulling rasters
  msColNames      <- paste( "ms", pars, sep = "" )
  ssColNames      <- paste( "ss", pars, sep = "" )
  DeltaColNames   <- pars

  # use makeRasterFromTable() to create raster bricks
  ssMRE <- makeRasterFromTable( MREtable, axes, ssColNames )
  msMRE <- makeRasterFromTable( MREtable, axes, msColNames )
  ssMARE <- makeRasterFromTable( MAREtable, axes, ssColNames )
  msMARE <- makeRasterFromTable( MAREtable, axes, msColNames )
  Delta <- makeRasterFromTable( MAREtable, axes, DeltaColNames )

  if( save )
  {
    outFile <- paste( tableRoot, "complexPubRasters.pdf", sep = "")
    outFile <- file.path("./project/figs/",outFile)
    pdf ( file = outFile, width = wdh, height = hgt )
  }
  
  par(mfrow = c(5,length(pars)), mar = c(1,1,1,1), oma = c(5,5,3,2) ) 
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colBreaksHess <- seq(0,200,length = 8)
  redScale    <- brewer.pal(8, "Reds")
  greenScale  <- brewer.pal(8, "Greens" )
  rgScale     <- c(redScale[8:1],greenScale)
  rrScale     <- c(redScale[8:1],redScale)

  titles <- c("B2017","Umsy","q","D2017")

  # First plot Delta
  for(k in 1:length(pars) )
  {
    plot( Delta[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rgScale, breaks = colBreaks,
        axes = FALSE, cex.main = 1.5 )
      mtext( text = titles[k], side = 3, line = 2, cex = 1.5)
      if( k == 1 ) mtext ( text = "(a)", side = 2, las = 1, line = 2 )
      text( Delta[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(Delta$y) - 0.5, labels = ssMARE$y, las = 1 )
      axis( side = 1, at = 1:length(Delta$x) - 0.5, labels = Delta$x ) 
  }  
  # Then MARE plots
  for(k in 1:length(pars) )
  {
    plot( ssMARE[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
        axes = FALSE )
      if( k == 1 ) mtext ( text = "(b)", side = 2, las = 1, line = 2 )
      text( ssMARE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(ssMARE$y) - 0.5, labels = ssMARE$y, las = 1 )
      axis( side = 1, at = 1:length(ssMARE$x) - 0.5, labels = ssMARE$x ) 
  }  
  for(k in 1:length(pars) )
  {
    plot( msMARE[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
        axes = FALSE )
      if( k == 1 ) mtext ( text = "(c)", side = 2, las = 1, line = 2 )
      text( msMARE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(msMARE$y) - 0.5, labels = msMARE$y, las = 1 )
      axis( side = 1, at = 1:length(msMARE$x) - 0.5, labels = msMARE$x ) 
  }  
  # Then MRE plots
  for(k in 1:length(pars) )
  {
    plot( ssMRE[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
        axes = FALSE )
      if( k == 1 ) mtext ( text = "(d)", side = 2, las = 1, line = 2 )
      text( ssMRE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(ssMRE$y) - 0.5, labels = ssMRE$y, las = 1 )
      axis( side = 1, at = 1:length(ssMRE$x) - 0.5, labels = ssMRE$x ) 
  }  
  for(k in 1:length(pars) )
  {
    plot( msMRE[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
        axes = FALSE )
      if( k == 1 ) mtext ( text = "(e)", side = 2, las = 1, line = 2 )
      text( msMRE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(msMRE$y) - 0.5, labels = msMRE$y, las = 1 )
      axis( side = 1, at = 1:length(msMRE$x) - 0.5, labels = msMRE$x ) 
  }  
  mtext(  side = 1, outer = TRUE, text = "Correlation in species effects",
          line = 2 )
  mtext(  side = 2, outer = TRUE, text = "Relative magnitude of shared effects",
          line = 3 )
  if( save ) dev.off()
}


plotComplexPubRasters2 <- function(  tableRoot = "RE_coarse_pub",
                                    axes = c("corr","kappaMult"),
                                    pars = c("BnT","Bmsy","Umsy","Dep"),
                                    wdh = 14, hgt = 18,
                                    breakRange = c(-1,1),
                                    tcex = 1.5,
                                    legend = FALSE,
                                    mpLabel   = "pubAM",
                                    hess = FALSE,
                                    nSp       = 5,
                                    tUtr      = 15,
                                    negCorr   = FALSE,
                                    devLabels = FALSE,
                                    save      = FALSE
                                  )
{
  # Load in MARE and MRE tables
  # browser()
  MAREfile  <- paste( tableRoot, "_MARE.csv", sep = "" )
  MREfile   <- paste( tableRoot, "_MRE.csv", sep = "" )

  MAREpath  <- file.path ( getwd(), "project/Statistics", MAREfile )
  MREpath   <- file.path ( getwd(), "project/Statistics", MREfile )

  MREtable  <- read.csv( MREpath, header = TRUE, stringsAsFactors = FALSE )
  MAREtable <- read.csv( MAREpath, header = TRUE, stringsAsFactors = FALSE )

  # reduce to the correct subset of sims
  MREtable <- MREtable  %>% 
              filter( mp == mpLabel,
                      nS == nSp,
                      tUtrough == tUtr,
                      lastNegCorr == negCorr ) %>%
              group_by( scenario, mp ) %>%
              summarise( ssBnT = mean( ssBnT ),
                         msBnT = mean( msBnT ),
                         ssBmsy = mean( ssBmsy ),
                         msBmsy = mean( msBmsy ),
                         ssUmsy = mean( ssUmsy ),
                         msUmsy = mean( msUmsy ),
                         ssq = mean( ssq ),
                         msq = mean( msq ),
                         ssDep = mean( ssDep ),
                         msDep = mean( msDep ),
                         corr = mean(corr),
                         kappaMult = mean(kappaMult) )
              
  
  MAREtable <-  MAREtable  %>% 
                filter( mp == mpLabel,
                        nS == nSp,
                        tUtrough == tUtr,
                        lastNegCorr == negCorr ) %>%
                group_by( scenario, mp ) %>%
                summarise(  ssBnT = mean( ssBnT ),
                            msBnT = mean( msBnT ),
                            ssBmsy = mean( ssBmsy ),
                            msBmsy = mean( msBmsy ),
                            ssUmsy = mean( ssUmsy ),
                            msUmsy = mean( msUmsy ),
                            ssq = mean( ssq ),
                            msq = mean( msq ),
                            ssDep = mean( ssDep ),
                            msDep = mean( msDep ),
                            corr = mean(corr),
                            kappaMult = mean(kappaMult) )
  # Compute Delta values              
  MAREtable <-  MAREtable %>%
                mutate( BnT     = log2(abs(ssBnT/msBnT)),
                        Umsy    = log2(abs(ssUmsy/msUmsy)),
                        Bmsy    = log2(abs(ssBmsy/msBmsy)),
                        q       = log2(abs(ssq/msq)),
                        Dep     = log2(abs(ssDep/msDep)) )

  # Now start pulling rasters
  msColNames      <- paste( "ms", pars, sep = "" )
  ssColNames      <- paste( "ss", pars, sep = "" )
  DeltaColNames   <- pars

  # use makeRasterFromTable() to create raster bricks
  ssMRE <- makeRasterFromTable( MREtable, axes, ssColNames )
  msMRE <- makeRasterFromTable( MREtable, axes, msColNames )
  ssMARE <- makeRasterFromTable( MAREtable, axes, ssColNames )
  msMARE <- makeRasterFromTable( MAREtable, axes, msColNames )
  Delta <- makeRasterFromTable( MAREtable, axes, DeltaColNames )

  if( save )
  {
    outFile <- paste( tableRoot, "complexPubRasters.pdf", sep = "")
    outFile <- file.path("./project/figs/",outFile)
    pdf ( file = outFile, width = wdh, height = hgt )
  }
  
  par(mfrow = c(1,length(pars)), mar = c(1,1,1,1), oma = c(5,5,3,2) ) 
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colBreaksHess <- seq(0,200,length = 8)
  redScale    <- brewer.pal(8, "Reds")
  greenScale  <- brewer.pal(8, "Greens" )
  rgScale     <- c(redScale[8:1],greenScale)
  rrScale     <- c(redScale[8:1],redScale)

  titles <- c("B2017","B0","Umsy","D2017")

  # First plot Delta
  for(k in 1:length(pars) )
  {
    plot( Delta[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rgScale, breaks = colBreaks,
        axes = FALSE, cex.main = 1.5 )
      mtext( text = titles[k], side = 3, line = 2, cex = 1.5)
      # if( k == 1 ) mtext ( text = "(a)", side = 2, las = 1, line = 2 )
      text( Delta[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(Delta$y) - 0.5, labels = ssMARE$y, las = 1 )
      axis( side = 1, at = 1:length(Delta$x) - 0.5, labels = Delta$x ) 
  }  
  # # Then MARE plots
  # for(k in 1:length(pars) )
  # {
  #   plot( ssMARE[[1]][[k]], legend = FALSE, 
  #       asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
  #       axes = FALSE )
  #     if( k == 1 ) mtext ( text = "(b)", side = 2, las = 1, line = 2 )
  #     text( ssMARE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
  #     axis( side = 2, at = 1:length(ssMARE$y) - 0.5, labels = ssMARE$y, las = 1 )
  #     axis( side = 1, at = 1:length(ssMARE$x) - 0.5, labels = ssMARE$x ) 
  # }  
  # for(k in 1:length(pars) )
  # {
  #   plot( msMARE[[1]][[k]], legend = FALSE, 
  #       asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
  #       axes = FALSE )
  #     if( k == 1 ) mtext ( text = "(c)", side = 2, las = 1, line = 2 )
  #     text( msMARE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
  #     axis( side = 2, at = 1:length(msMARE$y) - 0.5, labels = msMARE$y, las = 1 )
  #     axis( side = 1, at = 1:length(msMARE$x) - 0.5, labels = msMARE$x ) 
  # }  
  # # Then MRE plots
  # for(k in 1:length(pars) )
  # {
  #   plot( ssMRE[[1]][[k]], legend = FALSE, 
  #       asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
  #       axes = FALSE )
  #     if( k == 1 ) mtext ( text = "(d)", side = 2, las = 1, line = 2 )
  #     text( ssMRE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
  #     axis( side = 2, at = 1:length(ssMRE$y) - 0.5, labels = ssMRE$y, las = 1 )
  #     axis( side = 1, at = 1:length(ssMRE$x) - 0.5, labels = ssMRE$x ) 
  # }  
  # for(k in 1:length(pars) )
  # {
  #   plot( msMRE[[1]][[k]], legend = FALSE, 
  #       asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
  #       axes = FALSE )
  #     if( k == 1 ) mtext ( text = "(e)", side = 2, las = 1, line = 2 )
  #     text( msMRE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
  #     axis( side = 2, at = 1:length(msMRE$y) - 0.5, labels = msMRE$y, las = 1 )
  #     axis( side = 1, at = 1:length(msMRE$x) - 0.5, labels = msMRE$x ) 
  # }  
  mtext(  side = 1, outer = TRUE, text = "Inter-species Correlation",
          line = 2 )
  mtext(  side = 2, outer = TRUE, text = "Shared Environmental Effects",
          line = 3 )
  if( save ) dev.off()
}


pull <- function(x,y) 
{
  # pull()
  # Pulls a single column from a tbl_df object and returns
  # it as a vector. Usefule for working with dplyr tables
  # inputs:   x = tbl_df object
  #           y = column name/number
  x[,if(is.name(substitute(y))) deparse(substitute(y)) 
      else y, drop = FALSE][[1]]
}

# makeRasterFromTable()
# Creates a raster brick from a table, which is either a stat
# table created in .statTableXXX() or a pared down version 
# inside a function
# inputs:   table=data.frame of stats
#           axes=variables that make up the x and y axes of the raster
#           colnames=names of brick layers
# ouputs:   raster brick made up of table entries
makeRasterFromTable <- function ( table, 
                                  axes, 
                                  colnames )
{
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }

  x   <- unique(table %>% pull(axesCols[1]))
  x   <- x[order(x)]
  y   <- unique(table %>% pull(axesCols[2]))
  y   <- y[order(y)]


  # First, create an array to turn into the raster
  z <- array( NA, dim=c(length(y),length(x),length(colnames)),
                dimnames = list( y[ length( y ):1 ], x, colnames ) )
  for (xIdx in 1:length(x))
  {
    for (yIdx in 1:length(y))
    {
      # Get x and y values
      xVal <- x[xIdx]
      yVal <- y[yIdx]
      # Get rows
      xRows <- which ( table[,axesCols[1] ] == xVal)
      yRows <- which ( table[,axesCols[2] ] == yVal)
      rows <- intersect(xRows,yRows)
      if( length(rows) == 0 ) next
      z[length(y)-yIdx+1,xIdx,] <- as.numeric(table[rows,colnames])
    }
  }
  # create raster brick
  raster <- brick( z, xmn = 0, xmx = length(x), ymn = 0, ymx = length(y) )  
  # output raster brick and axis vectors
  out <- list()
  out$brick <- raster
  out$x <- x
  out$y <- y

  out
}

# plotRasters()
# Reads in a statistics table produced by .statTableXXX()
# and produces plots of performance contours.
# inputs:     tableName=charactre vector of file name root
#             axes=
plotRasters <- function ( tableName = "allSame_Fhist_msIncr_MARE",
                          axes = c("corr","kappaMult"),
                          pars = c("BnT","Umsy","q","Dep"),
                          wdh = 14, hgt = 18,
                          breakRange = c(-1,1),
                          compMean = TRUE,
                          tcex = 1.5,
                          legend = FALSE,
                          devLabels = FALSE,
                          hess = FALSE,
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
                        tcex = tcex,
                        breakRange = breakRange,
                        compMean = compMean,
                        legend = legend,
                        devLabels = devLabels,
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
                        model = "ms", tcex = tcex,
                        breakRange = breakRange,
                        compMean = compMean,
                        devLabels = devLabels,
                        hess = hess )
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
                        model = "ss", tcex = tcex,
                        breakRange = breakRange,
                        compMean = compMean,
                        devLabels = devLabels,
                        hess = hess )
      dev.off()
    }
  }
  # No return
}


# plotModContour()
# Plots the correlation  experiment MSE contours for
# estimation of BnT and Umsy (abundance and productivity) by species
# and for the complex.
plotModContour <- function (  tableName = "RE_coarse_Bmsy_MARE.csv", 
                              mpLabel   = "corrTRUE_SigPriorIG",
                              axes      = c("corr","kappaMult"),
                              pars      = c("BnT","Umsy","q"),
                              nSp       = 2,
                              model     = "ss",
                              compMean  = FALSE,
                              tUtr      = 15,
                              negCorr   = FALSE,
                              tcex      = 2,
                              breakRange= c(-1,1),
                              devLabels = FALSE,
                              hess      = FALSE )
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
  # Add hessianPD column if requested
  if ( hess ) pars <- c( pars, "HessPD" )

  x   <- unique(table[,axesCols[1] ])
  x   <- x[order(x)]
  y   <- unique(table[,axesCols[2] ])
  y   <- y[order(y)]

  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    z <- array( NA, dim=c(length(y),length(x),length(pars)),
                dimnames = list( y[ length( y ):1 ], x, pars ) )

    for (xIdx in 1:length(x))
    {
      for (yIdx in 1:length(y))
      {
        # Get x and y values
        xVal <- x[xIdx]
        yVal <- y[yIdx]
        # Get rows
        xRows <- which ( specTab[,axesCols[1] ] == xVal)
        yRows <- which ( specTab[,axesCols[2] ] == yVal)
        rows <- intersect(xRows,yRows)
        if( length(rows) == 0 ) next
        colNames <- paste(model,pars,sep="")
        z[length(y)-yIdx+1,xIdx,] <- as.numeric(specTab[rows,colNames])
      }
    } 
    # browser()
    specList[[s]] <- brick( z, xmn = 0, xmx = length(x), ymn = 0, ymx = length(y) )
  }

  names(specList) <- species
  
  par(mfrow = c((nS+1),length(pars)),  mar = c(1,1,1,1), oma = c(5,5,2,2) )
  for (s in 1:nS)
  {
    # Create a raster from the stat table info
    plotBrick <- specList[[s]]
    # Make titles for the plots
    specName <- species[s]

    titles <- names(plotBrick)
    if( hess ) l <- length(titles) else l <- length(titles) + 1

    colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
    colBreaksHess <- seq(0,200,length = 8)
    colScale  <- brewer.pal(8, "Reds")
    colScale      <- c(colScale[8:1],colScale)
    # plot the rasters, if first species plot main titles
    if ( s == 1 )
    {
      for(k in 1:(l-1))
      {
        plot( plotBrick[[k]], main = titles[k], legend = FALSE, 
              asp=NA, ylab="", xlab = "", col = colScale, breaks = colBreaks,
              axes = FALSE )
          if( k == 1 ) mtext ( text = species[s], side = 2, las = 0, line = 2 )
          text( plotBrick[[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
          axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
          axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
      }
      if( hess )
      {
        plot( plotBrick[[l]],main = titles[l], legend = FALSE, 
            asp=NA, ylab = "", xlab = "", col = colScale[9:16], 
            breaks = colBreaksHess,
            axes = FALSE )
          text( plotBrick[[l]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
          axis( side = 1, at = 1:length(x) - 0.5, labels = x )    
      }
      
    }
    else {
      for(k in 1:(l-1))
      {
        plot( plotBrick[[k]],main = "", legend = FALSE, 
              asp=NA, ylab = "", xlab = "", col = colScale, breaks = colBreaks,
              axes = FALSE )
          if( k == 1 ) mtext ( text = species[s], side = 2, las = 0, line = 2 )
          text( plotBrick[[k]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
          axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
      }
      if( hess )
      {
        plot( plotBrick[[l]],main = "", legend = FALSE, 
              asp=NA, ylab = "", xlab = "", col = colScale[9:16], 
              breaks = colBreaksHess, axes = FALSE )
          text( plotBrick[[l]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
          axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
      }
      
    }   
  }

  # Now create an average/sum raster brick for the complex
  # Some fudging due to the hessian logical flag
  if( hess ) m <- l else m <- l-1
  for( k in 1:m )
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
          asp=NA, ylab = "", xlab = "", col = colScale, breaks = colBreaks,
          axes = FALSE )
      if( k == 1 ) mtext ( text = "Complex", side = 2, las = 0, line = 2 )
      text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )  
      axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
      axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
  }
  if ( hess )
  {
    plot( plotBrick[[l]],main = "", legend = FALSE, 
          asp=NA, ylab = "", xlab = "", col = colScale, breaks = nSp * colBreaksHess )
      text( plotBrick[[l]], digits = 2, halo = T, cex = tcex )  
      axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
      axis( side = 1, at = 1:length(x) - 0.5, labels = x )   
  }
  
  # Labels
  mtext( side = 1, text = axes[1], outer = TRUE, cex = 2, line = 3, las = 0 )
  mtext( side = 2, text = axes[2], outer = TRUE, cex = 2, line = 3, las = 0)
  if( devLabels )
  {
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
plotCompContour <- function ( tableName   = "RE_coarse_pub_MARE.csv", 
                              mpLabel     = "pubAM",
                              axes        = c("corr","kappaMult"),
                              pars        = c("BnT","Umsy","q","Dep"),
                              nSp         = 2,
                              method      = "deviance",
                              compMean    = FALSE,
                              tUtr        = 15,
                              negCorr     = FALSE,
                              tcex        = 2,
                              breakRange  = c(-2,2),
                              devLabels   = FALSE,
                              legend      = FALSE  
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

  x   <- unique(table[,axesCols[1] ])
  x   <- x[order(x)]
  y   <- unique(table[,axesCols[2] ])
  y   <- y[order(y)]
  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    # browser()
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
        xRows <- which ( specTab[ ,axesCols[1] ] == xVal)
        yRows <- which ( specTab[ ,axesCols[2] ] == yVal)
        rows <- intersect(xRows,yRows)
        # Recover MSE comparison
        # browser()
        z[ length( y ) - yIdx + 1, xIdx, ] <- as.numeric(specTab[rows,pars])
      }
    } 
    specList[[s]] <- brick( z, xmn = 0, xmx = length(x), ymn = 0, ymx = length(y) )
  }
  names(specList) <- species

  # Start plotting

  par(mfrow = c(nS+1,length(pars)), mar = c(1,1,1,1), oma = c(5,5,2,2) )
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
      if(s == 1) plot( plotBrick[[k]], main = titles[k], legend = FALSE, 
                        breaks = colBreaks, col=colScale, asp=NA,
                        xlab = "", ylab = "", axes = FALSE )
      else plot( plotBrick[[k]],main = "", legend = FALSE, 
                  breaks = colBreaks, col=colScale, asp=NA,
                  xlab = "", ylab = "", axes = FALSE )
      if( k == 1 ) mtext ( text = species[s], side = 2, las = 0, line = 2 )
      text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )  
      axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
      axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
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
      xRows <- which ( compTable[,axesCols[1] ] == xVal)
      yRows <- which ( compTable[,axesCols[2] ] == yVal)
      rows <- intersect(xRows,yRows)
      # Recover MSE comparison
      z[length(y) - yIdx + 1,xIdx,] <- as.numeric(compTable[rows,pars])
    }
  } 
  compBrick <- brick( z, xmn = 0, xmx = length(x), ymn = 0, ymx = length(y) )

  for(k in 1:length(titles))
  {
    plot( compBrick[[k]],main = "", legend = FALSE, 
          asp=NA, col = colScale, breaks = colBreaks,
          xlab = "", ylab = "", axes = FALSE )
      if( k == 1 ) mtext ( text = "Complex", side = 2, las = 0, line = 2 )
    text( compBrick[[k]], digits = 2, halo = T, cex = tcex ) 
    axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
    axis( side = 1, at = 1:length(x) - 0.5, labels = x )   
  }

  # Plot the legends
  if ( legend )
  {
    plot( plotBrick[[1]], legend.only=TRUE, col=colScale, fill=colScale,
          legend.width=1, legend.shrink=0.75,
          axis.args = list( at = colBreaks,
                            labels = round(colBreaks,2),
                            cex.axis = 0.6),
          zlim = c(-3,3),
          legend.args = list( text="log2(ssMSE/msMSE)", 
                              side=4, font=2, 
                              line=3, cex=0.75, las=0) )
  
  }
  
  # Labels
  mtext( side = 1, text = axes[1], outer = TRUE, cex = 2, line = 3, las = 0 )
  mtext( side = 2, text = axes[2], outer = TRUE, cex = 2, line = 3, las = 0)
  if ( devLabels )
  {
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

plotBenvelopes <- function( sim = 1,
                              fYear = 1988,
                              est = FALSE,
                              labSize = 2, 
                              tickSize = 1.2,
                              devLabels = TRUE,
                              alpha = .8 )
{
  # First, load blob
  .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")

  # Recover blob elements for plotting
  nS <- blob$opMod$nS
  nT <- blob$opMod$nT

  # Species names
  specNames <- blob$ctrl$speciesNames

  # True OM quantities
  omBt  <- blob$om$Bt
  Bmsy  <- blob$opMod$Bmsy

  # Single species model
  ssBt  <- blob$am$ss$Bt

  # Multispecies model
  msBt  <- blob$am$ms$Bt
  msBt[ msBt < 0 ] <- NA

  years <- fYear:(fYear + nT - 1)

  for( s in 1:nS )
  {
    omBt[,s,]  <- omBt[,s,]/Bmsy[s]/2
    ssBt[,s,]  <- ssBt[,s,]/Bmsy[s]/2
    msBt[,s,]  <- msBt[,s,]/Bmsy[s]/2
  }  
  
  # Set colours
  ssCol <- "steelblue"
  msCol <- "salmon"
  # Create transparent fill colours for polys
  ssFill <- col2rgb(ssCol)/255
  ssFill <- rgb(ssFill[1],ssFill[2],ssFill[3],alpha=0.4)
  msFill <- col2rgb(msCol)/255
  msFill <- rgb(msFill[1],msFill[2],msFill[3],alpha=0.4)

  ## CHECK THAT THESE ARRAY CALCS ARE WORKING ##
  # Now, reduce to envelopes
  # OM values
  omBt  <- apply( X = omBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  # SS AM 
  ssBt  <- apply( X = ssBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  
  # MS AM 
  msBt  <- apply( X = msBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  
  # Set up plotting environment
  par ( mfrow = c(nS,1), mar = c(1,2,2,0), oma = c(4,4,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )
  # loop over species and plot biomass
  for( s in 1:nS )
  {
    plot    ( x = c(fYear,max(years)), y = c(0,1.1), type = "n",
              ylim = c(0,1.1), ylab = "", axes=FALSE, las = 1, xlab = ""  )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
      mtext( side = 2, text = specNames[s], line = 2.5, las = 0 )
      # OM
      polygon(  x = c(years,rev(years)),
                y = c(omBt[1,s,],rev(omBt[3,s,])),
                border = NA, col = "grey80" )
      if( est )
      {
        # SS envelope
        polygon( x = c( years, rev( years ) ),
                 y = c( ssBt[ 1, s, ], rev( ssBt[ 3, s, ] ) ),
                 border = NA, col = ssFill )
        # MS envelop
        polygon( x = c( years, rev( years ) ),
                 y = c( msBt[ 1, s, ], rev( msBt[ 3, s, ] ) ),
                 border = NA, col = msFill )
        # median
        lines( x = years, y = ssBt[ 2, s, ], col = ssCol, lwd = 2, lty = 2 )
        lines( x = years, y = msBt[ 2, s, ], col = msCol, lwd = 2, lty = 2 )

      }
      lines( x = years, y = omBt[ 2, s, ], lwd = 2, lty = 2 )
      abline( h = (1 - alpha / 2), lty = 3, lwd = .8 )
    
  }
  mtext( text = "Depletion", side = 2, line = 2.5, las = 0, outer = T, cex = 1.5 )
  mtext( text = "Year", outer = TRUE, side = 1, padj = 1.5, cex = 1.5, line = 2)
  if( devLabels )
    mtext ( text = c(stamp),side=1, outer = TRUE, 
            at = c(0.9),padj=2,col="grey50",cex=0.8 )
}


# plotqDists()
# Does what it says, plots to quartz the estimated and simulated q distributions
# for the complex in a given simulation and replicate.
# inputs:   sim = simulation number
#           rep = replicate number
# ouputs:   NULL
# side-eff: plots to quartz
# source: SDNJ
plotqDists <- function( sim = 1, rep = 1 )
{
  # load the simulation
  .loadSim(sim)

  # Now grab opMod q distribution pars
  qSurvM  <- blob$opMod$qSurvM
  qSurvSD <- blob$opMod$qSurvSD
  nSurv   <- blob$opMod$nSurv
  nStocks <- blob$opMod$nS
  # grab simulated q 
  simq    <- blob$om$q_os[rep,,]
  # grab both SS and MS estimated q
  estqSS  <- blob$am$ss$q_os[rep,,]
  estqMS  <- blob$am$ms$q_os[rep,,]
  # And MS estimated distribution
  estqSurvM   <- blob$am$ms$qbar[rep,]
  estqSurvSD  <- sqrt(blob$am$ms$tauq2[rep,])

  mpLabel <- blob$ctrl$mpLabel

  # Split plotting window
  par( mfrow = c(nSurv,1), mar = c(1,1,1,1), oma = c(3,3,1,1) )
  # Calculate range of q values
  qRange  <- range( simq, estqSS, estqMS )
  xLim    <- c( .8*qRange[1],1.2*qRange[2] )
  for( survIdx in 1:nSurv )
  {
    # Now make simulated and estimated distributions
    x       <- seq( xLim[1], xLim[2], length = 100 )
    qlNorm    <- dlnorm( x, meanlog = log(qSurvM[survIdx]), sdlog = qSurvSD[survIdx] )
    qlNormEst <- dlnorm( x, meanlog = log(estqSurvM[survIdx]), sdlog = estqSurvSD[survIdx] )
    # Initialise plots
    plot( x = xLim, y = c( 0, 1.2*max(qlNorm) ), type = "n", axes = F, xlab = "", ylab = "" )
      axis( side = 1 )
      panLab( x = .8, y = 0.8, txt = paste("Survey ", survIdx, sep = "") )
      # Plot distributions
      lines( x = x, y = qlNorm, lty = survIdx + 1, lwd = 3 )
      lines( x = x, y = qlNormEst, lty = survIdx + 1, lwd = 3, col = "grey60" )
      # Plot simulated survey mean q
      abline( v = qSurvM[survIdx], lty = survIdx + 1 )
      # Now plot the simulated values for each stock
      text( x = simq[survIdx,], y = (1:nStocks)*max(qlNorm)/(nStocks + 1),
            labels = 1:nStocks )
      # Joint SS and MS estimates with a thin line
      segments( x0 = estqSS[survIdx,], y0 = (1:nStocks)*max(qlNorm)/(nStocks + 1),
                x1 = estqMS[survIdx,], y1 = (1:nStocks)*max(qlNorm)/(nStocks + 1),
                lwd = .8, lty = 4 )
      # Plot SS and MS estimates
      points( x = estqSS[survIdx,], y = (1:nStocks)*max(qlNorm)/(nStocks + 1), pch = 16 )
      points( x = estqMS[survIdx,], y = (1:nStocks)*max(qlNorm)/(nStocks + 1), pch = 17 )
      abline( v = estqSurvM[survIdx], lty = survIdx + 1 + nSurv )
  }
  mtext( side = 1, text = "Catchability", outer = T, line = 1 )
  mtext( side = 1, line = 2, adj = 1, text = mpLabel )

  
} 

# Plots simulation tulip plots of biomass, catch and harvest
# rate
plotBCUenvelopes <- function( sim = 1,
                              fYear = 1988,
                              est = FALSE,
                              labSize = 2, 
                              tickSize = 1.2,
                              devLabels = TRUE )
{
  # First, load blob
  .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")

  # Recover blob elements for plotting
  nS <- blob$opMod$nS

  # Year indexing
  if(!is.null(blob$opMod$fYear)) 
  {
    nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
    fYear <- min(blob$opMod$fYear)
  } else nT <- blob$opMod$nT

  years <- fYear:(fYear + nT - 1)


  # Species names
  specNames <- blob$ctrl$speciesNames

  # True OM quantities
  omBt  <- blob$om$Bt
  Ct    <- blob$om$Ct
  Ut    <- blob$om$Ut

  # Single species model
  ssBt  <- blob$am$ss$Bt

  # Multispecies model
  msBt  <- blob$am$ms$Bt
  msBt[msBt < 0] <- NA

  years <- fYear:(fYear + nT - 1)

  # Estimated Ut
  ssUt <- Ct / ssBt
  msUt <- Ct / msBt

  ## Now scale quantities to eqbm values
  # Pull Bmsy and Umsy from the blob
  Umsy <- blob$opMod$Umsy
  Bmsy <- blob$opMod$Bmsy
  MSY  <- Umsy*Bmsy

  
  for( s in 1:nS )
  {
    omBt[,s,]  <- omBt[,s,]/Bmsy[s]/2
    Ct[,s,]    <- Ct[,s,]/MSY[s]
    Ut[,s,]    <- Ut[,s,]/Umsy[s]
    ssBt[,s,]  <- ssBt[,s,]/Bmsy[s]/2
    msBt[,s,]  <- msBt[,s,]/Bmsy[s]/2
    msUt[,s,]    <- msUt[,s,]/Umsy[s]
    ssUt[,s,]    <- ssUt[,s,]/Umsy[s]
  }  
  
  # Set colours
  ssCol <- "steelblue"
  msCol <- "salmon"
  # Create transparent fill colours for polys
  ssFill <- col2rgb(ssCol)/255
  ssFill <- rgb(ssFill[1],ssFill[2],ssFill[3],alpha=0.4)
  msFill <- col2rgb(msCol)/255
  msFill <- rgb(msFill[1],msFill[2],msFill[3],alpha=0.4)

  ## CHECK THAT THESE ARRAY CALCS ARE WORKING ##
  # Now, reduce to envelopes
  # OM values
  omBt  <- apply( X = omBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  Ct    <- apply( X = Ct, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  omUt  <- apply( X = Ut, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.5), na.rm = TRUE )
  # SS AM 
  ssBt  <- apply( X = ssBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  ssUt  <- apply( X = ssUt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )

  # MS AM 
  msBt  <- apply( X = msBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  msUt  <- apply( X = msUt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )

  # Set up plotting environment
  par ( mfrow = c(3,nS), mar = c(1,3,2,0), oma = c(3,3,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )
  # loop over species and plot biomass
  for( s in 1:nS )
  {
    plot    ( x = c(fYear,max(years)), y = c(0,1.1), type = "n",
              ylim = c(0,1.1), ylab = "", axes=FALSE, las = 1, xlab = "" ,
              main = specNames[s] )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
      # OM
      polygon(  x = c(years,rev(years)),
                y = c(omBt[1,s,],rev(omBt[3,s,])),
                border = NA, col = "grey80" )
      if( est )
      {
        # SS envelope
        polygon( x = c( years, rev( years ) ),
                 y = c( ssBt[ 1, s, ], rev( ssBt[ 3, s, ] ) ),
                 border = NA, col = ssFill )
        # MS envelop
        polygon( x = c( years, rev( years ) ),
                 y = c( msBt[ 1, s, ], rev( msBt[ 3, s, ] ) ),
                 border = NA, col = msFill )
        # median
        lines( x = years, y = ssBt[ 2, s, ], col = ssCol, lwd = 2, lty = 2 )
        lines( x = years, y = msBt[ 2, s, ], col = msCol, lwd = 2, lty = 2 )

      }
      lines( x = years, y = omBt[ 2, s, ], lwd = 2, lty = 2 )
    if( s == 1 ) mtext( text = "Depletion", side = 2, line = 2.5, las = 0 )
  }

  # Now catch
  for( s in 1:nS )
  {
    plot    ( x = c(fYear,max(years)), y = c(0,2.1), type = "n",
              ylim = c(0,2.1), ylab = "", axes=FALSE, las = 1, xlab = "" )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
      # OM
      polygon(  x = c(years,rev(years)),
                y = c(Ct[1,s,],rev(Ct[3,s,])),
                border = NA, col = "grey80" )
      lines( x = years, y = Ct[ 2, s, ], lwd = 2, lty = 2 )
    if( s == 1 ) mtext( text = "Ct/MSY", side = 2, line = 2.5, las = 0 )
  }

  # Lastly, exploitation rate
  for( s in 1:nS )
  {
    plot    ( x = c(fYear,max(years)), y = c(0,2.2), type = "n",
              ylim = c(0,2.2), ylab = "", axes=FALSE, las = 1, xlab = "" )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
      if( est )
      {
        # SS envelope
        polygon( x = c( years, rev( years ) ),
                 y = c( ssUt[ 1, s, ], rev( ssUt[ 3, s, ] ) ),
                 border = NA, col = ssFill )
        # MS envelope
        polygon( x = c( years, rev( years ) ),
                 y = c( msUt[ 1, s, ], rev( msUt[ 3, s, ] ) ),
                 border = NA, col = msFill )
        # medians
        lines( x = years, y = ssUt[ 2, s, ], lwd = 2, lty = 2, col = ssCol )
        lines( x = years, y = msUt[ 2, s, ], lwd = 2, lty = 2, col = msCol )

      }
      # OM
      lines( x = years, y = omUt[ s, ], lwd = 2, lty = 2 )
    if( s == 1 ) mtext( text = "Ut/Umsy", side = 2, line = 2.5, las = 0 )
  }
  mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5)
  if( devLabels )
    mtext ( text = c(stamp),side=1, outer = TRUE, 
            at = c(0.9),padj=2,col="grey50",cex=0.8 )
}

# plotBC()
# Plots a given replicate's true and estimated time series of biomass and 
# catch for all stocks
# inputs:   rep=replicate number for plotting; 
#           sim=number indicating blob to load (alpha order in project folder)
#           folder=name of folder/blob file (supercedes sim number)
# output:   NULL
# usage:    post-simulation run, plotting performance
plotBC <- function( rep = 1, sim=1, legend=TRUE,
                    data = FALSE, labSize = 2, tickSize = 1.2,
                    fYear = 1988, devLabels = TRUE, CIs = FALSE,
                    scale = FALSE )
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
  

  # Species names
  specNames <- blob$ctrl$speciesNames

  # True OM quantities
  omBt  <- blob$om$Bt[rep,,]
  Ct    <- blob$om$Ct[rep,,]
  It    <- blob$om$It[rep,,]
  Ut    <- blob$om$Ut[rep,,]
  
  # Single species model
  ssBt  <- blob$am$ss$Bt[rep,,]
  # ssq   <- blob$am$ss$q[rep,]

  # Multispecies model
  msBt  <- blob$am$ms$Bt[rep,,]
  msBt[msBt < 0] <- NA
  # msq   <- blob$am$ms$q[rep,]  

  # Year indexing
  if(!is.null(blob$opMod$fYear)) 
  {
    nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
    fYear <- min(blob$opMod$fYear)
  } else nT <- blob$opMod$nT
  # browser()
  years <- fYear:(fYear + nT - 1)

  # Estimated Ut
  ssUt <- Ct / ssBt
  msUt <- Ct / msBt

  if( CIs )
  {
    # arrays to hold CIs
    ssBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )
    msBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )

    msCIs <- blob$am$ms$CIs[[ rep ]]
    ssCIs <- blob$am$ss$CIs[[rep]]

    # browser()
    if( !is.na(msCIs) )
    {
      Btrows <- which( msCIs$par == "Bt" )
      msBio[ , , 1 ] <- matrix( msCIs[ Btrows , "uCI" ], nrow = nS, ncol = nT, byrow = FALSE )
      msBio[ , , 2 ] <- matrix( msCIs[ Btrows, "val" ], nrow = nS, ncol = nT, byrow = FALSE )
      msBio[ , , 3 ] <- matrix( msCIs[ Btrows, "lCI" ], nrow = nS, ncol = nT, byrow = FALSE ) 
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
  
  }

  # Pull Bmsy and Umsy from the blob
  Umsy <- blob$opMod$Umsy
  Bmsy <- blob$opMod$Bmsy
  MSY  <- Umsy*Bmsy

  if( scale )
  {
    # Now loop over species and scale
    for( s in 1:nS )
    {
      omBt[s,]  <- omBt[s,]/Bmsy[s]/2
      Ct[s,]    <- (Ct[s,]/MSY[s])
      Ut[s,]    <- Ut[s,]/Umsy[s]
      ssBt[s,]  <- ssBt[s,]/Bmsy[s]/2
      msBt[s,]  <- msBt[s,]/Bmsy[s]/2
      msUt[s,]    <- msUt[s,]/Umsy[s]
      ssUt[s,]    <- ssUt[s,]/Umsy[s]
      if(CIs)
      {
        ssBio[s,,]<- ssBio[s,,]/Bmsy[s]/2
        msBio[s,,]<- msBio[s,,]/Bmsy[s]/2  
      }
    }  
  }
  
  
  # Recover diagnostics for the fits
  hpdSS <- blob$am$ss$hesspd[rep,]
  grdSS <- blob$am$ss$maxGrad[rep,]
  hpdMS <- blob$am$ms$hesspd[rep]
  grdMS <- blob$am$ms$maxGrad[rep]

  # Set colours for each model
  ssCol <- "steelblue"
  msCol <- "salmon"

  # Set up plot window
  par ( mfrow = c(nS,1), mar = c(1,3,2,0), oma = c(3,4,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )
  # Plot biomass, actual and estimated, including 2 index series,
  # scaled by model estimated q
  maxBt <- 1.1*max(omBt,na.rm=T)
  if(CIs) maxBt <- max(maxBt,2)
  for ( s in 1:nS )
  {
    if( !scale ) maxBt <- 1.1*max ( omBt[s,] ,na.rm=TRUE)
    if( (!scale) & CIs ) maxBt <- max ( maxBt, msBio[s,,], ssBio[s,,], na.rm = T )
    plot    ( x = c(fYear,max(years)), y = c(0,maxBt), type = "n",
              ylim = c(0,maxBt), ylab = "", axes=FALSE, xlab = "" ,
              main = "" )
      # plot catch
      rect( xleft = years-.4, xright = years+.4,
            ybottom = 0, ytop = Ct[s,], col = "grey80", border = NA )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
    if (s == 1) panLegend ( x=0.2,y=1,legTxt=c("ss","ms"),
                            col=c(ssCol,msCol), lty = c(2,3), 
                            lwd = c(2,2), cex=c(1), bty="n" )
      panLab( x = .5, y = .9, txt = paste("Stock ", s, sep = "" ) )


    
    if( devLabels )
    {
      if ( hpdSS[s] & !is.na(hpdSS[s]) ) panLab (x=0.9,y=0.9,txt="h",col=ssCol,cex=1.1)
      if ( hpdMS & !is.na(hpdMS) ) panLab (x=0.9,y=0.85,txt="h",col=msCol,cex=1.1)  
    }
    if( data ) points  ( x = years, y = It[s,]/ssq[s], pch = 2, cex = 0.6, col="grey70" )
    if( data ) points  ( x = years, y = It[s,]/msq[s], pch = 5, cex = 0.6, col="grey70" )
    lines   ( x = years, y = omBt[s,], col = "black", lwd = 2)
    lines   ( x = years, y = ssBt[s,], col = ssCol, lwd = 2, lty = 2 )
    lines   ( x = years, y = msBt[s,], col = msCol, lwd = 2, lty = 3 )
    if( CIs )
    {
      lines ( x = years, y = ssBio[s,,1], col = ssCol, lwd = 1, lty = 2)
      lines ( x = years, y = ssBio[s,,3], col = ssCol, lwd = 1, lty = 2)
      lines ( x = years, y = msBio[s,,1], col = msCol, lwd = 1, lty = 3)
      lines ( x = years, y = msBio[s,,3], col = msCol, lwd = 1, lty = 3)
    }
  }
  mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5)
  mtext ( text = "Biomass (Kt)", outer = TRUE, side = 2, line = 2, las = 0 )
  if( devLabels )
    mtext ( text = c(stamp,repCount),side=1, outer = TRUE, 
            at = c(0.9,0.1),padj=2,col="grey50",cex=0.8 )
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
plotBCU <- function ( rep = 1, sim=1, legend=TRUE,
                      data = FALSE, labSize = 2, tickSize = 1.2,
                      fYear = 1988, devLabels = TRUE, CIs = FALSE,
                      scale = FALSE )
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
  

  # Species names
  specNames <- blob$ctrl$speciesNames

  # True OM quantities
  omBt  <- blob$om$Bt[rep,,]
  Ct    <- blob$om$Ct[rep,,]
  It    <- blob$om$It[rep,,]
  Ut    <- blob$om$Ut[rep,,]
  
  # Single species model
  ssBt  <- blob$am$ss$Bt[rep,,]
  # ssq   <- blob$am$ss$q[rep,]

  # Multispecies model
  msBt  <- blob$am$ms$Bt[rep,,]
  msBt[msBt < 0] <- NA
  # msq   <- blob$am$ms$q[rep,]  

  # Year indexing
  if(!is.null(blob$opMod$fYear)) 
  {
    nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
    fYear <- min(blob$opMod$fYear)
  } else nT <- blob$opMod$nT
  # browser()
  years <- fYear:(fYear + nT - 1)

  # Estimated Ut
  ssUt <- Ct / ssBt
  msUt <- Ct / msBt

  if( CIs )
  {
    # arrays to hold CIs
    ssBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )
    msBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )

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
  
  }

  # Pull Bmsy and Umsy from the blob
  Umsy <- blob$opMod$Umsy
  Bmsy <- blob$opMod$Bmsy
  MSY  <- Umsy*Bmsy

  if( scale )
  {
    # Now loop over species and scale
    for( s in 1:nS )
    {
      omBt[s,]  <- omBt[s,]/Bmsy[s]/2
      Ct[s,]    <- Ct[s,]/MSY[s]
      Ut[s,]    <- Ut[s,]/Umsy[s]
      ssBt[s,]  <- ssBt[s,]/Bmsy[s]/2
      msBt[s,]  <- msBt[s,]/Bmsy[s]/2
      msUt[s,]    <- msUt[s,]/Umsy[s]
      ssUt[s,]    <- ssUt[s,]/Umsy[s]
      if(CIs)
      {
        ssBio[s,,]<- ssBio[s,,]/Bmsy[s]/2
        msBio[s,,]<- msBio[s,,]/Bmsy[s]/2  
      }
    }  
  }
  
  
  # Recover diagnostics for the fits
  hpdSS <- blob$am$ss$hesspd[rep,]
  grdSS <- blob$am$ss$maxGrad[rep,]
  hpdMS <- blob$am$ms$hesspd[rep]
  grdMS <- blob$am$ms$maxGrad[rep]

  # Set colours for each model
  ssCol <- "steelblue"
  msCol <- "salmon"

  # Set up plot window
  par ( mfrow = c(3,nS), mar = c(1,3,2,0), oma = c(3,3,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )
  # Plot biomass, actual and estimated, including 2 index series,
  # scaled by model estimated q
  maxBt <- 1.1*max(omBt,na.rm=T)
  if(CIs) maxBt <- max(maxBt,2)
  for ( s in 1:nS )
  {
    if( !scale ) maxBt <- 1.1*max ( omBt[s,] ,na.rm=TRUE)
    if( (!scale) & CIs ) maxBt <- max ( maxBt, msBio[s,,], ssBio[s,,], na.rm = T )
    plot    ( x = c(fYear,max(years)), y = c(0,maxBt), type = "n",
              ylim = c(0,maxBt), ylab = "", axes=FALSE, las = 1, xlab = "" ,
              main = specNames[s] )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
    if (s == 1) panLegend ( x=0.2,y=1,legTxt=c("ss","ms"),
                            col=c(ssCol,msCol), lty = c(2,3), 
                            lwd = c(2,2), cex=c(1), bty="n" )
    if( scale )
    {
      ylabBio <- "(a)"
      labLas <- 1
    } else {
      ylabBio <- "Biomass (kt)"
      labLas <- 0
    }
    if( s == 1 ) mtext( text = ylabBio, side = 2, line = 2.5, las = labLas )
    if( devLabels )
    {
      if ( hpdSS[s] & !is.na(hpdSS[s]) ) panLab (x=0.9,y=0.9,txt="h",col=ssCol,cex=1.1)
      if ( hpdMS & !is.na(hpdMS) ) panLab (x=0.9,y=0.85,txt="h",col=msCol,cex=1.1)  
    }
    if( data ) points  ( x = years, y = It[s,]/ssq[s], pch = 2, cex = 0.6, col="grey70" )
    if( data ) points  ( x = years, y = It[s,]/msq[s], pch = 5, cex = 0.6, col="grey70" )
    lines   ( x = years, y = omBt[s,], col = "black", lwd = 2)
    lines   ( x = years, y = ssBt[s,], col = ssCol, lwd = 2, lty = 2 )
    lines   ( x = years, y = msBt[s,], col = msCol, lwd = 2, lty = 3 )
    if( CIs )
    {
      lines ( x = years, y = ssBio[s,,1], col = ssCol, lwd = 1, lty = 2)
      lines ( x = years, y = ssBio[s,,3], col = ssCol, lwd = 1, lty = 2)
      lines ( x = years, y = msBio[s,,1], col = msCol, lwd = 1, lty = 3)
      lines ( x = years, y = msBio[s,,3], col = msCol, lwd = 1, lty = 3)
    }
  }
  # Now plot catch
  maxCt <- max(Ct, na.rm=T)
  for ( s in 1:nS )
  {
    if( !scale ) maxCt <- max ( Ct[s,])
    plot  ( x = years, y = Ct[s,], col = "blue", lwd = 2, type = "l",
            ylim = c(0,maxCt), ylab = "", axes=FALSE, las = 1, xlab = "" )
      axis( side = 1 )
      axis( side = 2, las = 1 )
    if( scale )
    {
      ylabCt <- "(b)"
      labLas <- 1
    } else {
      ylabCt <- "Catch (kt)"
      labLas <- 0
    }
    if( s == 1 ) mtext( text = ylabCt, side = 2, line = 2.5, las = labLas )

  }
  # Now F
  maxUt <- max(Ut, na.rm=T)
  for ( s in 1:nS )
  {
    if(!scale) maxUt <- max ( Ut[s,])
    plot  ( x = years, y = Ut[s,], col = "black", lwd = 2, type = "l",
            ylim = c(0,maxUt), ylab = "", axes=FALSE, las = 1, xlab = "" )
      axis( side = 1 )
      axis( side = 2, las = 1 )
    if( scale )
    {
      ylabUt <- "(c)"
      labLas <- 1
    } else 
    {
      ylabUt <- "Exploitation Rate"
      labLas <- 0
    }
    if( s == 1 ) mtext( text = ylabUt, side = 2, line = 2.5, las = labLas )
    lines (x = years, y = ssUt[s,], col= ssCol, lwd = 2, lty = 2 )
    lines (x = years, y = msUt[s,], col= msCol, lwd = 2, lty = 3 )
  }
  mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5)
  if( devLabels )
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
plotStockPerfMultiSim <- function ( pars = c("BnT","Bmsy","q_os","Umsy","dep"), 
                                    sims=1, spec = 1, nSurv = 2,
                                    devLabels = TRUE,
                                    title = TRUE, plotMARE = FALSE )
{
  # Set plotting window
  par (mfrow = c(1,length(sims)), mar = c(1,0,2,1), oma = c(3,0,2,0) )

  # Create a wrapper function for generating quantiles
  quantWrap <- function ( entry = 1, x = ssRE, ... )
  {
    if (is.null(dim(x[[entry]]))) x[[entry]] <- matrix(x[[entry]],ncol=1)
    xDim <- dim(x[[entry]])
    margins <- 2:length(xDim)
    quants <- apply ( X = x[[entry]], FUN = quantile, MARGIN = margins, ... )
    if(class(quants) == "array")
      return( aperm(quants) )
    else return( t(quants) )
  }  

  blobList <- vector(mode = "list", length=length(sims))
  # Loop over sims and load their MLEs
  for( sIdx in 1:length(sims) )
  { 
    sim <- sims[sIdx]
    .loadSim(sim)

    blobList[[sIdx]]$mp       <- blob$ctrl$mp
    blobList[[sIdx]]$scenario <- blob$ctrl$scenario
    blobList[[sIdx]]$stamp    <- paste( blob$ctrl$scenario, blob$ctrl$mp, sep = ":")
    blobList[[sIdx]]$nO       <- blob$opMod$nSurv

    blobList[[sIdx]]$ssQuant  <- lapply ( X = seq_along(blob$am$ss$err.mle[pars]), 
                                          FUN = quantWrap, x=blob$am$ss$err.mle[pars], 
                                          probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    names( blobList[[sIdx]]$ssQuant ) <- names( blob$am$ss$err.mle[pars] )

    
    blobList[[sIdx]]$msQuant  <- lapply ( X = seq_along(blob$am$ms$err.mle[pars]), 
                                          FUN = quantWrap, x=blob$am$ms$err.mle[pars], 
                                          probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

    names( blobList[[sIdx]]$msQuant ) <- names( blob$am$ms$err.mle[pars] )


  }

  # Multisurvey parameters
  if( "q_os" %in% pars )
  {
    q_os <- paste( "q_", 1:nSurv, sep = "" )
    pars <- c( pars[pars != "q_os"], q_os ) 
  }

  if( "tau2_o" %in% pars )
  {
    tau2_o <- paste( "tau2_", 1:nSurv, sep = "" )
    pars <- c( pars[pars != "tau2_o"], tau2_o ) 
  }

  # Create an array of quantile values 
  # (dims = species,percentiles,parameters,models)
  quantiles <- array ( NA, dim = c(length(sims),3,length(pars),2))
  dimnames (quantiles) <- list (  character(length(sims)),
                                  c("2.5%","50%","97.5%"), 
                                  pars,
                                  c("ss","ms") )


  MPs   <- "Single Stock"
  for( sIdx in 1:length(sims) )
  {
    # populate table
    for (p in 1:length(pars))
    {
      parName <- pars[p]
      # Change behaviour if survey specific stuff is involved
      # SDNJ: Check backwards compatibility for master branch
      if( grepl( "q_", parName ) | grepl("tau2_", parName ) )
      {
        if( grepl( "q_", parName ) ) parName <- "q_os"
        if( grepl("tau2_", parName ) ) parName <- "tau2_o"  
        surveyNo <- as.numeric(str_split(pars[p],"_")[[1]][2])
        # Now fill in the survey specific deets
        # REs
        quantiles[sIdx,,pars[p],1] <- blobList[[sIdx]]$ssQuant[ parName ][[1]][spec,surveyNo,]
        quantiles[sIdx,,pars[p],2] <- blobList[[sIdx]]$msQuant[ parName ][[1]][spec,surveyNo,]
        # # MARE
        # MARE[s,pars[p],1] <- ssAQuant[[ parName ]][s,surveyNo]
        # MARE[s,pars[p],2] <- msAQuant[[ parName ]][s,surveyNo]
      } else {
        quantiles[sIdx,,pars[p],1] <- blobList[[sIdx]]$ssQuant[ parName ][[1]][spec,]
        quantiles[sIdx,,pars[p],2] <- blobList[[sIdx]]$msQuant[ parName ][[1]][spec,]
        # # MARE
        # MARE[s,pars[p],1] <- ssAQuant[[ parName ]][s]
        # MARE[s,pars[p],2] <- msAQuant[[ parName ]][s]
      }
    }

    MPs <- c(MPs,blobList[[sIdx]]$mp)
  }
  # Pull out RE dists
  med <- quantiles[ ,"50%",, ]
  q975 <- quantiles[ ,"97.5%",, ]
  q025 <- quantiles[ ,"2.5%",, ]

  # Set up plotting environment
  par( mfrow = c(length(pars),1), mar = c(1,2,1,2), oma = c(3,7,1,6) )
  xlim  <- c(-1.5,1.5)  
  # Loop over pars
  for( pIdx in 1:length(pars) )
  {
    par <- pars[pIdx]
    plot( x = xlim, y = 2*c(0,length(MPs)+1), type = "n", axes = F )
      if( pIdx == length(pars) ) axis( side = 1 )
      abline( v = 0, lty = 2, lwd = 1 )
      axis( side = 2 , at = 2*1:length(MPs), labels = MPs, las = 1, cex.axis = 0.8 )
      segments( x0 = q025[1,par,"ss"], y0 = 2, x1 = q975[1,par,"ss"], y1 = 2, lty=1, col="grey60", lwd = 3 )
      points( x = med[1,par,"ss"], y = 2, pch = 16 )
      segments( x0 = q025[,par,"ms"], y0 = 2*(2:length(MPs)), x1 = q975[,par,"ms"], y1 = 2*(2:length(MPs)), lty=1, col="grey60", lwd = 3 )
      points( x = med[,par,"ms"], y = 2*(2:length(MPs)), pch = 16 )
      mtext( side = 4, text = par, las = 1, line = 2 )
  }

  if( title )
  {
    mtext ( text = "Relative Error", side = 1, outer = TRUE, line = 2, cex = 1.5)
    mtext ( text = "AM configuration", side = 2, outer = TRUE, line = 5, cex = 1.5)
    mtext ( text = "Parameter", side = 4, outer = TRUE, line = 4, cex = 1.5, las = 0)
  }

}


# plotSimPerf()
# A function to plot simulation-estimation performance for a whole
# collection of replicates. Resulting plot is a dot and bar graph
# showing relative error distributions for nominated parameters
# inputs:   sim = numeric indicator of simulation in project folder
#           folder = optional character indicating sim folder name
#           pars = nominated estimated leading and derived parameters 
# outputs:  NULL
plotSimPerf <- function ( pars = c("Bmsy","Umsy","q_os","dep","BnT","tau2_o","totRE"), sim=1, 
                          est="MLE", devLabels = TRUE,
                          title = TRUE, plotMARE = FALSE )
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
  specNames <- blob$ctrl$speciesNames[1:nS]
  nO        <- ifelse(!is.null(blob$opMod$nSurv), blob$opMod$nSurv, 1)

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
    xDim <- dim(x[[entry]])
    margins <- 2:length(xDim)
    quants <- apply ( X = x[[entry]], FUN = quantile, MARGIN = margins, ... )
    if(class(quants) == "array")
      return( aperm(quants) )
    else return( t(quants) )
  }

  # browser()
  ssARE <- lapply(X = ssRE, FUN = abs)
  msARE <- lapply(X = msRE, FUN = abs)

  # generate quantiles
  ssQuant <- lapply ( X = seq_along(ssRE), FUN = quantWrap, x=ssRE, 
                      probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  msQuant <- lapply ( X = seq_along(msRE), FUN = quantWrap, x=msRE, 
                      probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

  ssAQuant <- lapply (  X = seq_along(ssARE), FUN = quantWrap, x=ssARE, 
                        probs = c(0.5), na.rm = TRUE)
  msAQuant <- lapply (  X = seq_along(msARE), FUN = quantWrap, x=msARE, 
                        probs = c(0.5), na.rm = TRUE)

  # get names of parameters
  ssPars <- names ( ssRE )
  msPars <- names ( msRE )
  names ( ssQuant) <- ssPars
  names ( msQuant) <- msPars

  names ( ssAQuant) <- ssPars
  names ( msAQuant) <- msPars

  # Multisurvey parameters
  if( "q_os" %in% pars )
  {
    q_os <- paste( "q_", 1:nO, sep = "" )
    pars <- c( pars[pars != "q_os"], q_os ) 
  }

  if( "tau2_o" %in% pars )
  {
    tau2_o <- paste( "tau2_", 1:nO, sep = "" )
    pars <- c( pars[pars != "tau2_o"], tau2_o ) 
  }

  # Create an array of quantile values 
  # (dims = species,percentiles,parameters,models)
  quantiles <- array ( NA, dim = c(nS,3,length(pars),2))
  dimnames (quantiles) <- list (  specNames,
                                  c("2.5%","50%","97.5%"), 
                                  pars,
                                  c("ss","ms") )

  MARE <- array ( NA, dim = c(nS,length(pars),2),
                  dimnames = list ( specNames, pars,
                                    c("ss","ms") ) )


  # populate table
  for ( s in 1:nS )
  {
    for (p in 1:length(pars))
    {
      parName <- pars[p]
      # Change behaviour if survey specific stuff is involved
      # SDNJ: Check backwards compatibility for master branch
      if( grepl( "q_", parName ) | grepl("tau2_", parName ) )
      {
        if( grepl( "q_", parName ) ) parName <- "q_os"
        if( grepl("tau2_", parName ) ) parName <- "tau2_o"  
        surveyNo <- as.numeric(str_split(pars[p],"_")[[1]][2])
        # Now fill in the survey specific deets
        # REs
        quantiles[s,,pars[p],1] <- ssQuant[ parName ][[1]][s,surveyNo,]
        quantiles[s,,pars[p],2] <- msQuant[ parName ][[1]][s,surveyNo,]
        # MARE
        MARE[s,pars[p],1] <- ssAQuant[[ parName ]][s,surveyNo]
        MARE[s,pars[p],2] <- msAQuant[[ parName ]][s,surveyNo]
      } else {
        quantiles[s,,pars[p],1] <- ssQuant[ parName ][[1]][s,]
        quantiles[s,,pars[p],2] <- msQuant[ parName ][[1]][s,]
        # MARE
        MARE[s,pars[p],1] <- ssAQuant[[ parName ]][s]
        MARE[s,pars[p],2] <- msAQuant[[ parName ]][s]
      }
    }
  }
  # Set plotting window
  par (mfrow = c(1,nS), mar = c(1,0,2,1), oma = c(3,0,2,0) )

  for ( s in 1:nS )
  {
    med <- quantiles[s,"50%",,]
    mare <- MARE[s,,]
    q975 <- quantiles[s,"97.5%",,]
    q025 <- quantiles[s,"2.5%",,]

    # Plot main dotchart
    plotTitle <- paste(specNames[s])
    FitCounter <- paste( "nReps = ", sum(blob$goodReps[,s]), sep = "" )
    if ( s == 1) dotchart ( x = t(med), xlim = c(-1.5,1.5), main = "",
                            pch = 16, cex = 1.2 )
    else dotchart ( x = t(med), xlim = c(-1.5,1.5), main = "",
                    pch = 16, cex = 1.2)
      mtext( side = 3, text = plotTitle, cex = 1.5, las = 0, line = 2 )
      mtext( side = 3, text = FitCounter, cex = 1, las = 0, line = 1 )
      # axis ( side = 1, at = seq(-1.5,1.5,by=0.5) )
    
    # Now add segments
    for ( p in length(pars):1)
    {
      parIdx <- length(pars) - p + 1
      plotY <- 2 * p + 2 * (p-1)
      segments( q025[pars[parIdx],"ss"],plotY-1,q975[pars[parIdx],"ss"],plotY-1,lty=1,col="grey60", lwd=3)  
      segments( q025[pars[parIdx],"ms"],plotY,q975[pars[parIdx],"ms"],plotY,lty=1,col="grey60", lwd=3)  
      if(plotMARE) 
      {
        points(x = mare[parIdx,"ss"],y = plotY-1, pch = 2, cex = 1.2 )
        points(x = mare[parIdx,"ms"],y = plotY, pch = 2, cex = 1.2 )
      }
    }

    abline ( v = 0, lty = 3, lwd = 0.5 )

  } 
  if( title )
  {
    mtext ( text = "Relative Error", side = 1, outer = TRUE, line = 2, cex = 1.5)
    mtext ( text = "Parameter", side = 2, outer = TRUE, line = 1, cex = 1.5)
  }
  
  if( devLabels )
    mtext ( text = c(stamp),side=1, outer = TRUE, 
            at = c(0.75),padj=1,col="grey50" )

}


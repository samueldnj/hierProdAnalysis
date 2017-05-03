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

checkCorr <- function(  tabName = "rKqExp.csv",
                        cols = c("ssq","ssUmsy","ssBmsy"),
                        spec = NULL )
{
  # first read in the table
  tabPath <- file.path(getwd(),"project/statistics",tabName)
  table <- read.csv( tabPath, header=TRUE, stringsAsFactors=FALSE ) 

  # restrict to correct number of species
  if( !is.null(spec) ) table <- table %>% filter( species == spec )

  # Now compute the correlation
  corrMtx <- cor( table[, cols], use = "pairwise.complete.obs" )

  corrMtx
}


groupPars <- c("qbar","tauq2","Umsybar","sigU2")

metaModels <- function( tabName = "allSame_RE_msIncr__MRE",
                        multiResp = c("BnT","Umsy","q","Dep","Bmsy"),
                        singleResp = groupPars,
                        spec = c("Stock1"),
                        expVars = c("kappaMult","corr","mp","nS"),
                        sig = .05, intercept = FALSE,
                        scaled = TRUE, saveOut = TRUE, interactions = TRUE )
{
  # Create an rDataFile output name
  rDataName <- tabName
  ssResp <- paste("ss", multiResp, sep = "" )
  msResp <- paste("ms", multiResp, sep = "" )
  
  # Now some lists to hold meta-models
  fitList <- vector( mode = "list", length = 3 )
  names(fitList) <- c("ss","ms","group")
  # responses that have single species and multispecies versions
  respList <- vector(mode = "list", length = length(multiResp))
  names( respList ) <- multiResp
  # Group parameters (MS prior mean/var, convergence stats )
  groupList <- vector(mode = "list", length = length(singleResp))
  fitList$ss <- respList
  fitList$ms <- respList
  fitList$group <- groupList
  
  # Now loop over response variables and fit models
  for( rIdx in 1:length(multiResp) )
  {
    resp <- multiResp[rIdx]
    ssResp <- paste( "ss", resp, sep = "" )
    msResp <- paste( "ms", resp, sep = "" )
    fitList$ss[[rIdx]] <- .DASEmodelSelection(  tabName = tabName,
                                                resp = ssResp,
                                                expVars = expVars,
                                                sig = sig,
                                                intercept = intercept,
                                                scaled = scaled,
                                                abs = FALSE,
                                                spec = spec,
                                                interactions = interactions )
    fitList$ms[[rIdx]] <- .DASEmodelSelection(  tabName = tabName,
                                                resp = msResp,
                                                expVars = expVars,
                                                sig = sig,
                                                intercept = intercept,
                                                scaled = scaled,
                                                abs = FALSE,
                                                spec = spec,
                                                interactions = interactions )
  }
  if(!is.null(singleResp))
  {
    for( rIdx in 1:length(singleResp) )
    {
      resp <- singleResp[rIdx]
      fitList$group[[rIdx]] <- .DASEmodelSelection( tabName = tabName,
                                                    resp = resp,
                                                    expVars = expVars,
                                                    sig = sig,
                                                    intercept = intercept,
                                                    scaled = scaled,
                                                    abs = FALSE,
                                                    interactions = interactions )
    }  
  } 
  savePath <- file.path( getwd(), "project", "Statistics", rDataName ) 
  dir.create( savePath )
  if(saveOut) save( fitList, 
                    file = file.path( savePath,
                                      paste(rDataName,".RData", sep = "" ) ) )

  # Now output the aic tables
  # First, the ss/ms pars
  for( par in multiResp)
  {
    # recover ss and ms AIC tables
    AICss <- fitList$ss[[par]]$AICrank
    AICms <- fitList$ms[[par]]$AICrank
    # Main effects first
    AICssFile <- paste(rDataName,par,"AICss.csv",sep = "")
    AICmsFile <- paste(rDataName,par,"AICms.csv",sep = "")
    
    write.csv( x = AICss, file = file.path(savePath,AICssFile ) )
    write.csv( x = AICms, file = file.path(savePath,AICmsFile ) )
  }
  # then the group pars
  for( par in singleResp)
  {
    # recover ss and ms AIC tables
    AIC <- fitList$group[[par]]$AICrank
    # Main effects first
    AICfile <- paste(rDataName,par,"AIC.csv",sep = "")
    
    write.csv( x = AIC, file = file.path(savePath,AICfile ) )
  }
  fitList
}

# .DASEmodelSelection()
# Uses the following procedure to select candidate models for
# analysing simulation model output.
# Basic idea, uses strong heredity of main effects: 
#   1. Fit all main effects, check for significance
#   2. Choose main effects model with lowest AIC value from among significant models, 
#       repeat 1 with all possible interactions and second order effects from that model
#   3. Choose among first order, fo+interactions, fo+second order and
#       fo+so+interactions models using AICc
# Models from out$ms$sel and out$ss$sel can be used in plotOMSpecResponse 
# to plot the glm fits next to the observed values.
.DASEmodelSelection <- function(  tabName = "rKq_msInc__MRE",
                                  resp = "ssq",
                                  spec = "Stock1",
                                  expVars = c("qOM","UmsyOM","BmsyOM","mp"),
                                  sig = 0.1,
                                  intercept = FALSE,
                                  abs = FALSE,
                                  scaled = TRUE,
                                  interactions = FALSE )
{
  # First, fit main effects
  mainEff <- .DASEexperiment( tabName = tabName,
                              resp = resp, spec = spec,
                              scaled = scaled,
                              baseRHS = NULL,
                              abs = abs,
                              newExp = expVars,
                              intercept = intercept )

  if( interactions )
  {
    # Then, look at the ANOVA table for the 
    # factors, since some may be qualitative and have 
    # insignificant levels (small effect sizes)
    mainDrop <- lapply( X = mainEff$GLM, FUN = drop1, test = "LRT" )
    mainDrop <- mainDrop[[length(mainDrop)]]

    insig <- which( mainDrop[-1,5] >= sig )
    insig <- rownames(mainDrop[-1,])[insig]
    
    dropExp <- names(insig)

    signifExp <- expVars[ !(expVars %in% dropExp) ]

    # Function to make interactions out of main effects
    makeIntEff <- function ( main = expVars )
    {
      interactions <- character(length = 0)
      for( i in 1:(length(main)-1) )
      {
        for( j in (i+1):length(main) )
        {
          interactions <- c( interactions, paste( main[i], main[j], sep = ":" ) )
        }
      }
      interactions
    }
    # Now paste main effects together
    signifRHS <- paste(signifExp, collapse = " + ")

    # make interactions and second order effects
    intEff <- makeIntEff ( signifExp )
    # remove qualitative factors for second order effects
    if( "mp" %in% signifExp ) signifExpSO <- signifExp[ signifExp != "mp" ]
    else signifExpSO <- signifExp
    soEff  <- paste( "I(", signifExpSO, "^2)", sep = "" )

    # And run the experiment again
    fullModels <- .DASEexperiment(  tabName = tabName, resp = resp,
                                    spec = spec, scaled = FALSE,
                                    baseRHS = signifRHS,
                                    newExp = c(intEff, soEff) )

    allModels <- c(mainEff$GLM[-length(mainEff$GLM)],fullModels$GLM)
  } else allModels <- mainEff$GLM

  # Now rank AICc values in each model
  aicRank     <- AICrank( allModels, sig = sig )

  aicRank
}


# Automate ranking by AIC values, append effect sizes
# to the output.
AICrank <- function ( modelList, sig )
{
  # Count models
  nModels   <- length( modelList )
  summList  <- lapply( X = modelList, FUN = summary )
  dropList  <- lapply( X = modelList, FUN = drop1, test = "LRT" )

  # Pull out model with all coefficients for effect size recording
  fullModel <- modelList[[ nModels ]]
  fullSumm  <- summList[[ nModels ]]
  coefFull  <- coef( fullModel )
  coefFSumm <- coef( fullSumm )
  nCoef     <- nrow( coefFSumm )

  # Create a data frame to hold the info
  rank.df           <- matrix(NA, ncol = 4, nrow = nModels )
  colnames(rank.df) <- c( "Number", "AICc", "deltaAICc", "maxPr" )
  rank.df           <- as.data.frame(rank.df)

  # Make effect size df
  effect.df             <- matrix ( 0, ncol = nrow( coefFSumm ) + 1, nrow = nModels )
  colnames(effect.df)   <- c("Number", dimnames(coefFSumm)[[1]])
  effect.df             <- as.data.frame(effect.df)
  
  # Add an idenitifying number (useful for looking through model objects)
  rank.df[,1] <- 1:nModels
  effect.df[,1] <- 1:nModels
  for(k in 1:nModels)
  {
    nPar <- nrow(summList[[k]]$coefficients)
    nObs <- length(summList[[k]]$deviance.resid)
    rank.df[k,"AICc"] <- summList[[k]]$aic + nPar*(nPar + 1) / (nObs - nPar - 1)
    rank.df[k,"maxPr"] <- max(dropList[[k]][-1,5])
    coefSumm <- coef(summList[[k]])
    coefSumm[,1:2] <- round(coefSumm[,1:2] * 100,2)
    effects <- paste( coefSumm[,1], " (", coefSumm[,2],")", sep = "" )
    for( eIdx in 1:nrow( coefSumm ) )
    {
      effName <- dimnames(coefSumm)[[1]][eIdx]
      effect.df[k,effName] <- effects[eIdx]
    }
  } 

  effect.df[is.na(effect.df)] <- 0

  # Join effect sizes to ranks, and remove insignificant models
  rank.df <-  rank.df %>% 
              left_join( y = effect.df, by = "Number" ) %>%
              filter( maxPr < sig )


  rank.df[,"deltaAICc"] <- rank.df[,"AICc"] - min(rank.df[,"AICc"])
  rank.df <- rank.df[order(rank.df[,"deltaAICc"]),]
  
  return( list( AICrank = rank.df,
                models = modelList[rank.df[,"Number"]] )
        )
}


# .DASEexperiment()
# Function to apply the DASE methodology to the DASEex.csv table, my
# first attempt at discovering interactions between model inputs. Following
# exposition in Ch2 of Kleijnen, 2008: Design and Analysis of Simulation Experiments
# inputs:   tabName = character giving the root of the statTable
#           par = character indicating response variable parameter
#           species = character name of species to look at for effect
#           scaled = logical indicating whether to scale factor levels identically
#           rhs = character containing rhs of formula object for linear meta-model
# outputs:  fit = lm output for metamodel
.DASEexperiment <- function(  tabName = "lowK_rKq_RE",
                              resp = "q",
                              spec = "Dover",
                              scaled = TRUE,
                              baseRHS = NULL,
                              abs = FALSE,
                              intercept = FALSE,
                              newExp = c( "qOM", "UmsyOM", "BmsyOM" ) )
{
  tabFile <- paste(tabName, ".csv", sep = "")
  tabPath <- file.path(getwd(),"project/Statistics",tabFile)
  table <- read.csv( tabPath, header=TRUE, stringsAsFactors=FALSE ) 

  if( abs )
  {
    table[,resp] <- abs( table[,resp])
  }

  # restrict to nominated species/stock
  if(!is.null(spec)) 
    table <-  table %>%
              filter( species == spec )

  # Functions to rearrange the table using dplyr
  medTop <- function(x)
  {
    if(is.numeric(x)) out <- median(x)
    if(length(unique(x)) == 1) out <- x[1]

    out
  }

  # Take median of entries (MRE,MARE values), will just replicate the table
  # if using MRE/MARE data
  medTab <- table %>% 
            group_by(scenario,mp,species) %>%
            summarise_all(.funs=medTop) %>%
            ungroup()
  
  scaleNum <- function(x)
  {
    if(length(unique(x)) == 1) return(x[1])
    ranX <- range(x)
    gradX <- 2 / (ranX[2] - ranX[1])
    out <- gradX * x + (1 - gradX*ranX[2]) 
    out
  }
  pull    <- function(x,n) {x[[n]]}
  pullMP  <- function(x)
  {
    if(length(x) == 1) out <- x
    if(length(x) == 2) out <- x[[1]]
    if(length(x) == 3) out <- paste(x[[1]],x[[2]],sep = "_")
    out
  }

  cutScenario <- function(x)
  {
    x <- str_split(x,"_")
    x <- sapply(X = x, FUN = pull, n = 2)
    x
  }

  cutMP <- function(x)
  {
    x <- str_split(x,"_")
    x <- sapply(X = x, FUN = pullMP )
    x
  }

  # medTab <- file.path(getwd(),"project/statistics",medTab)
  # medTab <- read.csv( medTab, header=TRUE, stringsAsFactors=FALSE )   
  # medTab <- medTab %>% filter( species == spec )


  # Scale inputs
  if( scaled )
  {
    table <-  table %>%
              mutate( qOM = scaleNum(qOM),
                      UmsyOM = scaleNum(UmsyOM),
                      BmsyOM = scaleNum(BmsyOM),
                      tau2OM = scaleNum(tau2OM),
                      kappaTrue = scaleNum(kappaTrue),
                      corr = scaleNum(corr),
                      SigmaTrue = scaleNum(SigmaTrue) )
  }

  # Function to fit GLMs given...
  fitGLM <- function( rhs = expVars[1],
                      resp = ssPar,
                      dat = table )
  {
    form    <- as.formula( paste( resp, "~", rhs, sep = "" ) )
    glmObj  <- glm( formula = form, data = dat, na.action = "na.omit" )

    return(glmObj)
  }

  # Function to make interactions out of main effects
  makeIntEff <- function ( main = expVars )
  {
    interactions <- character(length = 0)
    for( i in 1:(length(main)-1) )
    {
      for( j in (i+1):length(main) )
      {
        interactions <- c( interactions, paste( main[i], main[j], sep = ":" ) )
      }
    }
    interactions
  }

  # Generate the sets of possible effects
  # First, take the power set of main effects and collapse into formulae
  newExpCombos  <- powerset( newExp )
  newExpCombos  <- lapply( X = newExpCombos, FUN = paste, collapse = " + " )
  rhsList       <- paste( baseRHS, newExpCombos, sep = " + " )
  if( !intercept ) rhsList <- paste( rhsList, "0", sep = " + ")
  # Fit to all possible main effects for both models
  GLM       <- lapply( X = rhsList, FUN = fitGLM, resp = resp, dat = table )

  # A function to lapply to the glmObjects and compute the AIC
  medResids   <- function( glmObj, data, resp )
  {
    pred <- predict( glmObj, newdata = data )
    resids <- data[,resp] - pred
    resids
  }

  # Now run predict the medTab from the models
  resids    <- lapply( X = GLM, FUN = medResids, data = medTab, resp = resp )


  return( list (  GLM = GLM, resids = resids,
                  medTab = medTab ) )
}

#.statTables()
# Wrapper for .statTableXXX functions, produces stacked tables of stats 
# for a group of simulations.
# inputs:   sims=integer vector indicating simulations in ./project/
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
.statTables <- function (sims=1,tabNameRoot = "statTable")
{ 
  # MSE
  # .statTableMSE(sims,paste(tabNameRoot,"_MSE.csv",sep=""))
  # MRE
  .statTableMRE(sims,paste(tabNameRoot,"_MRE.csv",sep=""))
  # MARE
  .statTableMARE(sims,paste(tabNameRoot,"_MARE.csv",sep=""))
  # raw RE distributions
  .statTableRE(sims,paste(tabNameRoot,"_RE.csv", sep = ""))
}

#.statTableRE()
# Wrapper for .simStatRE, produces stacked tables of relative errors 
# for a group of simulations.
# inputs:   sims=integer vector indicating simulations in ./project/
# outputs:  statTable=table of statistics for a project/group
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
.statTableRE <- function (sims=1,tabName = "statTable.csv")
{ 
  # call function
  tableList <- lapply ( X = sims, FUN = .simStatRE )

  # now make the table and return
  statTable <-  do.call("rbind",tableList)
  savePath <- file.path(getwd(),"project","Statistics",tabName)
  write.csv ( statTable, file = savePath )
  statTable
}

#.statTableMSE()
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

#.statTableMRE()
# Wrapper for .simStats, produces stacked tables of stats for a group
# of simulations.
# inputs:   sims=integer vector indicating simulations in ./project/
# outputs:  statTable=table of statistics for a project/group
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
.statTableMRE <- function (sims=1,tabName = "statTable.csv")
{ 
  # call function
  tableList <- lapply ( X = sims, FUN = .simStatMRE )

  # now make the table and return
  statTable <-  do.call("rbind",tableList)
  savePath <- file.path(getwd(),"project","Statistics",tabName)
  write.csv ( statTable, file = savePath )
  statTable
}

#.statTableMARE()
# Wrapper for .simStats, produces stacked tables of stats for a group
# of simulations.
# inputs:   sims=integer vector indicating simulations in ./project/
# outputs:  statTable=table of statistics for a project/group
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
.statTableMARE <- function (sims=1,tabName = "statTable.csv")
{ 
  # call function
  tableList <- lapply ( X = sims, FUN = .simStatMARE )

  # now make the table and return
  statTable <-  do.call("rbind",tableList)
  savePath <- file.path(getwd(),"project","Statistics",tabName)
  write.csv ( statTable, file = savePath )
  statTable
}

# .simStatMRE()
# Produces a statistics table for leading pars in a simulation from
# the output produced by a runSimEst() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame of mean squared error BnT and Umsy
# usage:    in lapply to produce stats for a group of simulations
.simStatMARE <- function ( sim=1 )
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
  species <- blob$ctrl$speciesName[1:nS]
  nReps   <- blob$ctrl$nReps

  # get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  success <- blob$goodReps

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","species","kappaTrue",
                  "SigmaTrue", "kappaMult", "ssBnT","msBnT","ssUmsy","msUmsy",
                  "ssBmsy","msBmsy","ssDep","msDep",
                  "ssq","msq", "msHessPD", "ssHessPD","nReps",
                  "Umax", "tUpeak", "tUtrough", "tau2OM" )
  
  statTable <- matrix(NA,nrow=nS,ncol=length(colLabels))
  
  colnames(statTable)   <- colLabels
  statTable             <- as.data.frame(statTable)

  # get multiplier for shared effects
  kappaMult <- opMod$kappaMult
  if (is.null(opMod$kappaMult)) kappaMult <- 1

  # Start filling stat table
  # First, info and true pars
  statTable$scenario        <- blob$ctrl$scenarioName
  statTable$mp              <- blob$ctrl$mpLabel
  statTable$species         <- species[1:nS]
  statTable$kappaTrue       <- ifelse(is.null(opMod$kappa2),sqrt(opMod$pars$kappa2),sqrt(opMod$kappa2))*kappaMult
  statTable$SigmaTrue       <- ifelse(is.null(opMod$SigmaDiag),sqrt(opMod$pars$Sigma2[1:nS]),sqrt(opMod$SigmaDiag[1:nS]))
  statTable$kappaMult       <- kappaMult
  statTable$corr            <- ifelse(is.null(opMod$corrOffDiag),opMod$corrMult,opMod$corrOffDiag)
  statTable$nReps           <- length(success)
  statTable$Umax            <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUtrough        <- opMod$tUtrough
  statTable$tUpeak          <- opMod$tUpeak
  statTable$tau2OM          <- opMod$tau2[1:nS]
  statTable$lastNegCorr     <- ifelse(is.null(opMod$lastNegCorr),ifelse(nS==5,TRUE,FALSE),opMod$lastNegCorr)
  statTable$nS              <- opMod$nS
  statTable$UmsyOM          <- opMod$Umsy[1:nS]
  statTable$BmsyOM          <- opMod$Bmsy[1:nS]
  statTable$qOM             <- opMod$q[1:nS]
  statTable$s2q             <- blob$assess$s2lnq
  statTable$fixqbar         <- blob$assess$fixqbar
  statTable$fixtauq2        <- blob$assess$fixqbar
  statTable$s2Umsy          <- blob$assess$s2lnUmsy
  if (is.null(blob$assess$tauq2Prior))
  {
    statTable$tauq2P1        <- ifelse(is.null(blob$assess$tauq2IGa),blob$assess$tauq2IG[1],blob$assess$tauq2IGa)
    statTable$tauq2P2        <- ifelse(is.null(blob$assess$tauq2IGb),blob$assess$tauq2IG[2],blob$assess$tauq2IGb)
    statTable$sigU2P1        <- ifelse(is.null(blob$assess$sigU2IGa),blob$assess$sigU2IG[1],blob$assess$sigU2IGa)
    statTable$sigU2P2        <- ifelse(is.null(blob$assess$sigU2IGb),blob$assess$sigU2IG[2],blob$assess$sigU2IGb)  
  } else {
    statTable$tauq2P1         <- ifelse(is.null(blob$assess$tauq2P1),blob$assess$tauq2Prior[1],blob$assess$tauq2P1)
    statTable$tauq2P2         <- ifelse(is.null(blob$assess$tauq2P2),blob$assess$tauq2Prior[2],blob$assess$tauq2P2)
    statTable$sigU2P1         <- ifelse(is.null(blob$assess$sigU2P1),blob$assess$sigU2Prior[1],blob$assess$sigU2P1)
    statTable$sigU2P2         <- ifelse(is.null(blob$assess$sigU2P2),blob$assess$sigU2Prior[2],blob$assess$sigU2P2) 
  }
  statTable$sigmaPriorCode  <- blob$assess$sigmaPriorCode
  statTable$sigUPriorCode   <- blob$assess$sigUPriorCode
  statTable$tauqPriorCode   <- blob$assess$tauqPriorCode


  # Now errors
  for (s in 1:nS)
  {
    statTable[s,"ssBnT"]    <- median (abs(ss$err.mle$BnT[success,s]), na.rm=TRUE)
    statTable[s,"msBnT"]    <- median (abs(ms$err.mle$BnT[success,s]), na.rm=TRUE)
    statTable[s,"ssUmsy"]   <- median (abs(ss$err.mle$Umsy[success,s]), na.rm=TRUE)
    statTable[s,"msUmsy"]   <- median (abs(ms$err.mle$Umsy[success,s]), na.rm=TRUE)
    statTable[s,"ssBmsy"]   <- median (abs(ss$err.mle$Bmsy[success,s]), na.rm=TRUE)
    statTable[s,"msBmsy"]   <- median (abs(ms$err.mle$Bmsy[success,s]), na.rm=TRUE)
    statTable[s,"ssDep"]    <- median (abs(ss$err.mle$dep[success,s]), na.rm=TRUE)
    statTable[s,"msDep"]    <- median (abs(ms$err.mle$dep[success,s]), na.rm=TRUE)
    statTable[s,"ssq"]      <- median (abs(ss$err.mle$q[success,s]), na.rm=TRUE)
    statTable[s,"msq"]      <- median (abs(ms$err.mle$q[success,s]), na.rm=TRUE)
    statTable[s,"ssHessPD"] <- sum ( ss$hesspd[,s] , na.rm=TRUE)
    statTable[s,"msHessPD"] <- sum ( ms$hesspd , na.rm=TRUE)
  }
  
  # return
  statTable
}

# .simStatMRE()
# Produces a statistics table for leading pars in a simulation from
# the output produced by a runSimEst() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame of mean squared error BnT and Umsy
# usage:    in lapply to produce stats for a group of simulations
.simStatMRE <- function ( sim=1 )
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
  species <- blob$ctrl$speciesName[1:nS]
  nReps   <- blob$ctrl$nReps

  # get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  success <- blob$goodReps

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","species","kappaTrue",
                  "SigmaTrue", "kappaMult", "corr","ssBnT","msBnT","ssUmsy","msUmsy",
                  "ssBmsy","msBmsy","ssDep","msDep",
                  "ssq", "msq", "msq025", "msq975", "ssq025", "ssq975",
                  "msHessPD", "ssHessPD","nReps",
                  "Umax", "tUpeak", "tUtrough", 
                  "tau2OM","nS","lastNegCorr",
                  "UmsyOM", "BmsyOM", "qOM",
                  "qbar025", "qbar975", "tauq2025", "tauq2975" )
  
  statTable <- matrix(NA,nrow=nS,ncol=length(colLabels))
  

  colnames(statTable)   <- colLabels
  statTable             <- as.data.frame(statTable)

  # get multiplier for shared effects
  kappaMult <- opMod$kappaMult
  if (is.null(opMod$kappaMult)) kappaMult <- 1

  # Start filling stat table
  # First, info and true pars
  statTable$scenario        <- blob$ctrl$scenarioName
  statTable$mp              <- blob$ctrl$mpLabel
  statTable$species         <- species[1:nS]
  statTable$kappaTrue       <- ifelse(is.null(opMod$kappa2),sqrt(opMod$pars$kappa2),sqrt(opMod$kappa2))*kappaMult
  statTable$SigmaTrue       <- ifelse(is.null(opMod$SigmaDiag),sqrt(opMod$pars$Sigma2[1:nS]),sqrt(opMod$SigmaDiag[1:nS]))
  statTable$kappaMult       <- kappaMult
  statTable$corr            <- ifelse(is.null(opMod$corrOffDiag),opMod$corrMult,opMod$corrOffDiag)
  statTable$nReps           <- length(success)
  statTable$Umax            <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUtrough        <- opMod$tUtrough
  statTable$tUpeak          <- opMod$tUpeak
  statTable$tau2OM          <- opMod$tau2[1:nS]
  statTable$lastNegCorr     <- ifelse(is.null(opMod$lastNegCorr),ifelse(nS==5,TRUE,FALSE),opMod$lastNegCorr)
  statTable$nS              <- opMod$nS
  statTable$UmsyOM          <- opMod$Umsy[1:nS]
  statTable$BmsyOM          <- opMod$Bmsy[1:nS]
  statTable$qOM             <- opMod$q[1:nS]
  statTable$s2q             <- blob$assess$s2lnq
  statTable$fixqbar         <- blob$assess$fixqbar
  statTable$fixtauq2        <- blob$assess$fixqbar
  statTable$s2Umsy          <- blob$assess$s2lnUmsy
  if (is.null(blob$assess$tauq2Prior))
  {
    statTable$tauq2P1        <- ifelse(is.null(blob$assess$tauq2IGa),blob$assess$tauq2IG[1],blob$assess$tauq2IGa)
    statTable$tauq2P2        <- ifelse(is.null(blob$assess$tauq2IGb),blob$assess$tauq2IG[2],blob$assess$tauq2IGb)
    statTable$sigU2P1        <- ifelse(is.null(blob$assess$sigU2IGa),blob$assess$sigU2IG[1],blob$assess$sigU2IGa)
    statTable$sigU2P2        <- ifelse(is.null(blob$assess$sigU2IGb),blob$assess$sigU2IG[2],blob$assess$sigU2IGb)  
  } else {
    statTable$tauq2P1         <- ifelse(is.null(blob$assess$tauq2P1),blob$assess$tauq2Prior[1],blob$assess$tauq2P1)
    statTable$tauq2P2         <- ifelse(is.null(blob$assess$tauq2P2),blob$assess$tauq2Prior[2],blob$assess$tauq2P2)
    statTable$sigU2P1         <- ifelse(is.null(blob$assess$sigU2P1),blob$assess$sigU2Prior[1],blob$assess$sigU2P1)
    statTable$sigU2P2         <- ifelse(is.null(blob$assess$sigU2P2),blob$assess$sigU2Prior[2],blob$assess$sigU2P2) 
  }
  statTable$sigmaPriorCode  <- blob$assess$sigmaPriorCode
  statTable$sigUPriorCode   <- blob$assess$sigUPriorCode
  statTable$tauqPriorCode   <- blob$assess$tauqPriorCode



  # Now errors
  for (s in 1:nS)
  {
    statTable[s,"ssBnT"]    <- median (ss$err.mle$BnT[success,s], na.rm = TRUE )
    statTable[s,"msBnT"]    <- median (ms$err.mle$BnT[success,s], na.rm = TRUE )
    statTable[s,"ssUmsy"]   <- median (ss$err.mle$Umsy[success,s], na.rm = TRUE )
    statTable[s,"ssUmsy025"]<- quantile (ss$err.mle$Umsy[success,s], probs = 0.025, na.rm = TRUE )
    statTable[s,"ssUmsy975"]<- quantile (ss$err.mle$Umsy[success,s], probs = 0.975, na.rm = TRUE )
    statTable[s,"msUmsy"]   <- median (ms$err.mle$Umsy[success,s], na.rm = TRUE )
    statTable[s,"msUmsy025"]<- quantile (ms$err.mle$Umsy[success,s], probs = 0.025, na.rm = TRUE )
    statTable[s,"msUmsy975"]<- quantile (ms$err.mle$Umsy[success,s], probs = 0.975, na.rm = TRUE )
    statTable[s,"ssBmsy"]   <- median (ss$err.mle$Bmsy[success,s], na.rm = TRUE )
    statTable[s,"msBmsy"]   <- median (ms$err.mle$Bmsy[success,s], na.rm = TRUE )
    statTable[s,"ssDep"]    <- median (ss$err.mle$dep[success,s], na.rm = TRUE )
    statTable[s,"msDep"]    <- median (ms$err.mle$dep[success,s], na.rm = TRUE )
    statTable[s,"ssq"]      <- median (ss$err.mle$q[success,s], na.rm = TRUE )
    statTable[s,"ssq025"]   <- quantile (ss$err.mle$q[success,s], probs = 0.025, na.rm = TRUE )
    statTable[s,"ssq975"]   <- quantile (ss$err.mle$q[success,s], probs = 0.975, na.rm = TRUE )
    statTable[s,"msq"]      <- median (ms$err.mle$q[success,s], na.rm = TRUE )
    statTable[s,"msq025"]   <- quantile (ms$err.mle$q[success,s], probs = 0.025, na.rm = TRUE )
    statTable[s,"msq975"]   <- quantile (ms$err.mle$q[success,s], probs = 0.975, na.rm = TRUE )
    statTable[s,"ssHessPD"] <- sum ( ss$hesspd[,s] , na.rm = TRUE )
    statTable[s,"msHessPD"] <- sum ( ms$hesspd , na.rm = TRUE )
  }
  statTable$qbar     <- median( blob$am$ms$qbar[success], na.rm=TRUE)
  statTable$qbar025  <- quantile( blob$am$ms$qbar[success], probs = 0.025, na.rm = TRUE )
  statTable$qbar975  <- quantile( blob$am$ms$qbar[success], probs = 0.975, na.rm = TRUE )
  statTable$tauq2    <- median( blob$am$ms$tauq2[success], na.rm=TRUE)
  statTable$tauq2025 <- quantile( blob$am$ms$tauq2[success], probs = 0.025, na.rm = TRUE )
  statTable$tauq2975 <- quantile( blob$am$ms$tauq2[success], probs = 0.975, na.rm = TRUE )
  statTable$sigU2    <- median( blob$am$ms$sigU2[success], na.rm=TRUE)
  statTable$Umsybar  <- median( blob$am$ms$Umsybar[success],na.rm=TRUE)
  statTable$sigU2025 <- quantile( blob$am$ms$sigU2[success], probs = 0.025, na.rm = TRUE )
  statTable$sigU2975 <- quantile( blob$am$ms$sigU2[success], probs = 0.975, na.rm = TRUE )
  statTable$Ubar025  <- quantile( blob$am$ms$Umsybar[success], probs = 0.025, na.rm = TRUE )
  statTable$Ubar975  <- quantile( blob$am$ms$Umsybar[success], probs = 0.975, na.rm = TRUE )
  # return
  statTable
}

# .simStatMRE()
# Produces a statistics table for leading pars in a simulation from
# the output produced by a runSimEst() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame of mean squared error BnT and Umsy
# usage:    in lapply to produce stats for a group of simulations
.simStatRE <- function ( sim=1 )
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
  species <- blob$ctrl$speciesName[1:nS]
  nReps   <- blob$ctrl$nReps

  # get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  success <- blob$goodReps

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","species","kappaTrue",
                  "SigmaTrue", "kappaMult", "corr","ssBnT","msBnT","ssUmsy","msUmsy",
                  "ssBmsy","msBmsy","ssDep","msDep",
                  "ssq", "msq",
                  "msHessPD", "ssHessPD","nReps",
                  "Umax", "tUpeak", "tUtrough", 
                  "tau2OM","nS","lastNegCorr",
                  "UmsyOM", "BmsyOM", "qOM", "rep")
  
  statTable <- matrix(NA,nrow=nS*nReps,ncol=length(colLabels))
  

  colnames(statTable)   <- colLabels
  statTable             <- as.data.frame(statTable)

  # get multiplier for shared effects
  kappaMult <- opMod$kappaMult
  if (is.null(opMod$kappaMult)) kappaMult <- 1

  # Start filling stat table
  # First, info and true pars
  statTable$scenario        <- blob$ctrl$scenarioName
  statTable$mp              <- blob$ctrl$mpLabel
  statTable$species         <- species[1:nS]
  statTable$kappaTrue       <- ifelse(is.null(opMod$kappa2),sqrt(opMod$pars$kappa2),sqrt(opMod$kappa2))*kappaMult
  statTable$SigmaTrue       <- ifelse(is.null(opMod$SigmaDiag),sqrt(opMod$pars$Sigma2[1:nS]),sqrt(opMod$SigmaDiag[1:nS]))
  statTable$kappaMult       <- kappaMult
  statTable$corr            <- ifelse(is.null(opMod$corrOffDiag),opMod$corrMult,opMod$corrOffDiag)
  statTable$nReps           <- length(success)
  statTable$Umax            <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUtrough        <- opMod$tUtrough
  statTable$tUpeak          <- opMod$tUpeak
  statTable$tau2OM          <- opMod$tau2[1:nS]
  statTable$lastNegCorr     <- ifelse(is.null(opMod$lastNegCorr),ifelse(nS==5,TRUE,FALSE),opMod$lastNegCorr)
  statTable$nS              <- opMod$nS
  statTable$UmsyOM          <- opMod$Umsy[1:nS]
  statTable$BmsyOM          <- opMod$Bmsy[1:nS]
  statTable$qOM             <- opMod$q[1:nS]
  statTable$s2q             <- blob$assess$s2lnq
  statTable$fixqbar         <- blob$assess$fixqbar
  statTable$fixtauq2        <- blob$assess$fixqbar
  statTable$s2Umsy          <- blob$assess$s2lnUmsy
  if (is.null(blob$assess$tauq2Prior))
  {
    statTable$tauq2P1        <- ifelse(is.null(blob$assess$tauq2IGa),blob$assess$tauq2IG[1],blob$assess$tauq2IGa)
    statTable$tauq2P2        <- ifelse(is.null(blob$assess$tauq2IGb),blob$assess$tauq2IG[2],blob$assess$tauq2IGb)
    statTable$sigU2P1        <- ifelse(is.null(blob$assess$sigU2IGa),blob$assess$sigU2IG[1],blob$assess$sigU2IGa)
    statTable$sigU2P2        <- ifelse(is.null(blob$assess$sigU2IGb),blob$assess$sigU2IG[2],blob$assess$sigU2IGb)  
  } else {
    statTable$tauq2P1         <- ifelse(is.null(blob$assess$tauq2P1),blob$assess$tauq2Prior[1],blob$assess$tauq2P1)
    statTable$tauq2P2         <- ifelse(is.null(blob$assess$tauq2P2),blob$assess$tauq2Prior[2],blob$assess$tauq2P2)
    statTable$sigU2P1         <- ifelse(is.null(blob$assess$sigU2P1),blob$assess$sigU2Prior[1],blob$assess$sigU2P1)
    statTable$sigU2P2         <- ifelse(is.null(blob$assess$sigU2P2),blob$assess$sigU2Prior[2],blob$assess$sigU2P2) 
  }
  statTable$sigmaPriorCode  <- blob$assess$sigmaPriorCode
  statTable$sigUPriorCode   <- blob$assess$sigUPriorCode
  statTable$tauqPriorCode   <- blob$assess$tauqPriorCode

  # Now errors
  for (r in 1:nReps)
  {
    statTable[(r-1)*nS+(1:nS),"rep"]      <- r
    statTable[(r-1)*nS+(1:nS),"ssBnT"]    <- (ss$err.mle$BnT[r,] )
    statTable[(r-1)*nS+(1:nS),"msBnT"]    <- (ms$err.mle$BnT[r,] )
    statTable[(r-1)*nS+(1:nS),"ssUmsy"]   <- (ss$err.mle$Umsy[r,] )
    statTable[(r-1)*nS+(1:nS),"msUmsy"]   <- (ms$err.mle$Umsy[r,] )
    statTable[(r-1)*nS+(1:nS),"ssBmsy"]   <- (ss$err.mle$Bmsy[r,] )
    statTable[(r-1)*nS+(1:nS),"msBmsy"]   <- (ms$err.mle$Bmsy[r,] )
    statTable[(r-1)*nS+(1:nS),"ssDep"]    <- (ss$err.mle$dep[r,] )
    statTable[(r-1)*nS+(1:nS),"msDep"]    <- (ms$err.mle$dep[r,] )
    statTable[(r-1)*nS+(1:nS),"ssq"]      <- (ss$err.mle$q[r,] )
    statTable[(r-1)*nS+(1:nS),"msq"]      <- (ms$err.mle$q[r,] )
    statTable[(r-1)*nS+(1:nS),"ssHessPD"] <- (ss$hesspd[r,] )
    statTable[(r-1)*nS+(1:nS),"msHessPD"] <- ms$hesspd[r]
    statTable[(r-1)*nS+(1:nS),"qbar"]     <- ( blob$am$ms$qbar[r])
    statTable[(r-1)*nS+(1:nS),"tauq2"]    <- ( blob$am$ms$tauq2[r])
    statTable[(r-1)*nS+(1:nS),"sigU2"]    <- ( blob$am$ms$sigU2[r])
    statTable[(r-1)*nS+(1:nS),"Umsybar"]  <- ( blob$am$ms$Umsybar[r])
  }
  
  # return
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
  species <- blob$ctrl$speciesName[1:nS]
  nReps   <- blob$ctrl$nReps

  # get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  success <- blob$goodReps

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","species","kappaTrue",
                  "SigmaTrue", "kappaMult","ssBnT","msBnT","ssUmsy","msUmsy",
                  "ssBmsy","msBmsy","ssMSY","msMSY","ssDep","msDep",
                  "ssq","msq", "msHessPD", "ssHessPD", "nReps",
                  "Umax", "tUpeak", "tUtrough", "tau2OM" )
  
  statTable <- matrix(NA,nrow=nS,ncol=length(colLabels))
  
  colnames(statTable)   <- colLabels
  statTable             <- as.data.frame(statTable)

  # get multiplier for shared effects
  kappaMult <- opMod$kappaMult
  if (is.null(opMod$kappaMult)) kappaMult <- 1

  # Start filling stat table
  # First, info and true pars
  statTable$scenario        <- blob$ctrl$scenarioName
  statTable$mp              <- blob$ctrl$mpLabel
  statTable$species         <- species[1:nS]
  statTable$kappaTrue       <- ifelse(is.null(opMod$kappa2),sqrt(opMod$pars$kappa2),sqrt(opMod$kappa2))*kappaMult
  statTable$SigmaTrue       <- ifelse(is.null(opMod$SigmaDiag),sqrt(opMod$pars$Sigma2[1:nS]),sqrt(opMod$SigmaDiag[1:nS]))
  statTable$kappaMult       <- kappaMult
  statTable$corr            <- ifelse(is.null(opMod$corrOffDiag),opMod$corrMult,opMod$corrOffDiag)
  statTable$nReps           <- length(success)
  statTable$Umax            <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUtrough        <- opMod$tUtrough
  statTable$tUpeak          <- opMod$tUpeak
  statTable$tau2OM          <- opMod$tau2[1:nS]
  statTable$lastNegCorr     <- ifelse(is.null(opMod$lastNegCorr),ifelse(nS==5,TRUE,FALSE),opMod$lastNegCorr)
  statTable$nS              <- opMod$nS
  statTable$UmsyOM          <- opMod$Umsy[1:nS]
  statTable$BmsyOM          <- opMod$Bmsy[1:nS]
  statTable$qOM             <- opMod$q[1:nS]
  statTable$s2q             <- blob$assess$s2lnq
  statTable$fixqbar         <- blob$assess$fixqbar
  statTable$fixtauq2        <- blob$assess$fixqbar
  statTable$s2Umsy          <- blob$assess$s2lnUmsy
  if (is.null(blob$assess$tauq2Prior))
  {
    statTable$tauq2P1        <- ifelse(is.null(blob$assess$tauq2IGa),blob$assess$tauq2IG[1],blob$assess$tauq2IGa)
    statTable$tauq2P2        <- ifelse(is.null(blob$assess$tauq2IGb),blob$assess$tauq2IG[2],blob$assess$tauq2IGb)
    statTable$sigU2P1        <- ifelse(is.null(blob$assess$sigU2IGa),blob$assess$sigU2IG[1],blob$assess$sigU2IGa)
    statTable$sigU2P2        <- ifelse(is.null(blob$assess$sigU2IGb),blob$assess$sigU2IG[2],blob$assess$sigU2IGb)  
  } else {
    statTable$tauq2P1         <- ifelse(is.null(blob$assess$tauq2P1),blob$assess$tauq2Prior[1],blob$assess$tauq2P1)
    statTable$tauq2P2         <- ifelse(is.null(blob$assess$tauq2P2),blob$assess$tauq2Prior[2],blob$assess$tauq2P2)
    statTable$sigU2P1         <- ifelse(is.null(blob$assess$sigU2P1),blob$assess$sigU2Prior[1],blob$assess$sigU2P1)
    statTable$sigU2P2         <- ifelse(is.null(blob$assess$sigU2P2),blob$assess$sigU2Prior[2],blob$assess$sigU2P2) 
  }
  statTable$sigmaPriorCode  <- blob$assess$sigmaPriorCode
  statTable$sigUPriorCode   <- blob$assess$sigUPriorCode
  statTable$tauqPriorCode   <- blob$assess$tauqPriorCode


  # Now errors
  for (s in 1:nS)
  {
    statTable[s,"ssBnT"]    <- mean ( (ss$Bt[success,s,nT] - om$Bt[success,s,nT])^2, na.rm=TRUE)
    statTable[s,"msBnT"]    <- mean ( (ms$Bt[success,s,nT] - om$Bt[success,s,nT])^2, na.rm=TRUE)
    statTable[s,"ssUmsy"]   <- mean ( (ss$Umsy[success,s] - opMod$Umsy[s])^2, na.rm=TRUE)
    statTable[s,"msUmsy"]   <- mean ( (ms$Umsy[success,s] - opMod$Umsy[s])^2, na.rm=TRUE)
    statTable[s,"ssBmsy"]   <- mean ( (ss$Bmsy[success,s] - opMod$Bmsy[s])^2, na.rm=TRUE)
    statTable[s,"msBmsy"]   <- mean ( (ms$Bmsy[success,s] - opMod$Bmsy[s])^2, na.rm=TRUE)
    statTable[s,"ssMSY"]    <- mean ( (ss$msy[success,s] -  opMod$Umsy[s]*opMod$Bmsy[s])^2, na.rm=TRUE)
    statTable[s,"msMSY"]    <- mean ( (ms$msy[success,s] -  opMod$Umsy[s]*opMod$Bmsy[s])^2, na.rm=TRUE)
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
  sYear <- opMod$sYear
  fYear <- opMod$fYear
  nT    <- fYear - min(sYear) + 1

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
    # browser()
    ss$err.mle$Bmsy[success,s]   <- (ss$Bmsy[success,s] - opMod$Bmsy[s])/opMod$Bmsy[s]
    ss$err.mle$Umsy[success,s]   <- (ss$Umsy[success,s] - opMod$Umsy[s])/opMod$Umsy[s]
    ss$err.mle$kappa2[success,s] <- (ss$kappa2[success,s] - (opMod$kappa2*(opMod$kappaMult^2)+opMod$SigmaDiag[s]))/(opMod$kappa2*(opMod$kappaMult^2)+opMod$SigmaDiag[s])
    ss$err.mle$tau2[success,s]   <- (ss$tau2[success,s] - opMod$tau2[s])/opMod$tau2[s]
    ss$err.mle$q[success,s]      <- (ss$q[success,s] - opMod$q[s])/opMod$q[s]
    ss$err.mle$dep[success,s]    <- (ss$dep[success,s,nT] - om$dep[success,s,nT])/om$dep[success,s,nT]
    ss$err.mle$BnT[success,s]    <- (ss$Bt[success,s,nT] - om$Bt[success,s,nT])/om$Bt[success,s,nT]
    ss$err.mle$totRE[success,s]  <- ss$err.mle$kappa2[success,s]

    # Now fill in ms MLE relative errors
    # some are only estimated once (instead of nS times)
    if (s == 1)
    {
      ms$err.mle$kappa2[success,]    <- (ms$kappa2[success,] - opMod$kappa2*(opMod$kappaMult^2))/opMod$kappa2/(opMod$kappaMult^2)
    }
    # Now the rest of the pars
    ms$err.mle$Sigma2[success,s] <- (ms$Sigma2[success,s] - opMod$SigmaDiag[s])/opMod$SigmaDiag[s]
    ms$err.mle$tau2[success,s]   <- (ms$tau2[success,s] - opMod$tau2[s])/opMod$tau2[s]
    ms$err.mle$Bmsy[success,s]   <- (ms$Bmsy[success,s] - opMod$Bmsy[s])/opMod$Bmsy[s]
    ms$err.mle$Umsy[success,s]   <- (ms$Umsy[success,s] - opMod$Umsy[s])/opMod$Umsy[s]    
    ms$err.mle$q[success,s]      <- (ms$q[success,s] - opMod$q[s])/opMod$q[s]
    ms$err.mle$dep[success,s]    <- (ms$dep[success,s,nT] - om$dep[success,s,nT])/om$dep[success,s,nT]
    ms$err.mle$BnT[success,s]    <- (ms$Bt[success,s,nT] - om$Bt[success,s,nT])/om$Bt[success,s,nT]

    # calculate total RE variance
    totVarFit <- ms$kappa2 + ms$Sigma2[,s]
    totVarOM  <- opMod$kappa2*(opMod$kappaMult^2) + opMod$SigmaDiag[s]
    ms$err.mle$totRE[success,s]  <- (totVarFit[success] - totVarOM) / totVarOM
  }

  # Append these to blob
  blob$am$ss <- ss
  blob$am$ms <- ms

  # Now save the good replicates
  blob$goodReps <- success

  blob
}



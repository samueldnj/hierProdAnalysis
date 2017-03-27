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
.DASEmodelSelection <- function(  tabName = "rKqExp_MRE.csv",
                                  resp = "q",
                                  spec = "Dover",
                                  expVars = c("qOM","UmsyOM","BmsyOM"),
                                  sig = 0.05, maxDeltaAICc = 20 )
{
  # First, fit main effects
  mainEff <- .DASEexperiment( tabName = tabName,
                              resp = resp, spec = spec,
                              scaled = FALSE,
                              baseRHS = NULL,
                              newExp = expVars )

  # Then, look at the main effects summary
  mainSummSS <- lapply( X = mainEff$ssGLM, FUN = summary )
  mainSummMS <- lapply( X = mainEff$msGLM, FUN = summary )

  coeffSS <- coef(mainSummSS[[length(mainSummSS)]])
  coeffMS <- coef(mainSummSS[[length(mainSummSS)]])  

  ssInsig <- which( coeffSS[,4] >= sig )
  msInsig <- which( coeffMS[,4] >= sig )

  dropSS <- names(ssInsig)
  dropMS <- names(msInsig)

  signifSS <- expVars[ !(expVars %in% dropSS) ]
  signifMS <- expVars[ !(expVars %in% dropMS) ]

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
  # make interactions and second order effects
  intSS <- makeIntEff ( signifSS )
  soSS  <- paste( "I(", signifSS, "^2)", sep = "" )

  # Now paste main effects together
  signifSSrhs <- paste(signifSS, collapse = " + ")
  # And run the experiment again
  fullSS <- .DASEexperiment(  tabName = tabName, resp = resp,
                              spec = spec, scaled = FALSE,
                              baseRHS = signifSSrhs,
                              newExp = c(intSS,soSS) )

  # Run for fullMS if different main effects
  if( length(setdiff(signifMS,signifSS)) != 0 )
  {
    signifMSrhs <- paste(signifMS, collapse = " + ")
    # make interactions and second order effects
    intMS <- makeIntEff ( signifMS )
    soMS  <- paste( "I(", signifMS, "^2)", sep = "" )


    fullMS <- .DASEexperiment(  tabName = tabName, resp = resp,
                              spec = spec, scaled = FALSE,
                              baseRHS = signifMSrhs,
                              newExp = c(intMS,soMS) )
  } else
    fullMS <- fullSS
  
  fullSummSS <- lapply( X = fullSS$ssGLM, FUN = summary )
  fullSummMS <- lapply( X = fullMS$msGLM, FUN = summary )

  # Now rank AICc values in each model
  aicRankSS <- AICrank( fullSummSS )
  aicRankMS <- AICrank( fullSummMS )

  # Now pull the models that fit the AIC and significance criteria
  modelNumSS  <- aicRankSS[ which( (  aicRankSS$maxPr < sig) &
                                      aicRankSS$deltaAICc < maxDeltaAICc ), "Number" ]
  modelNumMS  <- aicRankSS[ which( (  aicRankSS$maxPr < sig) &
                                      aicRankSS$deltaAICc < maxDeltaAICc ), "Number" ]

  ms <- list( sel = fullMS$msGLM[modelNumMS],
              aicRank = aicRankMS,
              fullMS = fullMS 
            )

  ss <- list( sel = fullMS$ssGLM[modelNumSS],
              aicRank = aicRankSS,
              fullSS = fullSS 
            )

  out <- list( ms = ms, ss = ss )

  out
}


# Automate checking for AIC values
AICrank <- function ( summList )
{
  # Count models
  nModels <- length( summList )

  # Create a data frame to hold the info
  rank.df <- matrix(NA, ncol = 4, nrow = nModels )
  colnames(rank.df) <- c( "Number", "AICc", "deltaAICc", "maxPr" )
  rank.df <- as.data.frame(rank.df)

  rank.df[,1] <- 1:nModels
  for(k in 1:nModels)
  {
    nPar <- nrow(summList[[k]]$coefficients)
    nObs <- length(summList[[k]]$deviance.resid)
    rank.df[k,"AICc"] <- summList[[k]]$aic + nPar*(nPar + 1) / (nObs - nPar - 1)
    rank.df[k,"maxPr"] <- max(coef(summList[[k]])[,4])
  } 

  rank.df[,"deltaAICc"] <- rank.df[,"AICc"] - min(rank.df[,"AICc"])
  rank.df <- rank.df[order(rank.df[,"deltaAICc"]),]
  
  rank.df
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
.DASEexperiment <- function(  tabName = "lowK_rKq_RE.csv",
                              resp = "q",
                              spec = "Dover",
                              scaled = TRUE,
                              baseRHS = NULL,
                              newExp = c( "qOM", "UmsyOM", "BmsyOM" ) )
{
  tabPath <- file.path(getwd(),"project/statistics",tabName)
  table <- read.csv( tabPath, header=TRUE, stringsAsFactors=FALSE ) 

  # restrict to correct number of species
  table  <-   table %>%
              filter( species == spec,
                      BmsyOM < 100, BmsyOM > 5 )
  # Calculate a table with median bias
  medTab <- table %>%
            group_by( qOM, UmsyOM, BmsyOM ) %>%
            summarise(  nReps = sum(!is.na(ssq)),
                        varssq = var(ssq,na.rm =T),
                        varmsq = var(msq,na.rm =T),
                        varssUmsy = var(ssUmsy,na.rm =T),
                        varmsUmsy = var(msUmsy,na.rm =T),
                        varssBmsy = var(ssBmsy,na.rm =T),
                        varmsBmsy = var(msBmsy,na.rm =T), 
                        ssq = median(ssq,na.rm =T),
                        msq = median(msq,na.rm =T),
                        ssUmsy = median(ssUmsy,na.rm =T),
                        msUmsy = median(msUmsy,na.rm =T),
                        ssBmsy = median(ssBmsy,na.rm =T),
                        msBmsy = median(msBmsy,na.rm =T)
                      )

  # medTab <- file.path(getwd(),"project/statistics",medTab)
  # medTab <- read.csv( medTab, header=TRUE, stringsAsFactors=FALSE )   
  # medTab <- medTab %>% filter( species == spec )


  # Scale inputs
  if( scaled )
  {
    # Calculate ranges of inputs
    qRange      <- range(table$qOM)
    URange      <- range(table$UmsyOM)
    BRange      <- range(table$BmsyOM)
    tauRange    <- range(table$tau2OM)
    kappaRange  <- range(table$kappaTrue)
    corrRange   <- range(table$corr)
    SigmaRange  <- range(table$SigmaTrue)

    # Now calculate gradients of inputs
    qGrad       <- 2 / (qRange[2] - qRange[1])
    UGrad       <- 2 / (URange[2] - URange[1])
    BGrad       <- 2 / (BRange[2] - BRange[1])
    tauGrad     <- 2 / (tauRange[2] - tauRange[1])
    kappaGrad   <- 2 / (kappaRange[2] - kappaRange[1])
    corrGrad    <- 2 / (corrRange[2] - corrRange[1])
    SigmaGrad   <- 2 / (SigmaRange[2] - SigmaRange[1])

    table <-  table %>%
              mutate( qOM = qGrad*qOM + (1 - qGrad*qRange[2]),
                      UmsyOM = UGrad*UmsyOM + (1 - UGrad*URange[2]),
                      BmsyOM = BGrad*BmsyOM + (1 - BGrad*BRange[2]),
                      tau2OM = tauGrad*tau2OM + (1 - tauGrad*tauRange[2]),
                      kappaTrue = kappaGrad*kappaTrue + (1 - kappaGrad*kappaRange[2]),
                      corr = corrGrad*corr + (1 - corrGrad*corrRange[2]),
                      SigmaTrue = SigmaGrad*SigmaTrue + (1 - SigmaGrad*SigmaRange[2]))  
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

  # Create LHS of LM formula object
  ssResp <- paste( "ss", resp, sep = "" )
  msResp <- paste( "ms", resp, sep = "" )

  # Generate the sets of possible effects
  # First, take the power set of main effects and collapse into formulae
  newExpCombos  <- powerset( newExp )
  newExpCombos  <- lapply( X = newExpCombos, FUN = paste, collapse = " + " )
  rhsList       <- paste( baseRHS, newExpCombos, sep = " + " )
  # Fit to all possible main effects for both models
  ssGLM       <- lapply( X = rhsList, FUN = fitGLM, resp = ssResp, dat = table )
  msGLM       <- lapply( X = rhsList, FUN = fitGLM, resp = msResp, dat = table )

  # A function to lapply to the glmObjects and compute the AIC
  medResids   <- function( glmObj, data, resp )
  {
    pred <- predict( glmObj, newdata = data )
    resids <- data[,resp] - pred
    resids
  }

  # Now run predict the medTab from the models
  ssResids    <- lapply( X = ssGLM, FUN = medResids, data = medTab, resp = ssResp )
  msResids    <- lapply( X = msGLM, FUN = medResids, data = medTab, resp = msResp )


  return( list (  ssGLM = ssGLM, ssResids = ssResids,
                  msGLM = msGLM, msResids = msResids,
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
  .statTableMSE(sims,paste(tabNameRoot,"_MSE.csv",sep=""))
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
  statTable$species         <- species
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
  statTable$s2q             <- blob$assess$s2lnq
  statTable$tauq2           <- median(blob$am$ms$tauq2,na.rm=TRUE)
  statTable$qbar            <- median(blob$am$ms$qbar,na.rm=TRUE)
  statTable$s2Umsy          <- blob$assess$s2lnUmsy
  statTable$sigU2           <- median(blob$am$ms$sigU2,na.rm=TRUE)
  statTable$Umsybar         <- median(blob$am$ms$Umsybar,na.rm=TRUE)
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
  statTable$species         <- species
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
  statTable$UmsyOM          <- opMod$Umsy
  statTable$BmsyOM          <- opMod$Bmsy
  statTable$qOM             <- opMod$q
  statTable$s2q             <- blob$assess$s2lnq
  statTable$fixqbar         <- blob$assess$fixqbar
  statTable$fixtauq2        <- blob$assess$fixqbar
  statTable$qOM             <- blob$opMod$q[1:nS]
  statTable$UmsyOM          <- blob$opMod$Umsy[1:nS]
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
  statTable$species         <- species
  statTable$kappaTrue       <- ifelse(is.null(opMod$kappa2),sqrt(opMod$pars$kappa2),sqrt(opMod$kappa2))*kappaMult
  statTable$SigmaTrue       <- ifelse(is.null(opMod$SigmaDiag),sqrt(opMod$pars$Sigma2[1:nS]),sqrt(opMod$SigmaDiag[1:nS]))
  statTable$kappaMult       <- kappaMult
  statTable$corr            <- ifelse(is.null(opMod$corrOffDiag),opMod$corrMult,opMod$corrOffDiag)
  statTable$nReps           <- nReps
  statTable$Umax            <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUtrough        <- opMod$tUtrough
  statTable$tUpeak          <- opMod$tUpeak
  statTable$tau2OM          <- opMod$tau2[1:nS]
  statTable$lastNegCorr     <- ifelse(is.null(opMod$lastNegCorr),ifelse(nS==5,TRUE,FALSE),opMod$lastNegCorr)
  statTable$nS              <- opMod$nS
  statTable$UmsyOM          <- opMod$Umsy
  statTable$BmsyOM          <- opMod$Bmsy
  statTable$qOM             <- opMod$q
  statTable$s2q             <- blob$assess$s2lnq
  statTable$fixqbar         <- blob$assess$fixqbar
  statTable$fixtauq2        <- blob$assess$fixqbar
  statTable$qOM             <- blob$opMod$q[1:nS]
  statTable$UmsyOM          <- blob$opMod$Umsy[1:nS]
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
  statTable$species         <- species
  statTable$kappaTrue       <- ifelse(is.null(opMod$kappa2),sqrt(opMod$pars$kappa2),sqrt(opMod$kappa2))*kappaMult
  statTable$SigmaTrue       <- ifelse(is.null(opMod$SigmaDiag),sqrt(opMod$pars$Sigma2[1:nS]),sqrt(opMod$SigmaDiag[1:nS]))
  statTable$kappaMult       <- kappaMult
  statTable$corr            <- ifelse(is.null(opMod$corrOffDiag),opMod$corrMult,opMod$corrOffDiag)
  statTable$nReps           <- nReps
  statTable$Umax            <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUtrough        <- opMod$tUtrough
  statTable$tUpeak          <- opMod$tUpeak
  statTable$tau2OM          <- opMod$tau2[1:nS]
  statTable$lastNegCorr     <- ifelse(is.null(opMod$lastNegCorr),ifelse(nS==5,TRUE,FALSE),opMod$lastNegCorr)
  statTable$nS              <- opMod$nS
  statTable$s2q             <- blob$assess$s2lnq
  statTable$tauq2           <- median(blob$am$ms$tauq2,na.rm=TRUE)
  statTable$qbar            <- median(blob$am$ms$qbar,na.rm=TRUE)
  statTable$s2Umsy          <- blob$assess$s2lnUmsy
  statTable$sigU2           <- median(blob$am$ms$sigU2,na.rm=TRUE)
  statTable$Umsybar         <- median(blob$am$ms$Umsybar,na.rm=TRUE)
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



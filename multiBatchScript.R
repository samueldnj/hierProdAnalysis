# --------------------------------------------------------------------------
# multiBatchScript.R
# 
# Runs a list of batch files in turn.
# 
# Author: Samuel DN Johnson
# Date: 27 March, 2017
# 
# --------------------------------------------------------------------------

# Source the simulation framework
source("control.R")

# List batch file names, base control file names and experiment
# prefix names (for folder naming)
# vectors must length k.
batchFiles    <- c("infoScenarios.bch","infoScenarios.bch","pubBase4MPs.bch")
baseCtlFiles  <- c("simCtlFileRP.txt","simCtlFile.txt","simCtlFile.txt")
expPrefix     <- c("allSame_infoScenarios_RP","allSame_infoScenarios","pubBase")
plots         <- c(TRUE,TRUE,TRUE)

# Now loop over the experiments

for( i in 1:length(batchFiles))
{
  # Make batch design
  makeBatch(  batchCtlFile = batchFiles[i],
              baseCtlFile = baseCtlFiles[i] )
  # Run batch job in parallel
  .runBatchJob( par = T, prefix = expPrefix[i] )

  # Create stat tables
  sims <- grep(pattern = "sim", x = list.files("./project/"), value = T)
  nSims <- length(sims)
  .statTables(1:nSims,expPrefix[i],par=T)
  if(plots)
  {
    dumpPerfMetrics(  tabNameRoot = expPrefix[i], stockLabel = "Stock1",
                      vars = c("Umsy","BnT","Bmsy","Dep","q_1","q_2"),
                      varLabels = expression(U[MSY], B[T], B[MSY], B[T]/B[0], q[11], q[21]),
                      MPs = c("noJointPriors","qPriorOnly","UmsyPriorOnly","qUpriors" ),
                      MPlabels = expression("None", q, r, q/r ) ) 
    dumpBCsim(  simPath = file.path("./project"),
                prefix = expPrefix[i],
                MPs = c("noJointPriors","qPriorOnly","UmsyPriorOnly","qUpriors" ),
                MPlabels = expression("Single Stock","None", q, r, q/r ) )  
    dumpStockPerf(  simPath = file.path("./project"),
                    prefix = expPrefix[i],
                    MPs = c("noJointPriors","qPriorOnly","UmsyPriorOnly","qUpriors" ),
                    MPlabels = c( noJointPriors = "None", 
                                  qPriorOnly = expression(q), 
                                  UmsyPriorOnly = expression(r), 
                                  qUpriors = expression(q/r) ) )
  }
  

  # Now copy the project folder to dropbox
  copyDest <- file.path("/Volumes/home/thesisStuff/cwMSexperiments/TMB",paste(expPrefix[i],Sys.Date(),sep = "_") )
  dir.create( copyDest )

  # Copy project folder contents recursively to copyDest
  cat( "Copying project folder contents to ", copyDest, "\n", sep = "" )
  x <- file.copy( from = file.path( getwd(),"project"), to = copyDest, 
                  recursive = TRUE )

  if(!x)
  {
    cat(  "Error copying project folder contents to remote server,\n", 
          "using local dir instead" )
    copyDest <- file.path("../",paste(expPrefix[i],Sys.Date(),sep = "_") )
    dir.create( copyDest )

    # Copy project folder contents recursively to copyDest
    cat( "Copying project folder contents to ", copyDest, "\n", sep = "" )
    x <- file.copy( from = file.path( getwd(),"project"), to = copyDest, 
                    recursive = TRUE )
  }

  # Now remove the simulation folders and the contents of the batch sub-dir
  # from the project folder
  # sims <- grep(pattern = "sim", x = list.files("./project/"), value = T)
  simsPath <- file.path( getwd(),"project",sims)
  batchFldrContents <- list.files( file.path(getwd(), "project", "Batch") )
  batchContentsPath <- file.path( getwd(), "project", "Batch", batchFldrContents )

  # Copy out sims to dropbox, tidy up
  cat("Removing simulations from ./project/ \n", sep="")
  for(k in 1:length(simsPath))
    system(command=paste("rm -d -R ",simsPath[k],sep=""))

  cat("Removing batch files from ./project/batch folder\n", sep="")
  for(k in 1:length(batchContentsPath))
    system(command=paste("rm -d -R ",batchContentsPath[k],sep=""))

  cat("Experiment and tidy-up complete, ready to start next experiment.\n")
}


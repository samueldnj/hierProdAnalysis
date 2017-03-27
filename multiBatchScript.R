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

# Number of experiments:
nExp <- 2

# List batch file names, base control file names and experiment
# prefix names (for folder naming)
# vectors must length k.
batchFiles    <- c("msIncrements.bch", "msIncrementsNoneq.bch")
baseCtlFiles  <- c("simCtlFileAllSameBase.txt","simCtlFileAllSameBase.txt")
expPrefix     <- c("msIncrements", "msIncrements_noneq")

# Now loop over the experiments

for( i in 1:nExp)
{
  # Make batch design
  makeBatch(  batchCtlFile = batchFiles[i],
              baseCtlFile = baseCtlFiles[i] )
  # Run batch job in parallel
  .runBatchJob( par = T, prefix = expPrefix[i] )

  # Now copy the project folder to dropbox
  copyDest <- file.path("~/Dropbox/Work/cwMSexperiments",expPrefix[i])
  dir.create( copyDest )

  # Copy project folder contents recursively to copyDest
  cat( "Copying project folder contents to ", copyDest, "\n", sep = "" )
  file.copy(  from = file.path( getwd(),"/project"), to = copyDest, 
              recursive = TRUE )

  # Now remove the simulation folders and the contents of the batch sub-dir
  # from the project folder
  sims <- grep(pattern = "sim", x = prjContents)
  simsPath <- file.path( getwd(),"project",sims)
  batchFldrContents <- list.files( getwd(), "project", "batch" )
  batchContentsPath <- file.path( getwd(), "project", "batch", batchFldrContents )

  cat("Removing simulations from ./project/ \n", sep="")
  system(command=paste("rm -d -R ",simsPath,sep=""))

  cat("Removing batch files from ./project/batch folder\n", sep="")
  system(command=paste("rm -d -R ",batchContentsPath,sep=""))

  cat("Experiment and tidy-up complete, ready to start next experiment\n")
}


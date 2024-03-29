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
batchFiles    <- c("testColony.bch")
baseCtlFiles  <- c("simCtlFile.txt")
expPrefix     <- c("testColony")
plots         <- c(FALSE)

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

  # Now copy the project folder to Landmark NAS
  # copyDest <- file.path("/Volumes/home/thesisStuff/cwMSexperiments/TMB",paste(expPrefix[i],Sys.Date(),sep = "_") )
  # dir.create( copyDest )

  # # Copy project folder contents recursively to copyDest
  # cat( "Copying project folder contents to ", copyDest, "\n", sep = "" )
  # x <- file.copy( from = file.path( getwd(),"project"), to = copyDest, 
  #                 recursive = TRUE )

  # if(!x)
  # {
  #   cat(  "Error copying project folder contents to remote server,\n", 
  #         "using local dir instead" )
  #   copyDest <- file.path("../",paste(expPrefix[i],Sys.Date(),sep = "_") )
  #   dir.create( copyDest )

  #   # Copy project folder contents recursively to copyDest
  #   cat( "Copying project folder contents to ", copyDest, "\n", sep = "" )
  #   x <- file.copy( from = file.path( getwd(),"project"), to = copyDest, 
  #                   recursive = TRUE )
  # }

  # # Now remove the simulation folders and the contents of the batch sub-dir
  # # from the project folder
  # # sims <- grep(pattern = "sim", x = list.files("./project/"), value = T)
  # simsPath <- file.path( getwd(),"project",sims)
  # batchFldrContents <- list.files( file.path(getwd(), "project", "Batch") )
  # batchContentsPath <- file.path( getwd(), "project", "Batch", batchFldrContents )

  # # Copy out sims to dropbox, tidy up
  # cat("Removing simulations from ./project/ \n", sep="")
  # for(k in 1:length(simsPath))
  #   system(command=paste("rm -d -R ",simsPath[k],sep=""))

  # cat("Removing batch files from ./project/batch folder\n", sep="")
  # for(k in 1:length(batchContentsPath))
  #   system(command=paste("rm -d -R ",batchContentsPath[k],sep=""))

  # cat("Experiment and tidy-up complete, ready to start next experiment.\n")
}


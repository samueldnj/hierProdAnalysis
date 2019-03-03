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

  .runBatchJob( par = TRUE, prefix = expPrefix[i]  )
  
}
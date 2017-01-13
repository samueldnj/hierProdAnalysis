# Script to copy finished parallel batch project folders from
# unfinished parallel batch runs (error in a node)

# First, create the batch folder numbers
batchNumbers <- c()
remove <- FALSE

batchFolderNames <- paste("parBat4Corr",batchNumbers,sep="")
simFolderNames <- paste("sim_",batchFolderNames,sep="")
batchFolderNames <- file.path(getwd(),batchFolderNames)

for (i in 1:length(batchFolderNames))
{
  # Find the sim output folder in the project file
  batchDir <- file.path(batchFolderNames[i],"project",simFolderNames[i])
  # Set up the destination
  destination <- file.path(getwd(),"project")
  cat("\n", "Moving simulation ",i," sim folder to: ","\n",
    paste(getwd(),"/project/",sep=""),"\n", sep="")

  # Now copy the completed simulation
  file.copy(from=batchDir,to=destination,recursive=TRUE)
  if(remove)
  {
    cat("Removing folder ", batchFolderNames[i], "\n", sep="")
    system(command=paste("rm -d -R ",batchFolderNames[i],sep=""))
  }
  options(warn=1)  
}
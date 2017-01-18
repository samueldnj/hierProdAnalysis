# Little script to create a batch file for combinatorial combinations
# of scenarios or MPs in mseR style batch files

# Set up batch experiments by creating lists
# of parameters and different values at which to test them.
# Function below also requires abbreviations of the 

parList <-  list( "opMod$tUpeak" = seq(5,10,by=5),
                  "opMod$tUtrough " = seq(15,20,by=5),
                  "opMod$Umax" = seq(0.8,2,by=0.4) )
labels <- c("tUpk","tUtr","Umax")
prefix <- "5S-"
firstNo <- 1
# Function to create batch files from the lists above
batCreate <- function ( parList, outFile, label )
{
  exp <- expand.grid ( parList )
  for ( j in 1:nrow ( exp ))
  {
    # Create a label
    lab <- character()
    for ( k in 1:ncol ( exp ) ) 
      lab <- paste ( lab, label[k], as.numeric(exp[j,k]), sep = "" )
    lab <- paste(prefix,lab,sep="")
    # Print MP name
      if (j == 1) cat ( "# Scenario ", firstNo + j - 1, " : ", lab,
            "\n", sep = "", file = outFile, append = FALSE)
      else cat (  "# Scenario ", firstNo + j - 1, " : ", lab, 
                  "\n", sep = "", file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
      cat ( "scenario$scenario", firstNo + j - 1, "$ctrl$scenarioName '", lab, "'\n", sep = "", append = TRUE, 
            file = outFile )
      for ( k in 1:ncol ( exp ) )
        cat ( "scenario$scenario", firstNo + j - 1, "$", names(exp)[k], " ", exp[j,k], "\n", sep ="",
              file = outFile, append = TRUE )
      cat( "scenario$scenario", firstNo+j-1, "$ctrl$speciesNames c('Dover','English','Rock','Petrale')\n", sep ="",append=TRUE, file=outFile)
      cat( "scenario$scenario", firstNo+j-1, "$opMod$nS 4\n", sep="", append=TRUE, file=outFile )
      cat( "scenario$scenario", firstNo+j-1, "$opMod$q c(0.9,1.2,0.95,1.08)\n", sep="", append=TRUE, file=outFile )
      cat( "scenario$scenario", firstNo+j-1, "$opMod$tau2 c(0.01,0.01,0.01,0.01)\n", sep="", append=TRUE, file=outFile )
      cat( "scenario$scenario", firstNo+j-1, "$opMod$pars lisread('pars4Spec.txt')\n", sep="", append=TRUE, file=outFile )
      cat( "scenario$scenario", firstNo+j-1, "$opMod$parFile 'pars4Spec.txt'\n", sep="", append=TRUE, file=outFile )
      cat( "scenario$scenario", firstNo+j-1, "$assess$mBmsy c(40,80,30,60)\n", sep="", append=TRUE, file=outFile )
      cat( "scenario$scenario", firstNo+j-1, "$assess$sBmsy c(60,120,45,90)\n", sep="", append=TRUE, file=outFile )
      cat( "scenario$scenario", firstNo+j-1, "$assess$lnq_s c(0,0,0,0)\n", sep="", append=TRUE, file=outFile )

      # cat( "scenario$scenario", firstNo+j-1, "$ctrl$speciesNames c('Dover','English','Rock','Petrale','Arrowtooth')\n", sep ="",append=TRUE, file=outFile)
      # cat( "scenario$scenario", firstNo+j-1, "$opMod$nS 5\n", sep="", append=TRUE, file=outFile )
      # cat( "scenario$scenario", firstNo+j-1, "$opMod$q c(0.9,1.2,0.95,1.08,0.8)\n", sep="", append=TRUE, file=outFile )
      # cat( "scenario$scenario", firstNo+j-1, "$opMod$tau2 c(0.01,0.01,0.01,0.01,0.01)\n", sep="", append=TRUE, file=outFile )
      # cat( "scenario$scenario", firstNo+j-1, "$opMod$pars lisread('pars5Spec.txt')\n", sep="", append=TRUE, file=outFile )
      # cat( "scenario$scenario", firstNo+j-1, "$opMod$parFile 'pars5Spec.txt'\n", sep="", append=TRUE, file=outFile )
      # cat( "scenario$scenario", firstNo+j-1, "$assess$mBmsy c(40,80,30,60,100)\n", sep="", append=TRUE, file=outFile )
      # cat( "scenario$scenario", firstNo+j-1, "$assess$sBmsy c(60,120,45,90,150)\n", sep="", append=TRUE, file=outFile )
      # cat( "scenario$scenario", firstNo+j-1, "$assess$lnq_s c(0,0,0,0,0)\n", sep="", append=TRUE, file=outFile )

      cat ( "#\n", file = outFile, append = TRUE )
  }
  cat ( "# File Ends <not run>.", file = outFile, append = TRUE )
}


batCreate(parList,"autoBat.txt",labels)
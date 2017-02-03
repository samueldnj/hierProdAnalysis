# Little script to create a batch file for combinatorial combinations
# of scenarios or MPs in mseR style batch files

# Set up experiments to run batches over by creating lists
# of parameters and different values at which to test them.
# Function below also requires abbreviations of the 

scenParList <- list ( "opMod$nS" = c(2,4,5),
                      "opMod$corrOffDiag" = seq(0,0.8,by=0.2),
                      "opMod$kappaMult" = seq(2,2.6,by=0.2),
                      "opMod$lastNegCorr" = c(T,F)
                )

scenLabels <- c( "nS", "c", "m", "-c" )


# Function to create batch files from the lists above
scenarioCreate <- function ( parList, outFile, label )
{
  exp <- expand.grid ( parList )
  for ( j in 1:nrow ( exp ))
  {
    # Create a label
    lab <- character()
    for ( k in 1:ncol ( exp ) ) 
      lab <- paste ( lab, label[k], as.character(exp[j,k]), "_", sep = "" )
    # Print MP name
      if (j == 1) cat ( "# Scenario ", j, " : ", lab,
            "\n", sep = "", file = outFile, append = FALSE)
      else cat (  "# Scenario ", j, " : ", lab, 
                  "\n", sep = "", file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
      cat ( "scenario$scenario", j, "$ctrl$scenarioName '", lab, "'\n", sep = "", append = TRUE, 
            file = outFile )
      for ( k in 1:ncol ( exp ) )
        cat ( "scenario$scenario", j, "$", names(exp)[k], " ", exp[j,k], "\n", sep ="",
              file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
  }
  cat ( "# File Ends <not run>.", file = outFile, append = TRUE )
}

mpParList <- list ( "assess$msCorr" = c(T,F),
                    "assess$SigmaPriorCode" = c(0,1),
                    "assess$wishType" = c("diag","corr") )

mpLabels <- c("corr", "SigPrior", "wishType" )

mpCreate <- function ( parList, label, outFile )
{
  exp <- expand.grid ( parList )
  for ( j in 1:nrow ( exp ))
  {
    # Create a label
    lab <- character()
    for ( k in 1:ncol ( exp ) ) 
      lab <- paste ( lab, label[k], as.character(exp[j,k]),"_", sep = "" )
    # Print MP name
      if (j == 1) cat ( "# MP ", j, " : ", lab,
            "\n", sep = "", file = outFile, append = FALSE)
      else cat (  "# MP ", j, " : ", lab, 
                  "\n", sep = "", file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
      cat ( "mp$mp", j, "$ctrl$mpLabel '", lab, "'\n", sep = "", append = TRUE, 
            file = outFile )
      for ( k in 1:ncol ( exp ) )
        cat ( "mp$mp", j, "$", names(exp)[k], " ", exp[j,k], "\n", sep ="",
              file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
  }
  cat ( "# File Ends <not run>.", file = outFile, append = TRUE )
}


scenarioCreate(scenParList,"autoBat.txt",scenLabels)
# mpCreate ( mpParList,mpLabels,"mpAutoBat.txt")

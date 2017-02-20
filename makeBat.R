# Little script to create a batch file for combinatorial combinations
# of scenarios or MPs in mseR style batch files

# Set up experiments to run batches over by creating lists
# of parameters and different values at which to test them.
# Function below also requires abbreviations of the 

# - RE.bch (new ctl file structure)
scenParList <- list ( "opMod$nS" = c(2,4,5),
                      "opMod$corrOffDiag" = seq(0,0.8,by=0.2),
                      "opMod$kappaMult" = seq(0,2,by=0.5)
                )

scenLabels <- c( "nS", "c", "m" )


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
        cat ( "scenario$scenario", j, "$", names(exp)[k], " ", as.character(exp[j,k]), "\n", sep ="",
              file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
  }
  cat ( "# File Ends <not run>.", file = outFile, append = TRUE )
}

mpParList <- list ( "assess$s2lnUmsy" = c(0.0025,0.01,0.09,0.36),
                    "assess$sigU2P2" = c(0.01,0.1,0.2,0.3) )

mpLabels <- c("s2U", "betaU" )

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
      if (j == 1) cat ( "# MP ", j+16, " : ", lab,
            "\n", sep = "", file = outFile, append = FALSE)
      else cat (  "# MP ", j+16, " : ", lab, 
                  "\n", sep = "", file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
      cat ( "mp$mp", j+16, "$ctrl$mpLabel '", lab, "'\n", sep = "", append = TRUE, 
            file = outFile )
      for ( k in 1:ncol ( exp ) )
        cat ( "mp$mp", j+16, "$", names(exp)[k], " ", exp[j,k], "\n", sep ="",
              file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
  }
  cat ( "# File Ends <not run>.", file = outFile, append = TRUE )
}


# scenarioCreate(scenParList,"autoBat.txt",scenLabels)
mpCreate ( mpParList,mpLabels,"mpAutoBat.txt")

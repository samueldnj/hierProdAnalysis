# Little script to create a batch file for combinatorial combinations
# of scenarios or MPs in mseR style batch files

# Set up experiments to run batches over by creating lists
# of parameters and different values at which to test them.
# Function below also requires abbreviations of the 

parList <- list( "opMod$corrMult" = seq(0,1,by=0.2),
                 "opMod$kappaMult" = seq(0.2,2,by=0.2)
                           )
labels <- c("cM","vM")


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




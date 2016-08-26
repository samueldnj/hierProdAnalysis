# --------------------------------------------------------------------------
# control.R
# 
# Control script for the coastwide (single stock) simulation-estimation
# experiments comparing single species and multi-species assessment
# models.
# 
# Author: Samuel Johnson
# Date: 24 August, 2016
# 
# Sourcing scripts:
# 	simulation.R
# 	tools.R
#
# Estimator executables (ADMB):
#		ssProd
#		msProd
#
# Control file
#		controlFile.txt
#
# ToDo: 1. 	wrap the current code in functions to produce data for nS species
# 					and fit the single species model to each, returning posterior
# 					distributions - settle on a control list and data structure
# 					convention
#				2. 	simulate a covariance structure for epst and zetat (should we 
# 					just combine these?)
# 			3.  Code multispecies model (msProd) in ADMB - work off ssProd
# 					3a.	Estimate covariance in REs? Can we use ADMB-RE? Cholesky
# 							decomp is useful here...
# 			4.	Code TMB models (probably a separate script)
# 			5.	Ability to run multiple replicates with different seeds
# 			6. 	Tools to save single sets (control list settings) in project
# 					subfolder, blob style
# 			7. 	Plots!!!
# --------------------------------------------------------------------------

## Okay, let's start by sourcing the other scripts
source ( "simulation.R" )
source ( "tools.R" )

# Now read in control file
control <- lisread ( "controlFile.txt" )

# First, let's do a single species OM and get the AM to fit
# Set species number
s <- 1
# Create true biomass trends
bioList <- logProdModel ( msy = control$msy[s], Fmsy = control$Fmsy[s],
													nT = control$nT, Ft = rep (0.05, control$nT),
													sigma = control$sigma[s] )

# Now create observations and save to a list object
obsList <- list (nT = control$nT, Ct = bioList$Ct )
obsList $ It <- obsModel ( 	Bt = bioList$Bt, q = control$q[s], 
														tau = control$tau[s] )

# Now we want to turn obsList into an ADMB dat file (wrap this in a function)
obsList $ dumm <- 999
activeFileRoot <- "ssProd"
datFile <- paste (activeFileRoot, ".dat", sep = "")
cat ( "## Single Species Production Model data file, created ", 
			format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
			sep = "", file = datFile, append = FALSE )
lapply ( 	X = seq_along(obsList), FUN = writeADMB, x = obsList, 
					activeFile=datFile )

# Create a parameter list for writing a pin file
parList <- list ()
parList$lnMSY <- log( control$msy[s] )
parList$lnFMSY <- log( control$Fmsy[s] )
# parList$lnq <- log( control$q[s] )
parList$epst <- bioList$epst
# parList$lnsigma <- log ( control$sigma[s] )
# parList$lntau <- log ( control$tau[s] )
parList$mMSY <- control$msy[s]
parList$sMSY <- control$msy[s]
parList$mFMSY <- control$Fmsy[s]
parList$sFMSY <- control$Fmsy[s]
parList$alpha.sigma <- control$alpha.sigma[s]
parList$beta.sigma <- control$beta.sigma[s]
parList$alpha.tau <- control$alpha.tau[s]
parList$beta.tau <- control$beta.tau[s]
parList$mlnq <- log(control$q[s])
parList$slnq <- log(3) + log ( control$q[s])

# Write to pin file
pinFile <- paste ( activeFileRoot, ".pin", sep = "" )
cat ( "## Single Species Production Model pin file, created ", 
			format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
			sep = "", file = pinFile, append = FALSE )
lapply ( 	X = seq_along(parList), FUN = writeADMB, x = parList, 
					activeFile=pinFile )

# Now run the model
path <- getwd()
exec <- file.path(path,"ssProd")
mleCall <- paste ( exec, " -ainp ", pinFile, " -ind ", datFile, 
                      " -maxfn 5000", sep = "" )
system ( command = mleCall, wait =TRUE )

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
# 			3.  Estimate covariance in REs? Can we use ADMB-RE? Cholesky
# 							decomp is useful here...
# 			4.	Code TMB models (probably a separate script)
# 			5.	Ability to run multiple replicates with different seeds
# 			6. 	Tools to save single sets (control list settings) in project
# 					subfolder, blob style
# 			7. 	Plots!!!
# --------------------------------------------------------------------------

# Clean up
rm (list = ls())

## Okay, let's start by sourcing the other scripts
source ( "simulation.R" )
source ( "tools.R" )

# Now read in control file
ctlList <- lisread ( "controlFile.txt" )

# Create om
om <- opModel(seed = 1)

# Create data objects for EMs
datPar <- makeDataLists ( om, ctlList)

# Call EMs
ssRep <- list()
for (s in 1:ctlList$nS ) 
{
	cat ( "Fitting species ", s, "\n", sep = "")
	ssRep[[s]] <- callProcADMB ( 	dat = datPar$ssDat[[s]], par = datPar$ssPar[[s]],
																lab=s, fitTrials = 10, activeFileRoot = "ssProd",
																mcTrials = 20000, mcSave = 10, maxfn = 10000)
}
cat ( "Fitting ", ctlList$nS," species simultaneously.\n", sep = "")
msRep <- callProcADMB ( dat = datPar$msDat, par = datPar$msPar,
												fitTrials = 10, activeFileRoot = "msProd",
												mcTrials = 20000, mcSave = 10, maxfn = 10000 )



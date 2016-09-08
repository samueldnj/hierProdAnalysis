# --------------------------------------------------------------------------
# control.R
# 
# Control script for the coastwide (single stock) simulation-estimation
# experiments comparing single species and multi-species assessment
# models.
# 
# Author: Samuel D N Johnson
# Date: 24 August, 2016
# 
# Sourcing scripts:
# 	simulation.R
# 	tools.R
#   stats.R
#   plots.R
#
# Estimator executables (ADMB):
#		ssProd
#		msProd
#
# Control file
#		controlFile.txt
#
# ToDo: 2. 	simulate a covariance structure for epst and zetat (should we 
# 					just combine these?)
# 			3.  Estimate covariance in REs? Can we use ADMB-RE? Cholesky
# 							decomp is useful here...
# 			4.	Code TMB models (probably a separate script)
# 			6. 	Tools to save single sets (control list settings) in project
# 					subfolder, blob style
# 			7. 	Plots!!!
# --------------------------------------------------------------------------

# Clean up
rm (list = ls())

## Okay, let's start by sourcing the other scripts
source ( "simulation.R" )
source ( "tools.R" )
source ( "stats.R" )

# Now read in control file
ctlList <- lisread ( "controlFile.txt" )

# Test simEstProc
blob <- simEstProc ( ctl = ctlList, quiet = TRUE )

blob <- makeRelErrorDists ( blob )


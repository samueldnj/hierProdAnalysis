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
#   packages.R
# 	simulation.R
# 	tools.R
#   stats.R
#   plots.R
#
# Estimator executables (ADMB):
#		ssProdCV
#		msProdCV
#
# Control file
#		controlFile.txt
#
# ToDo: 4.	Code TMB models (probably a separate script)
#       5.  Missing values in msProd and ssProd
# 
# --------------------------------------------------------------------------

# Clean up
rm (list = ls())

## Okay, let's start by sourcing the other scripts
source ( "packages.R" )
source ( "simulation.R" )
source ( "tools.R" )
source ( "stats.R" )
source ( "plots.R" )


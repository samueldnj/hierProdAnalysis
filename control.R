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
#   init.R
# 	simulation.R
# 	tools.R
#   stats.R
#   plots.R
#
# Estimator executables (ADMB):
#		msProd
#
# Control file
#		simCtlFile.txt
#
# ToDo: - Shared prior on Umsy?
#       - 
# 
# --------------------------------------------------------------------------

# Clean up
rm (list = ls())
gc()

## Okay, let's start by sourcing the other scripts
source ( "init.R" )
source ( "simulation.R" )
source ( "tools.R" )
source ( "stats.R" )
source ( "plots.R" )


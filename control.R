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
# ToDo: 5.  Missing values in msProd and ssProd
# 
# --------------------------------------------------------------------------

## source other scripts
source ( "init.R" )
source ( "simulation.R" )
source ( "tools.R" )
source ( "stats.R" )
source ( "plots.R" )


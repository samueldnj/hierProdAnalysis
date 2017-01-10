# --------------------------------------------------------------------------
# packages.R
# 
# Script for loading packages for coastwide multispecies hierarchical
# assessment simulation estimation procedure
# 
# Author: Samuel D N Johnson
# Date: October 19, 2016
#
# --------------------------------------------------------------------------

# Load packages
library ( "coda" )
library ( "dplyr" )
library ( "TMB" )

# compile and load msProd objective function.
compile ("msProd.cpp")
dyn.load(dynlib("msProd"))
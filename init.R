# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# init.R
# 
# Checks if required packages are installed. If not, installs them.
# Then loads all required packages.
# 
# 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

cranPackages <- c("coda",
                  "dplyr"
                  "reshape2",
                  "TMB",
                  "raster",
                  "grid",
                  "RColorBrewer",
                  "HapEstXXR",
                  "parallel",
                  "stringr",
                  "wesanderson",
                  "scales",
                  "beepr" )


for( pkg in cranPackages )
  while(!require(pkg, character.only = TRUE) )
    install.packages( pkg )


githubPackages <- c(ggsidekick = "seananderson/ggsidekick")

for( pkgIdx in 1:length(githubPackages) )
  while(!require(names(githubPackages)[pkgIdx], character.only = TRUE))
    devtools::install_github(githubPackages[pkgIdx])


# compile and load msProd objective function.
compile ("msProd.cpp")
dyn.load(dynlib("msProd"))
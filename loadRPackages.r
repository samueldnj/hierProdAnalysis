packNames <- c( "animation", "coda", "bitops", "data.table","deSolve",
                "devtools", "evaluate", "geoR", "geoRglm", "ggplot2","httr", "lme4","maptools",
                "matrixcalc", "MCMCpack", "memoise","msm", "multicore", "mvtnorm",
                "OpenCL", "PBSmapping", "PBSmodelling", "RandomFields",  "randomForest",
                "RCurl", "rehsape2", "dplyr","TMB","raster","grid","RColorBrewer",
                "ref", "rlecuyer", "RMark", "RODBC", "scatterplot3d", "snow",
                "snowfall", "sp", "splancs", "stringr", "whisker", "XML",
                "HapEstXXR"
              )

for( i in 1:length(packNames) )
  install.packages( packNames[i] )

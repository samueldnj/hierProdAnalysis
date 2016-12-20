mcPar <- read.table("msProdParMC.dat",header=TRUE)
mcPar1 <- as.mcmc(mcPar[seq(nrow(mcPar)/2+1,nrow(mcPar),by=4),])
mcPar2 <- as.mcmc(mcPar[seq(nrow(mcPar)/2+2,nrow(mcPar),by=4),])

plot(mcPar1)

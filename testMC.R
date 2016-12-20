mcPar <- read.table("msProdParMC.dat",header=TRUE)
mcPar1 <- as.mcmc(mcPar[seq(nrow(mcPar)/2+1,nrow(mcPar),by=2),])
mcPar2 <- as.mcmc(mcPar[seq(nrow(mcPar)/2+2,nrow(mcPar),by=2),])

plot(mcPar1)

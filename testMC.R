mcPar <- read.table("msProdParMC.dat",header=TRUE)
mcPar1 <- as.mcmc(mcPar[seq(nrow(mcPar)/2+1,nrow(mcPar),by=3),])
mcPar2 <- as.mcmc(mcPar[seq(nrow(mcPar)/2+2,nrow(mcPar),by=3),])
mcPar3 <- as.mcmc(mcPar[seq(nrow(mcPar)/2+3,nrow(mcPar),by=3),])


plot(mcPar1)

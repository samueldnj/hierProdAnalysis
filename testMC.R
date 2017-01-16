mcPar <- read.table("msProdParMC.dat",header=TRUE)
nS <- 1
mcPar1 <- as.mcmc(mcPar[seq(nrow(mcPar)/2+1,nrow(mcPar),by=nS),])
# mcPar2 <- as.mcmc(mcPar[seq(nrow(mcPar)/2+2,nrow(mcPar),by=nS),])

plot(mcPar1)

# chains <- mcmc.list(mcPar1,mcPar2)

mcChain <- read.table("mcout.dat",header=TRUE)
mcChain <- mcChain[(nrow(mcChain)/2+1):nrow(mcChain),]
library(coda)
mcChain <- as.mcmc(mcChain)


samples <- nrow(mcChain)
chain1 <- as.mcmc(mcChain[1:(samples/2),])
chain2 <- as.mcmc(mcChain[(samples/2+1):samples,])

chains <- mcmc.list(chain1,chain2)




plot(mcChain)

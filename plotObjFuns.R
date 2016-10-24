# This is a little script designed to help plot objective functions around
# neighbourhoods of given points, hopefully allowing us to find the
# issue with convergence, and identify correlated variables

ssProdObjFun <- function ( par, dat )
{ 
  # Recover dat stuff
  nT          <- dat$nT
  Ct          <- dat$Ct
  It          <- dat$It

  # prod model pars
  Bmsy        <- exp(par$lnBmsy)
  Fmsy        <- exp(par$lnFmsy)
  tau2        <- exp(par$lnTau2)
  sigma2      <- exp(par$lnsigma2)
  q           <- exp(par$lnq)
  mBmsy       <- par$mBmsy
  sBmsy       <- par$sBmsy
  mFmsy       <- par$mFmsy
  sFmsy       <- par$sFmsy
  alpha_sigma <- par$alpha_sigma
  beta_sigma  <- par$beta_sigma
  alpha_tau   <- par$alpha_tau
  beta_tau    <- par$beta_tau
  mlnq        <- par$mlnq
  slnq        <- par$slnq
  rho         <- par$rho
  epst        <- par$epst

  # First, run state dynamics
  epstCorr <- numeric(length=nT)
  epstCorr[1] <- epst[1]

  Bt <- numeric(length=nT)

  Bt[1] <- 2*Bmsy * exp(epstCorr[1])

  for (t in 1:(nT-1))
  {
    epstCorr[t+1] <- rho*epstCorr[t] + epst[t+1]

    Bt[t+1] <- Bt[t] + 2*Fmsy*Bt[t] * (1 - Bt[t]/Bmsy/2) - Ct[t]
    Bt[t+1] <- Bt[t+1] * exp(epstCorr[t+1])
    Bt[t+1] <- max(1,Bt[t+1])
  }

  # Compute expected indices
  It_bar <- q * Bt

  # Take residuals
  SSRobs <- sum((log(It) - log(It_bar))^2)

  # Now calculate likelihoods
  obsNLL  <- nT*log(tau2)/2 + SSRobs/2/tau2

  procNLL <- nT*log(sigma2)/2 + sum(epst^2)/2/sigma2

  BmsyPrior   <- (Bmsy - mBmsy)^2/2/sBmsy/sBmsy
  FmsyPrior   <- (Fmsy - mFmsy)^2/2/sFmsy/sBmsy
  sigma2Prior <- (alpha_sigma + 1 ) * log(sigma2) + beta_sigma/sigma2
  tau2Prior   <- (alpha_tau + 1 ) * log(tau2) + beta_tau/tau2
  lnqPrior    <- (log(q) - mlnq)^2/2/slnq/slnq

  objFun <- obsNLL + procNLL + BmsyPrior + FmsyPrior + sigma2Prior + tau2Prior + lnqPrior

 return(objFun) 
}

plot2DLikeSS <- function (par,dat,
                          par1="lnBmsy",par2="lnFmsy",
                          lim1=c(5,6),lim2=c(-3,-2),
                          nGrid=100 )
{
  # Centre at given values
  par1idx <- which (grepl(par1,names(par)))[1]
  par2idx <- which (grepl(par2,names(par)))[1]
  par1Centre <- par[[par1idx]]
  par2Centre <- par[[par2idx]]


  # Create grids
  par1Seq <- seq ( from = lim1[1], to=lim1[2], length = nGrid )
  par2Seq <- seq ( from = lim2[1], to=lim2[2], length = nGrid )

  likeGrid <- matrix(NA, nrow = length(par1Seq), ncol=length(par2Seq) )

  for (i in 1:length(par1Seq))
  {
    for(j in 1:length(par2Seq))
    {
      par[[par1idx]] <- par1Seq[i]
      par[[par2idx]] <- par2Seq[j]
      likeGrid[i,j] <- ssProdObjFun(par,dat)
    }
  }

  image(x = par2Seq, y = par1Seq, z = likeGrid, col = terrain.colors(12),
        xlab=par2,ylab=par1)
}

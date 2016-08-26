# --------------------------------------------------------------------------
# simulation.R
# 
# Script for the simulation procedure to produce data for ADMB and TMB
# estimators.
# 
# Author: Samuel Johnson
# Date: 24 August, 2016
# 
# --------------------------------------------------------------------------


# logProdModel()
# Function to simulate a logistic based surplus production model,
# assuming that population is initially at equilibrium (2*Bmsy).
# inputs:     msy=maximum sustainable yield; Fmsy=fishing mortality for MSY
#             nT=length of simulation; Ct=nT-vector of catch (opt)
# 						Ft=nT-vector of fishing mortality (opt)
#             epst=nT-vector of proc errors (could be RW, or IID)
#             sigma = process error sd (uncorr)
# outputs:    Bt=nT-vector of modeled biomass; Ct=nT-vector of catch
# 						Ut=vector of exploitation rates; 
# usage:      Bexp(eps) = indices of abundance for estimation w/ logN errors eps
logProdModel <- function ( msy = 1, Fmsy = 0.1, nT = 50, Ft = NULL, Ct = NULL,
                           epst = rnorm ( n = nT ), sigma = 0.2 )
{
  # First, initialise a vector to hold biomass
  Bt <- numeric ( length = nT )

  # Multiply errors by variance
  epst <- epst * sigma

  # Compute bias correction
  sigma2 <- sigma * sigma

  if ( is.null ( Ft ) & is.null ( Ct ) ) 
  {
    cat ( "You're missing harvest, stupid.\n" )
    return ()
  }

  if ( is.null ( Ct ) ) 
  {
  	Ut <- 1 - exp( -Ft )
  	Ct <- numeric ( length = nT )
  }

  # # Compute U_MSY
  # Umsy <- (1 - exp(-Fmsy) )

  # Compute BMSY
  Bmsy <- msy / Fmsy

  # Populate the vector, with a special case for t = 1
  Bt [ 1 ] <- 2 * Bmsy * exp ( epst[1] - sigma2 / 2. )

  # Loop over remaining years
  for ( t in 2:nT )
  {
    if ( Ct [ t-1 ] == 0 ) Ct [ t-1 ] <- Ut[t-1] * Bt [ t-1 ]
    Bt [ t ] <-   (	Bt [ t-1 ] +                       # prev year's B
                  	2. * Fmsy * Bt [ t-1 ] * 
                  	( 1. - Bt [ t-1 ] / Bmsy / 2. ) -  # recruitment
                  	Ct [ t-1 ]                         # Catch
                  ) * exp ( epst [ t ] - sigma2 / 2. ) # Process error
}  # End loop for biomass

  # compute depletion value for comparison to model fit
  dep <- Bt [ nT ] / Bmsy / 2.

  # Return B, C, F and U series and
  # depletion value
  outList <- list ()
  outList $ Bt <- Bt
  outList $ Ct <- Ct
  outList $ Ft <- Ft
  outList $ Ut <- Ut
  outList $ epst <- epst
  outList $ dep <- dep
  return ( outList )
}

# Why keep the observation model separate???? 
# Right now just for modularity.

# obsModel()
# A function which creates observations (indices of biomass) from
# given true series of biomass and observation model parameters.
# inputs:			Bt=true series of biomass; q=survey catchability parameter
#							nT=length of time series; deltat=nT-vector of standard errors
# 						tau=obs error standard deviation (uncorrelated)
# outputs:		It=nT-vector of survey observations
# usage:			creates survey observations which are supplied to the estimator
obsModel <- function ( 	Bt, q = 0.05, nT = length (Bt), 
												deltat = rnorm(nT), tau = 0.4 )
{
	It <- q * Bt * exp ( deltat - tau*tau/2.0 )
	return ( It )
}



// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// msProd_TMB.cpp
// 
// A multi-stock surplus production (Schaefer) state-space model with 
// joint prior distributions applied to catchability (q) 
// and productivity (U_msy)
// 
// Author: Samuel Johnson
// Date: 1 November, 2016
// Purpose: The assessment model in a simulation study of multi-stock
// robin hood methods.
// 
// Updates:
//    May 13, 2017: Added multiple surveys, ready for application
//                  to real DERPA data
// 
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

#include <TMB.hpp>       // Links in the TMB libraries
#include <iostream>

// posfun
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

// invLogit
template<class Type>
Type invLogit(Type x, Type scale, Type trans){
  return scale/(Type(1.0) + exp(-Type(1.0)*x)) - trans;
}

// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /*data section*/
  // Data Structures
  DATA_ARRAY(It);               // CPUE data
  DATA_ARRAY(Ct);               // Catch data
  // Model dimensions
  int nO = It.dim(0);           // No. of surveys
  int nS = It.dim(1);           // No. of species
  int nT = It.dim(2);           // No of time steps

  // Model switches
  DATA_INTEGER(SigmaPriorCode); // 0 => IG on diagonal element, 1 => IW on cov matrix
  DATA_INTEGER(tauqPriorCode);  // 0 => IG on tauq2, 1 => normal
  DATA_INTEGER(sigUPriorCode);  // 0 => IG on sigU2, 1 => normal
  DATA_INTEGER(lnqPriorCode);   // 0 => hyperprior, 1 => multilevel
  DATA_INTEGER(lnUPriorCode);   // 0 => hyperprior, 1 => multilevel 
  DATA_IVECTOR(initT);          // first year of reconstructed catch hist
  DATA_IVECTOR(initBioCode);    // initial biomass at 0 => unfished, 1=> fished


  /*parameter section*/
  // Leading Parameters
  PARAMETER_VECTOR(lnBmsy);             // Biomass at MSY
  PARAMETER_VECTOR(lnUmsy);             // Optimal exploitation rate            
  PARAMETER_VECTOR(lntau2_o);           // survey obs error var
  PARAMETER_ARRAY(lnq_os);              // Survey-Species spec catchability
  PARAMETER_VECTOR(lnBinit);            // Non-equilibrium initial biomass
  // Priors
  PARAMETER_VECTOR(lnqbar_o);           // survey mean catchability
  PARAMETER_VECTOR(lntauq_o);           // survey catchability sd 
  PARAMETER(mq);                        // hyperprior mean catchability
  PARAMETER(sq);                        // hyperprior catchability sd
  PARAMETER(lnUmsybar);                 // shared prior mean Umsy
  PARAMETER(lnsigUmsy);                 // shared prior Umsy sd
  PARAMETER(mUmsy);                     // hyperprior mean Umsy
  PARAMETER(sUmsy);                     // hyperprior Umsy sd
  PARAMETER_VECTOR(mBmsy);              // prior mean eqbm biomass
  PARAMETER_VECTOR(sBmsy);              // prior eqbm biomass var
  PARAMETER_VECTOR(tau2IGa);            // IG parameters for tau2 prior
  PARAMETER_VECTOR(tau2IGb);            // IG parameters for tau2 prior
  PARAMETER_VECTOR(tauq2Prior);         // Hyperparameters for tauq2 prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(sigU2Prior);         // Hyperparameters for sigU2 prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(kappa2IG);           // IG parameters for kappa2 prior
  PARAMETER_VECTOR(Sigma2IG);           // IG parameters for Sigma2 prior
  PARAMETER_MATRIX(wishScale);          // IW scale matrix for Sigma prior
  PARAMETER(nu);                        // IW degrees of freedom for Sigma prior    
  // Random Effects
  PARAMETER_VECTOR(eps_t);              // year effect  
  PARAMETER(lnkappa2);                  // year effect uncorrelated variance
  PARAMETER_ARRAY(zeta_st);             // species effect
  PARAMETER(lnSigmaDiag);               // Species effect cov matrix diag
  PARAMETER_VECTOR(SigmaDiagMult);      // Sigma diagonal mults
  PARAMETER_VECTOR(logitSigmaOffDiag);  // Species effect corr chol factor off diag  
  PARAMETER(logit_gammaYr);             // AR1 auto-corr on year effect (eps)
  

  // State variables
  array<Type>       Bt(nS,nT);
  // Leading parameters
  vector<Type>      Bmsy    = exp(lnBmsy);
  vector<Type>      Umsy    = exp(lnUmsy);
  vector<Type>      tau2_o  = exp(lntau2_o);
  vector<Type>      Binit   = exp(lnBinit);
  array<Type>       q_os(nO,nS);
  for( int o = 0; o < nO; o++ )
  {
    for( int s = 0; s < nS; s++ )
    {
      q_os(o,s) = exp(lnq_os(o,s));
    }
  }
  
  // Prior hyperpars
  vector<Type>      qbar_o  = exp(lnqbar_o);
  vector<Type>      tauq_o  = exp(lntauq_o);
  vector<Type>      tauq2_o = exp(Type(2) * lntauq_o);
  Type              Umsybar = exp(lnUmsybar);
  Type              sigUmsy = exp(lnsigUmsy);
  Type              sigUmsy2= exp(Type(2) * lnsigUmsy);
  // Random Effects
  Type              kappa2  = exp(lnkappa2);
  Type              gammaYr = Type(1.96) / (Type(1.) + exp(Type(-2.)*logit_gammaYr) ) - Type(0.98);
  vector<Type>      SigmaOffDiag(nS*(nS-1)/2);
                    SigmaOffDiag  = Type(1.96) / (Type(1.) + exp(Type(-2.)*logitSigmaOffDiag) ) - Type(0.98);     
  vector<Type>      SigmaDiag = exp(lnSigmaDiag) * SigmaDiagMult;
  matrix<Type>      SigmaCorr(nS,nS);
  matrix<Type>      SigmaD(nS,nS);
  matrix<Type>      Sigma(nS,nS);
  vector<Type>      omegat(nT);
  // Scalars
  Type              nll     = 0.0;   // objective function (neg log likelihood)
  Type              nllRE   = 0.0;   // process error likelihood
  Type              pospen  = 0.0;   // posfun penalty
  // Derived variables
  vector<Type>      msy = Bmsy * Umsy;
  array<Type>       Ut(nS,nT);
  vector<Type>      DnT;

  // Procedure Section //
  // First create the correlation matrix for the species effects
  UNSTRUCTURED_CORR_t<Type> specEffCorr(SigmaOffDiag);
  SigmaCorr = specEffCorr.cov();
  // Standard deviation component
  SigmaD.fill(0.0);
  SigmaD.diagonal() = sqrt(SigmaDiag);
  // Now combine
  Sigma = (SigmaD * SigmaCorr);
  Sigma *= SigmaD;

  // Generate Population Dynamics
  // initialise first year effect at 0, then loop to fill in
  omegat(0) = 0;
  for( int t=1; t<nT; t++ )
  {
    omegat(t) = gammaYr * omegat(t-1) + (1 - gammaYr) * eps_t(t-1);  
  }
  
  Bt.fill(-1);
  // Now loop over species, reconstruct history from initT
  for( int s = 0; s < nS; s++ )
  {
    // initialise population, if initBioCode(s)==0 this will be eqbm
    if( initBioCode(s) == 0 ) Bt(s,initT(s)) = Type(2) * Bmsy(s);
    if( initBioCode(s) == 1 ) Bt(s,initT(s)) = Binit(s);
    for( int t = initT(s)+1; t < nT; t++ )
    {
      Bt(s,t) = Bt(s,t-1) + Bt(s,t-1)*Umsy(s) * (Type(2.0) - Bt(s,t-1)/Bmsy(s)) - Ct(s,t-1);
      Bt(s,t) = posfun(Bt(s,t),Type(2e-3),pospen);
      Bt(s,t) *= exp(omegat(t) + zeta_st(s,t-1));
      nllRE += Type(1000.0)*pospen;      
    }
  }
  // Loop again to add species and year effects
  for (int t=1; t<nT; t++)
  {
    // Add year effect contribution to objective function
    nllRE += Type(0.5)*(lnkappa2 + pow( eps_t( t - 1 ), 2 ) / kappa2 ) + eps_t(t-1);
    // Add correlated species effects contribution to likelihood
    if (nS > 1) nllRE += VECSCALE(specEffCorr,sqrt(SigmaDiag))(zeta_st.col(t-1));
  }
  // add REs to joint nll
  nll += nllRE;

  // Concentrate species specific obs error likelihood?
  // Sum of squares vector
  array<Type>   ss_os(nO,nS);
  vector<Type>  tau2hat_o(nO);
  array<Type>   validObs(nO,nS);
  array<Type>   qhat(nO,nS);
  array<Type>   zSum(nO,nS);
  // Fill with 0s
  validObs.fill(0.0);
  zSum.fill(0.0);
  ss_os.fill(0.0);
  Type nllObs = 0.0;
  // qhat.fill(0.0);
  // tau2hat_o.fill(0.0);
  // Compute observation likelihood
  // Loop over surveys
  for( int o=0; o < nO; o++ )
  {
    // species
    for( int s = 0; s < nS; s++ )
    {
      // times
      for( int t = initT(s); t < nT; t++ )
      {
        // only add a contribution if the data exists (Ist < 0 is missing)
        if (It(o,s,t) > 0) 
        {
          validObs(o,s) += int(1);
          Type res = log( It( o, s, t ) ) - log( Bt( s, t ) );
          zSum(o,s) += res;
          ss_os(o,s) += pow( res - lnq_os(o,s), 2 );
          nllObs += Type(0.5) * ( lntau2_o(o) + pow( res - lnq_os(o,s), 2 ) / tau2_o(o) );
        }       
      }
    }
    tau2hat_o( o ) = ss_os.matrix().row(o).sum() / ( validObs.matrix().row(o).sum() - 1 );
    for( int s = 1; s < nS; s++ )
    {
      if( nS == 1) qhat(o,s) = exp( zSum(o,s) / validObs(o,s) );
      if( nS > 1 ) qhat(o,s) = exp( ( zSum(o,s) / tau2_o(o) + lnqbar_o(o)/tauq_o(o) ) / ( validObs(s) / tau2hat_o(o) + 1 / tauq2_o(o) ) );  
    }   
  }
  nll += nllObs;

  // Add priors
  Type nllBprior = 0.0;
  Type nllqPrior = 0.0;
  Type nllUprior = 0.0;

  // eqbm biomass
  for (int s=0; s<nS; s++ )
  {
    nllBprior += Type(0.5) * pow( Bmsy(s) - mBmsy(s), Type(2) ) / sBmsy(s) / sBmsy(s); 
    if(initBioCode(s) == 1) nllBprior +=  pow( Binit(s) - mBmsy(s), 2 ) / pow(sBmsy(s)/Type(2.0),Type(2.0)) ; 
  }

  // multispecies shared priors
  if (nS > 1)
  {
    // 1st level priors
    for( int s = 0; s < nS; s++ )
    {
      // In here, loop over surveys
      for( int o = 0; o < nO; o++ )
      {
        // catchability
        // Shared Prior
        if( lnqPriorCode == 1 )
        {
          nllqPrior +=  lntauq_o(o) + Type(0.5) * pow( lnq_os(o,s) - lnqbar_o(o), 2 ) / tauq2_o(o) ;  
        }
        // No shared prior (uses the same prior as the SS model)
        if( lnqPriorCode == 0 )
        {
          nllqPrior += Type(0.5) * pow( q_os(o,s) - mq, 2 ) / sq / sq;  
        }
      }
      // productivity
      // Shared Prior
      if( lnUPriorCode == 1 )
      {
        nllUprior += lnsigUmsy + Type(0.5) * pow(Umsy(s) - Umsybar, 2 ) / sigUmsy2 ;
      }
      // No shared prior (uses SS model prior)
      if( lnUPriorCode == 0 )
      {
        nllUprior += Type(0.5) * pow( Umsy(s) - mUmsy, 2 ) / sUmsy / sUmsy;
      }
    }  
    // Hyperpriors
    // catchability
    if( lnqPriorCode == 1 )
    {
      for( int o = 0; o < nO; o++ ) nllqPrior += Type(0.5) * pow( qbar_o(o) - mq, 2 ) / sq / sq;  
    }
    // productivity
    if( lnUPriorCode == 1 )
    {
      nllUprior += Type(0.5) * pow( Umsybar - mUmsy, 2 ) / sUmsy / sUmsy;
    }
    // End multispecies shared priors

    // If not using shared priors, use the hyperpriors for all
    // species specific
  } 
  if( nS == 1 ) 
  {
    // Now for single species model
    for (int o=0; o<nO; o++)
    {
      // catchability
      nllqPrior += Type(0.5) * pow( exp(lnq_os(o,0)) - mq, 2 ) / sq / sq;
    }
    // productivity
    nllUprior += Type(0.5) * pow( Umsy(0) - mUmsy, 2 ) / sUmsy / sUmsy;
  }
  nll += nllBprior +  nllqPrior + nllUprior;
  
  // Variance IG priors
  // Obs error var
  Type nllVarPrior = 0.0;
  for( int o = 0; o < nO; o++ )
  {
    nllVarPrior += (tau2IGa(o)+Type(1))*lntau2_o(o)+tau2IGb(o)/tau2_o(o);  
  }
  // year effect deviations var
  nllVarPrior += (kappa2IG(0)+Type(1))*lnkappa2 + kappa2IG(1)/kappa2;
  // Now multispecies priors
  if (nS > 1)
  {
    // shared q prior variance
    if( tauqPriorCode == 0 )
    {
      for( int o = 0; o < nO; o++)
      {
        nllVarPrior += (tauq2Prior(0)+Type(1))*Type(2)*lntauq_o(o)+tauq2Prior(1)/tauq2_o(o);    
      }
    }
    if( tauqPriorCode == 1 )
    {
      for( int o = 0; o < nO; o++)
      {
        nllVarPrior += Type(0.5) * pow( tauq2_o(o) - tauq2Prior(0), 2) / tauq2Prior(1);
      }
    }
    // shared U prior variance
    if( sigUPriorCode == 0 )
    {
      nllVarPrior += (sigU2Prior(0)+Type(1))*Type(2.)*lnsigUmsy+sigU2Prior(1)/sigUmsy2;  
    }
    if( sigUPriorCode == 1 )
    {
      nllVarPrior += Type(0.5) * pow( sigUmsy2 - sigU2Prior(0), 2) / sigU2Prior(1);
    }
    
    // Apply Sigma Prior
    if( SigmaPriorCode == 0 ) // Apply IG to estimated SigmaDiag element
    {
      nllVarPrior += (Sigma2IG(0)+Type(1))*lnSigmaDiag+Sigma2IG(1)/exp(lnSigmaDiag);
    }
    if( SigmaPriorCode == 1 ) // Apply IW prior to Sigma matrix
    {
      matrix<Type> traceMat = wishScale * Sigma.inverse();
      Type trace = 0.0;
      for (int s=0;s<nS;s++) trace += traceMat(s,s);
      nllVarPrior += Type(0.5) *( (nu + nS + 1) * atomic::logdet(Sigma) + trace);
    }
  }
  nll += nllVarPrior;

  // Derive some output variables
  Ut  = Ct / Bt;
  DnT = Bt.col(nT-1)/Bmsy/2;

  // Reporting Section //
  // Variables we want SEs for
  ADREPORT(Bt);
  ADREPORT(q_os);
  ADREPORT(Bmsy);
  ADREPORT(Umsy);
  ADREPORT(msy);
  ADREPORT(Binit);
  ADREPORT(tau2_o);
  ADREPORT(kappa2);
  if (nS > 1 )
  {
    ADREPORT(SigmaDiag);
    ADREPORT(Umsybar);
    ADREPORT(sigUmsy2);
    ADREPORT(qbar_o);
    ADREPORT(tauq2_o)
  }
  ADREPORT(gammaYr);
  
  // Everything else //
  REPORT(Bt);
  REPORT(Ut);
  REPORT(Binit);
  REPORT(DnT);
  REPORT(q_os);
  REPORT(msy);
  REPORT(Bmsy);
  REPORT(Umsy);
  REPORT(msy);
  REPORT(tau2_o);
  REPORT(kappa2);
  REPORT(nT);
  REPORT(nO);
  REPORT(nS);
  REPORT(tau2hat_o);
  REPORT(qhat);
  if (nS > 1)
  {
    REPORT(Sigma);
    REPORT(SigmaCorr);
    REPORT(SigmaDiag);
    REPORT(Umsybar);
    REPORT(sigUmsy2);
    REPORT(qbar_o);
    REPORT(tauq2_o);
  }
  REPORT(gammaYr);
  REPORT(nll);
  REPORT(nllRE);
  REPORT(nllObs);
  REPORT(nllqPrior);
  REPORT(nllUprior);
  REPORT(nllBprior);
  REPORT(nllVarPrior);
  
  return nll;
}

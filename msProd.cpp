// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// msProd_TMB.cpp
// 
// TMB version of of the multispecies production model
// 
// Author: Samuel Johnson
// Date: August 25, 2014                    
// Purpose: To fit a multi-species RH production model to commercial catch
// and survey CPUE data.
// 
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

#include <TMB.hpp>                                // Links in the TMB libraries
#include <iostream>
// #include <convenience.hpp>

// posfun
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}


// invLogit
// template<class Type>
// Type invLogit(Type x, Type scale, Type trans){
//   return scale/(Type(1.0) + exp(-Type(1.0)*x)) - trans;
// }

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
  // indexing variables
  int nS=It.dim(0);
  int nT=It.dim(1);

  /*parameter section*/
  // Leading Parameters
  PARAMETER_VECTOR(lnBmsy);             // Biomass at MSY
  PARAMETER_VECTOR(lnUmsy);             // Optimal exploitation rate
  PARAMETER(lntau2);                    // Obs error variance
  PARAMETER_VECTOR(tau2mult);           // species mults for observation error
  PARAMETER_VECTOR(lnq);                // Species specific catchability
  // Priors
  PARAMETER(lnqbar);                    // prior mean catchability
  PARAMETER(lntauq2);                   // prior catchability variance 
  PARAMETER(mlnq);                      // hyperprior mean catchability
  PARAMETER(s2lnq);                     // hyperprior catchability variance
  PARAMETER(lnUmsybar);                 // shared prior mean Umsy
  PARAMETER(lnsigUmsy2);                // shared prior Umsy variance
  PARAMETER(mlnUmsy);                   // hyperprior mean Umsy
  PARAMETER(s2lnUmsy);                  // hyperprior Umsy variance
  PARAMETER_VECTOR(mlnBmsy);            // prior mean eqbm biomass
  PARAMETER_VECTOR(s2lnBmsy);           // prior eqbm biomass var
  PARAMETER_VECTOR(tau2IG);             // IG parameters for tau2 prior
  PARAMETER_VECTOR(tauq2IG);            // IG parameters for tauq2 prior
  PARAMETER_VECTOR(sigU2IG);            // IG parameters for tauq2 prior
  PARAMETER_VECTOR(kappa2IG);           // IG parameters for kappa2 prior
  PARAMETER_VECTOR(Sigma2IG);           // IG parameters for Sigma2 prior

  // Random Effects
  PARAMETER_VECTOR(eps_t);              // year effect
  PARAMETER(lnkappa2);                  // year effect uncorrelated variance
  PARAMETER(logit_gammaYr);             // AR1 auto-corr on year effect
  PARAMETER_ARRAY(zeta_st);             // species effect
  PARAMETER(lnSigmaDiag);               // Species effect cov matrix diag
  PARAMETER_VECTOR(SigmaDiagMult);      // Sigma diagonal mults
  PARAMETER_VECTOR(logitSigmaOffDiag);  // Species effect corr chol factor off diag  
  
  
  // State variables
  array<Type> lnBt(nS,nT);
  // Leading parameters
  vector<Type> Bmsy = exp(lnBmsy);
  vector<Type> Umsy = exp(lnUmsy);
  Type Umsybar      = exp(lnUmsybar);
  Type tau2         = exp(lntau2);
  vector<Type> q    = exp(lnq);
  Type qbar         = exp(lnqbar);
  // Prior hyperpars
  Type tauq2        = exp(lntauq2);
  Type sigUmsy2     = exp(lnsigUmsy2);
  // Random Effects
  Type kappa2       = exp(lnkappa2);
  Type gammaYr      = Type(2.) / (Type(1.) + exp(Type(-5.)*logit_gammaYr) ) - Type(1.);
  vector<Type> SigmaOffDiag(nS*(nS-1)/2);  
  SigmaOffDiag      = Type(2.) / (Type(1.) + exp(Type(-5.)*logitSigmaOffDiag) ) - Type(1.);     
  // Scalars
  Type obj    = 0.0;   // objective function (neg log likelihood)
  Type pospen = 0.0;   // posfun penalty
  // Derived Parameters
  array<Type>   Ut(nS,nT);
  vector<Type>  msy = Bmsy * Umsy;

  // Procedure Section //
  // First create the correlation matrix for the species effects
  UNSTRUCTURED_CORR_t<Type> specEffCorr(SigmaOffDiag);
  vector<Type> SigmaDiag(nS);
  SigmaDiag = SigmaDiagMult*exp(lnSigmaDiag);    

  // Generate Population Dynamics
  array<Type> Bt(nS,nT);
  Bt.col(0) = Type(2)*Bmsy;
  for (int t=1; t<nT; t++)
  {
    Bt.col(t) = Bt.col(t-1) + Type(2)*exp(log(Bt.col(t-1))+lnUmsy) * (Type(1) - exp(log(Bt.col(t-1))-lnBmsy))/Type(2) - Ct.col(t-1);
    Bt.col(t) *= exp(eps_t(t-1) + zeta_st.col(t-1));
    for (int s=0; s<nS;s++)
    {
      Bt(s,t) = posfun(Bt(s,t), Type(1.), pospen);
      obj += 10*pospen;  
    } 
    // Add year effect contribution to objective function
    if (t == 1) obj += Type(0.5)*(lnkappa2 + pow(eps_t(t-1),2)/kappa2) + eps_t(t-1);
    else obj += Type(0.5)*(lnkappa2 + pow(eps_t(t-1) - gammaYr*eps_t(t-2),2)/kappa2) + eps_t(t-1);
    // Add correlated species effects contribution to likelihood
    if (nS > 1) obj += VECSCALE(specEffCorr,sqrt(SigmaDiag))(zeta_st.col(t-1));
  }

  // Compute observation likelihood
  for (int t=0; t<nT; t++)
  {
    for(int s=0; s<nS; s++)
    {
      if (It(s,t) > 0) 
      {
        Type res = log(It(s,t)) - log(Bt(s,t));
        obj += Type(0.5)*( lntau2 + pow(res-lnq(s),2)/tau2);
      }
    }
  }
  // Add priors
  // eqbm biomass
  for (int s=0; s<nS; s++ )
  {
    obj += pow(lnBmsy(s) - mlnBmsy(s),2)/s2lnBmsy(s);  
  }
  // multispecies shared priors
  if (nS > 1)
  {
    // 1st level priors
    for (int s=0; s<nS;s++)
    {
      // catchability
      obj += Type(0.5)*( lntauq2 + pow(lnq(s) - lnqbar,2)/tauq2);
      // productivity
      obj += Type(0.5)*( lnsigUmsy2 + pow(lnUmsy(s) - lnUmsybar,2)/sigUmsy2 );
    }  
    // Hyperpriors
    // catchability
    obj -= dnorm(lnqbar,mlnq,sqrt(s2lnq),true);
    // productivity
    obj -= dnorm(lnUmsybar,mlnUmsy,sqrt(s2lnUmsy),true);
    // End multispecies shared priors
  } else {
    // Now for single species model
    for (int s=0; s<nS; s++)
    {
      obj -= dnorm(lnq(s),mlnq,sqrt(s2lnq),true);
      // productivity
      obj -= dnorm(lnUmsy(s),mlnUmsy,sqrt(s2lnUmsy),true);  
    }
    
  }
  
  
  // Variance IG priors
  // Obs error var
  obj += (tau2IG(0)+Type(1))*lntau2+tau2IG(1)/tau2;
  // shared q prior variance
  obj += (tauq2IG(0)+Type(1))*lntauq2+tauq2IG(1)/exp(lntauq2);
  // shared U prior variance
  obj += (sigU2IG(0)+Type(1))*lnsigUmsy2+sigU2IG(1)/exp(lnsigUmsy2);
  // year effect deviations var
  obj += (kappa2IG(0)+Type(1))*lnkappa2+kappa2IG(1)/kappa2;
  // Specific effect variance
  obj += (Sigma2IG(0)+Type(1))*lnSigmaDiag+Sigma2IG(1)/exp(lnSigmaDiag);

  

  // Derive some output variables
  Ut  = Ct / Bt;
  vector<Type> DnT = Bt.col(nT-1)/Bmsy/2;

  // Reporting variables
  ADREPORT(Bt);
  ADREPORT(q);
  ADREPORT(Bmsy);
  ADREPORT(Umsy);
  ADREPORT(msy);
  ADREPORT(Umsybar);
  ADREPORT(qbar);
  ADREPORT(tau2);
  ADREPORT(kappa2);
  if (nS > 1 )
  {
    ADREPORT(SigmaDiag);    
    ADREPORT(SigmaOffDiag); 
  }
  ADREPORT(gammaYr);
  ADREPORT(Ut);

  // still haven't figured out how to use report
  // REPORT(specEffCorr.cov());

  return obj;
}



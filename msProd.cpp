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
  Type gammaYr      = Type(2.) / (Type(1.) + exp(Type(-2.)*logit_gammaYr) ) - Type(1.);
  vector<Type> SigmaOffDiag(nS*(nS-1)/2);  
  SigmaOffDiag      = Type(2.) / (Type(1.) + exp(Type(-2.)*logitSigmaOffDiag) ) - Type(1.);     
  // Scalars
  Type obj    = 0.0;   // objective function (neg log likelihood)
  Type pospen = 0.0;   // posfun penalty
  // Derived Parameters
  array<Type>   Ut(nS,nT);
  vector<Type>  msy = Bmsy * Umsy;

  /*procedure section*/
  // First, create correlation component of Sigma
  UNSTRUCTURED_CORR_t<Type> specEffCorr(SigmaOffDiag);
  vector<Type> SigmaDiag(nS);
  SigmaDiag = SigmaDiagMult*exp(lnSigmaDiag);    
  // Compute state values, add residues to likelihood
  for (int t=0;t<nT;t++)
  {
    vector<Type> Bpred(nS);
    if (t == 0) Bpred = Type(2) * Bmsy;
    else Bpred = exp(lnBt.col(t-1)) + Type(2)*Umsy*exp(lnBt.col(t-1))*(Type(1) - Type(0.5)*exp(lnBt.col(t-1))/Bmsy) - Ct.col(t-1);
    for(int s=0;s<nS;s++) 
    {
      Bpred(s) = posfun(Bpred(s), Ct(s,t), pospen);
      obj += 10*pospen;
    }
    vector<Type> epst(nS);
    epst.fill(eps_t(t));
    lnBt.col(t) = log(Bpred) + epst + zeta_st.col(t);
    // Add correlated species effects contribution to likelihood
    if (nS > 1) obj += VECSCALE(specEffCorr,sqrt(SigmaDiag))(zeta_st.col(t));
  }
  // Add AR1 year effect contribution to likelihood
  obj += SCALE(AR1(gammaYr),sqrt(kappa2))(eps_t);

  // Compute observation likelihood
  // First create a diagonal covariance matrix for the observations
  matrix<Type> Tau(nS,nS);
  Tau.fill(0.0);
  Tau.diagonal() = tau2 * tau2mult;
  // Now create MV normal distribution with Tau as cov mtx
  MVNORM_t<Type> obs_nll(Tau);
  for (int t = 0; t< nT; t++)
  {
    obj += obs_nll(log(It.col(t)) - (lnq + lnBt.col(t)));
  }
  // Shared priors
  for (int s=0; s<nS;s++)
  {
    // catchability
    obj -= dnorm(lnq(s),lnqbar,sqrt(tauq2),true);
    // productivity
    obj -= dnorm(lnUmsy(s),lnUmsybar,sqrt(sigUmsy2),true);
    // equlibrium biomass
    obj -= dnorm(lnBmsy(s),mlnBmsy(s),sqrt(s2lnBmsy(s)),true);
  }
  // Hyperpriors
  // catchability
  obj -= dnorm(lnqbar,mlnq,sqrt(s2lnq),true);
  // productivity
  obj -= dnorm(lnUmsybar,mlnUmsy,sqrt(s2lnUmsy),true);

  array<Type> Bt(nS,nT);
  Bt = exp(lnBt);
  Ut = Ct / Bt;
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



// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// msProd_TMB.cpp
// 
// TMB version of the ADMB model in msProd.tpl
// 
// Author: Samuel Johnson
// Date: August 25, 2014                    
// Purpose: To fit a multi-species RH production model to commercial catch
// and survey CPUE data.
// 
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

#include <TMB.hpp>                                // Links in the TMB libraries
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
  /*data section*/
  DATA_ARRAY(It);               // CPUE data
  DATA_ARRAY(Ct);               // Catch data

  /*parameter section*/
  // Leading Parameters
  PARAMETER_VECTOR(lnBmsy);             // Biomass at MSY
  PARAMETER_VECTOR(lnUmsy);             // Optimal exploitation rate
  PARAMETER_VECTOR(lntau2);             // Species specific Obs error variance
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
  // Random Effects
  PARAMETER_VECTOR(eps_t);              // year effect
  PARAMETER(lnkappa2);                  // year effect uncorrelated variance
  PARAMETER(logit_gammaYr);             // AR1 auto-corr on year effect
  PARAMETER_ARRAY(zeta_st);             // species effect
  PARAMETER_VECTOR(lnSigmaDiag);        // Species effect cov matrix diag
  PARAMETER_VECTOR(lnSigmaOffDiag);     // Species effect cov chol off diag

  // derived vars //
  // indexing variables
  int nS=It.dim(0);
  int nT=It.dim(1);
  // State variables
  array<Type> lnBt(nS,nT);
  // Leading parameters
  vector<Type> Bmsy = exp(lnBmsy);
  vector<Type> Umsy = exp(lnUmsy);
  vector<Type> tau2 = exp(lntau2);
  vector<Type> q    = exp(lnq);
  // Prior hyperpars
  Type tauq2    = exp(lntauq2);
  Type sigUmsy2 = exp(lnsigUmsy2);
  // Random Effects
  Type kappa2  = exp(lnkappa2);
  Type gammaYr = invLogit(logit_gammaYr,Type(2.),Type(1.));
  //matrix<Type> Sigma(nS,nS);
  // Scalars
  Type obj    = 0.0;   // objective function (neg log likelihood)
  Type pospen = 0.0;   // posfun penalty

  /*procedure section*/
  // First, create correlation component of Sigma
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> Sigma(exp(lnSigmaOffDiag));

  // Compute state values, add residues to likelihood
  for (int t=0;t<nT;t++)
  {
    vector<Type> Bpred(nS);
    if (t == 0) Bpred = Type(2) * Bmsy;
    else Bpred = exp(lnBt.col(t-1)) + Type(2)*Umsy*exp(lnBt.col(t-1))*(Type(1) - Type(0.5)*exp(lnBt.col(t-1))/Bmsy) - Ct.col(t-1);
    for(int s=0;s<nS;s++) 
    {
      Bpred(s) = posfun(Bpred(s), Ct(s,t) + 1, pospen);
      obj += pospen;
    }
    vector<Type> epst(nS);
    epst.fill(eps_t(t));
    lnBt.col(t) = log(Bpred) + epst + zeta_st.col(t);
    // Add species effect contribution to likelihood
    obj += VECSCALE(Sigma,exp(Type(0.5)*lnSigmaDiag))(zeta_st.col(t));
  }
  // Add year effect contribution to likelihood
  obj += SCALE(AR1(gammaYr),sqrt(kappa2))(eps_t);

  // Compute observation likelihood
  // First create a diagonal covariance matrix for the observations
  matrix<Type> Tau(nS,nS);
  Tau.fill(0.0);
  Tau.diagonal() = tau2;
  // Now create MV normal distribution with Tau as cov mtx
  MVNORM_t<Type> obs_nll(Tau);
  for (int t = 0; t< nT; t++)
  {
    // This should actually be an MVNORM call, let's try with the vector
    obj += obsnll(log(It.col(t)) - (lnq + lnBt.col(t)));
  }
  // Shared priors
  for (int s=0; s<nS;s++)
  {
    // catchability
    obj -= dnorm(lnq(s),lnqbar,sqrt(tauq2),true);
    // productivity
    obj -= dnorm(lnUmsy(s),lnUmsybar,sqrt(sigUmsy2),true);
  }
  // Hyperpriors
  // catchability
  obj -= dnorm(lnqbar,mlnq,sqrt(s2lnq),true);
  // productivity
  obj -= dnorm(lnUmsybar,mlnUmsy,sqrt(s2lnUmsy),true);

  vector<Type> Bt = exp(lnBt);
  // Reporting variables
  ADREPORT(Bmsy);
  ADREPORT(Umsy);
  ADREPORT(q);
  ADREPORT(tau2);
  ADREPORT(kappa2);
  REPORT(Sigma.cov());
  ADREPORT(lnSigmaDiag);
  ADREPORT(lnSigmaOffDiag);
  ADREPORT(Bt);

  return obj;
}

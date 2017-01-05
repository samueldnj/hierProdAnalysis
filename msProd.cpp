  #include <admodel.h>
  
  int mcHeader=0;
  int mcTrials=0;
  ofstream mcoutpar("msProdParMC.dat");
  ofstream mcoutbio("msProdBioMC.dat");
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <msProd.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nT.allocate("nT");
  nS.allocate("nS");
  katch.allocate(1,nS,1,nT,"katch");
  It.allocate(1,nS,1,nT,"It");
  phz_Bmsy.allocate("phz_Bmsy");
  phz_Umsy.allocate("phz_Umsy");
  phz_tau.allocate("phz_tau");
  phz_kappa.allocate("phz_kappa");
  phz_Sigma.allocate("phz_Sigma");
  phz_lnq.allocate("phz_lnq");
  phz_mlnq.allocate("phz_mlnq");
  phz_s2lnq.allocate("phz_s2lnq");
  phz_varPriors.allocate("phz_varPriors");
  phz_omegat.allocate("phz_omegat");
  phz_zetat.allocate("phz_zetat");
  phz_AR.allocate("phz_AR");
  phz_chol.allocate("phz_chol");
  dumm.allocate("dumm");
 phz_tau2qs = phz_lnq;
    if(dumm!=999)
    {
      cout<<"Error reading data.\n Fix it!! \n dumm = " << dumm << endl;
      ad_exit(1);
    }
    cEntries=nS*(nS-1)/2;
    verbose=0;
    if (nS == 1)
    {
      // if only one species, estimate only the shared REs
      phz_omegat  = phz_zetat;
      phz_kappa   = phz_Sigma;
      phz_zetat   = -1;
      phz_chol    = -1;
      phz_Sigma   = -1;
      phz_tau2qs  = -1;
    }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  lnBmsy.allocate(1,nS,phz_Bmsy,"lnBmsy");
  lnUmsy.allocate(1,nS,phz_Umsy,"lnUmsy");
  Bmsy.allocate(1,nS,"Bmsy");
  #ifndef NO_AD_INITIALIZE
    Bmsy.initialize();
  #endif
  Umsy.allocate(1,nS,"Umsy");
  #ifndef NO_AD_INITIALIZE
    Umsy.initialize();
  #endif
  lnTau2.allocate(phz_tau,"lnTau2");
  lnkappa2.allocate(phz_kappa,"lnkappa2");
  lnSigma2.allocate(phz_Sigma,"lnSigma2");
  tau2.allocate("tau2");
  #ifndef NO_AD_INITIALIZE
  tau2.initialize();
  #endif
  kappa2.allocate("kappa2");
  #ifndef NO_AD_INITIALIZE
  kappa2.initialize();
  #endif
  Sigma2.allocate("Sigma2");
  #ifndef NO_AD_INITIALIZE
  Sigma2.initialize();
  #endif
  tau.allocate("tau");
  #ifndef NO_AD_INITIALIZE
  tau.initialize();
  #endif
  kappa.allocate("kappa");
  #ifndef NO_AD_INITIALIZE
  kappa.initialize();
  #endif
  Sigma.allocate("Sigma");
  #ifndef NO_AD_INITIALIZE
  Sigma.initialize();
  #endif
  gamma.allocate(0,0.99,phz_AR,"gamma");
  c.allocate(1,cEntries,-1,1,phz_chol,"c");
  chol.allocate(1,nS,1,nS,"chol");
  #ifndef NO_AD_INITIALIZE
    chol.initialize();
  #endif
  mBmsy.allocate(1,nS,-1,"mBmsy");
  sBmsy.allocate(1,nS,-1,"sBmsy");
  mUmsy.allocate(1,nS,-1,"mUmsy");
  sUmsy.allocate(1,nS,-1,"sUmsy");
  alpha_kappa.allocate(1,10,phz_varPriors,"alpha_kappa");
  beta_kappa.allocate(0.01,4,phz_varPriors,"beta_kappa");
  alpha_Sigma.allocate(1,10,phz_varPriors,"alpha_Sigma");
  beta_Sigma.allocate(0.01,4,phz_varPriors,"beta_Sigma");
  alpha_tau.allocate(1,10,phz_varPriors,"alpha_tau");
  beta_tau.allocate(0.01,4,phz_varPriors,"beta_tau");
  mlnq.allocate(-3.,3.,phz_mlnq,"mlnq");
  s2lnq.allocate(-1,"s2lnq");
  lnqbar.allocate(-5,2,phz_lnq,"lnqbar");
  lnTau2_qs.allocate(-5,2,phz_tau2qs,"lnTau2_qs");
  omegat.allocate(1,nT,-3.,3.,phz_omegat,"omegat");
  zetat.allocate(1,nS,1,nT,phz_zetat,"zetat");
  epst.allocate(1,nT,"epst");
  #ifndef NO_AD_INITIALIZE
    epst.initialize();
  #endif
  zetatc.allocate(1,nS,1,nT,"zetatc");
  #ifndef NO_AD_INITIALIZE
    zetatc.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  lnqhatSpec.allocate(1,nS,"lnqhatSpec");
  #ifndef NO_AD_INITIALIZE
    lnqhatSpec.initialize();
  #endif
  q.allocate(1,nS,"q");
  Bt.allocate(1,nS,1,nT,"Bt");
  #ifndef NO_AD_INITIALIZE
    Bt.initialize();
  #endif
  msy.allocate(1,nS,"msy");
  #ifndef NO_AD_INITIALIZE
    msy.initialize();
  #endif
  UnT_bar.allocate(1,nS,"UnT_bar");
  #ifndef NO_AD_INITIALIZE
    UnT_bar.initialize();
  #endif
  dep_bar.allocate(1,nS,"dep_bar");
  #ifndef NO_AD_INITIALIZE
    dep_bar.initialize();
  #endif
  z.allocate(1,nS,1,nT,"z");
  #ifndef NO_AD_INITIALIZE
    z.initialize();
  #endif
  zSum.allocate(1,nS,"zSum");
  #ifndef NO_AD_INITIALIZE
    zSum.initialize();
  #endif
  ss.allocate(1,nS,"ss");
  #ifndef NO_AD_INITIALIZE
    ss.initialize();
  #endif
  validObs.allocate(1,nS,"validObs");
  #ifndef NO_AD_INITIALIZE
    validObs.initialize();
  #endif
  obsLike.allocate(1,nS,"obsLike");
  #ifndef NO_AD_INITIALIZE
    obsLike.initialize();
  #endif
  procLike.allocate(1,nS,"procLike");
  #ifndef NO_AD_INITIALIZE
    procLike.initialize();
  #endif
  totalLike.allocate("totalLike");
  #ifndef NO_AD_INITIALIZE
  totalLike.initialize();
  #endif
  vCov.allocate(1,nS,1,nS,"vCov");
  #ifndef NO_AD_INITIALIZE
    vCov.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  pospen.allocate("pospen");
  #ifndef NO_AD_INITIALIZE
  pospen.initialize();
  #endif
  totalPrior.allocate("totalPrior");
  #ifndef NO_AD_INITIALIZE
  totalPrior.initialize();
  #endif
  BmsyPrior.allocate(1,nS,"BmsyPrior");
  #ifndef NO_AD_INITIALIZE
    BmsyPrior.initialize();
  #endif
  UmsyPrior.allocate(1,nS,"UmsyPrior");
  #ifndef NO_AD_INITIALIZE
    UmsyPrior.initialize();
  #endif
  kappaPrior.allocate("kappaPrior");
  #ifndef NO_AD_INITIALIZE
  kappaPrior.initialize();
  #endif
  tauPrior.allocate("tauPrior");
  #ifndef NO_AD_INITIALIZE
  tauPrior.initialize();
  #endif
  SigmaPrior.allocate("SigmaPrior");
  #ifndef NO_AD_INITIALIZE
  SigmaPrior.initialize();
  #endif
  lnqPrior.allocate(1,nS,"lnqPrior");
  #ifndef NO_AD_INITIALIZE
    lnqPrior.initialize();
  #endif
}

void model_parameters::userfunction(void)
{
  f =0.0;
  // Initialise obj function var
  f = 0.;
  // initialise derived variables
  dep_bar = 0.;
  UnT_bar = 0.;
  // Run state dynamics function
  stateDynamics();
  if (verbose) cout << "State Dynamic complete" << endl;
  // apply observation model, compute conditional value for lnq
  obsModel();
  if (verbose) cout << "Observation model complete" << endl;
  // Compute likelihood for dynamics
  calcLikelihoods();
  if (verbose) cout << "Likelihood computation complete" << endl;
  //cout << "nll = " << nll << endl;
  // Compute priors
  calcPriors();
  if (verbose)  cout << "Prior computation complete" << endl;
  // cout << "prior = " << prior << endl;
  // cout << "lnmsy = " << lnmsy << endl;
  // cout << "lnFmsy = " << lnFmsy << endl;
  // cout << "epst = " << epst << endl;
  // Output MCMC data if running mcmc trials
  if(mceval_phase())
  { 
    mcDumpOut();
  }
  // compute objective function value
  f += totalLike;              // Likelihood and proc error penalties
  f += totalPrior;             // add total of prior densities
  f += fpen;                   // positive function value penalties
}

void model_parameters::stateDynamics(void)
{
  // exponentiate leading parameters
  Umsy  = mfexp ( lnUmsy );         // optimal exploitation rate
  //cout << "Umsy = " << Umsy << endl;
  Bmsy  = mfexp ( lnBmsy );         // optimal eqbm biomass
  //cout << "Bmsy = " << Bmsy << endl;
  msy   = elem_prod( Bmsy, Umsy);   // optimal eqbm harvest
  //cout << "msy = " << msy << endl;
  // initialise penalisers
  fpen = 0.; pospen = 0.;
  // initialise RE variables
  // omegat.initialize();
  // zetat.initialize();
  chol.initialize();
  epst.initialize();
  zetatc.initialize();
  tau2.initialize();
  // variance pars //
  // Obs error variance
  tau2    = mfexp(lnTau2);
  tau     = mfexp(lnTau2*0.5);
  //cout << "lnTau2 = " << lnTau2 << endl;
  // shared effects
  if (phz_kappa > 0) 
  {
    kappa2  = mfexp(lnkappa2);
    kappa   = mfexp(lnkappa2*0.5);
    //cout << "lnkappa2 = " << lnkappa2 << endl;
  }
  //cout << "SS variances are good" << endl;
  // Species specific effects
  if (nS > 1 & phz_zetat > 0 ) 
  { 
    Sigma2  = mfexp(lnSigma2);
    Sigma   = mfexp(0.5*lnSigma2);
    //cout << "lnSigma2 = " << lnSigma2 << endl; 
    // Restrict shared effect to have mean 0
    //omegat -= mean(omegat);
    // Construct Cholesky factor of correlation mtx
    chol(1,1) = 1.;
    int k=1;
    for (int i=2;i<=nS;i++)
    {
      chol(i,i)=1.;
      for(int j=1;j<=i-1;j++)
      {
        chol(i,j)=c(k++);
      }
      chol(i)(1,i) /= norm(chol(i)(1,i));
    }
    // Compute correlated zetat values
    zetatc = chol * zetat;
  }
  //cout << "Starting Pop Dyn" << endl;
  // Run a loop for population dynamics of each stock
  // Assume stocks are at equilibrium initially ( B1 = B0e^delta1 )
  for (int s=1;s<=nS;s++)
  {
    // start AR(1) process
    epst(1) = omegat(1);
    Bt(s,1) = ( 2. * Bmsy(s) ) * mfexp(epst(1));
    if (phz_zetat > 0) Bt(s,1) *= mfexp(zetatc(s,1));
    // loop over time periods to build dynamic model for t > 1
    for(int t=1; t<nT; t++)
    { 
      // continue AR(1) process
      epst(t+1) = gamma*epst(t) + omegat(t+1);
      pospen = 0.;
      // Run state dynamics equation, add process error
      Bt(s,t+1) = ( Bt(s,t) + 
                    2. * Umsy(s) * Bt(s,t) * (1 - Bt(s,t)/Bmsy(s)/2.0 ) - 
                    katch(s,t) );
      Bt(s,t+1) *= mfexp ( epst(t+1) );
      if (phz_zetat > 0) Bt(s,t+1) *= mfexp(zetatc(s,t+1));
      Bt(s,t+1) = posfun ( Bt(s,t+1), katch(s,t+1)+10, pospen );
      // Increment function penaliser variable
      fpen += 1000. * pospen;
    }
  }
  //cout << "Bt = " <<  Bt <<  endl;
  //cout << "epst = " <<  epst <<  endl;
  // compute derived management quantities
  dep_bar = elem_div(column(Bt,nT),Bmsy) / 2;      // depletion estimate
  UnT_bar = elem_div(column(katch,nT),column(Bt,nT)); // Most recent exp rate
}

void model_parameters::obsModel(void)
{
  zSum=0.0;
  validObs.initialize();
  lnqhatSpec.initialize();
  ss.initialize();
  //cout << "Inside Observation Model" << endl;
  for ( int s=1;s<=nS;s++)
  {
    //cout << "Starting obs model for species " << s << endl;
    for ( int t = 1; t<=nT;t++)
    {
      //cout << "Species " << s << " time " << t << endl;
      if (It(s,t) > 0)
      {
        z(s,t) = log(It(s,t)) - log(Bt(s,t));
        zSum(s) += z(s,t);
        validObs(s) += 1;
      }
    }
    //cout << "z = " << z(s) << endl;
    // Compute conditional posterior mode of lnq_hat (BDA3 eq 11.10)
    if (s == 1) lnqhatSpec(s) = (mlnq/s2lnq + zSum(s)/tau2)/(1/s2lnq + validObs(s)/tau2);
    if (s > 1)
    {
      lnqhatSpec(s) = (zSum(s)/tau2 + lnqbar/mfexp(lnTau2_qs))/(validObs(s)/tau2 + 1/mfexp(lnTau2_qs));
    }
    //lnqhat(s) = zSum(s)/validObs(s); // mean residual for non-hierarchical model
    //cout << "lnqhat = " << lnqhat(s) << endl;
    for (int t=1;t<=nT;t++)
    {
      //cout << "Computing squared residual for time " << t << endl;
      if(It(s,t)>0.0)
      {
        ss(s) += pow (z(s,t) - lnqhatSpec(s),2.0);
      }
    }
    //cout << "Species " << s << " obs model complete" << endl;
  }
}

void model_parameters::calcLikelihoods(void)
{
  // Initialise NLL variables
  totalLike.initialize();
  obsLike.initialize(); 
  //omegaLike.initialize();
  procLike.initialize();
  // compute observation likelihood
  obsLike = 0.5*validObs*log(tau2) + 0.5*ss/tau2;
  totalLike += sum(obsLike);
  // Shared effect estimated
  if ((phz_omegat > 0) ) //& (phz_zetat<0) )
  {
    procLike(1) = 0.5*nT*lnkappa2 + 0.5*norm2(omegat)/(kappa2) + sum(omegat);
    totalLike += sum(procLike);
  }
  // species specific effect estimated
  if ( (phz_zetat > 0) ) //& (phz_omegat < 0 ) )
  {
    for (int s=1; s<=nS; s++)
    {
      procLike(s) = 0.5*nT*lnSigma2 + 0.5*norm2(zetat(s))/(Sigma2) + sum(zetat(s));
    }
    totalLike += sum(procLike);
  }
}

void model_parameters::calcPriors(void)
{
  // Initialise prior vars
  BmsyPrior.initialize();   // optimal B
  UmsyPrior.initialize();   // optimal U
  lnqPrior.initialize();    // shared catchability prior (vector valued)
  kappaPrior.initialize();  // shared RE variance prior
  SigmaPrior.initialize();  // species specific RE variance prior (vector valued)
  tauPrior.initialize();    // obs error variance prior (vector valued)
  totalPrior.initialize();  // sum of priors
  // First, the priors that are calculated in every phase
  // msy
  if ( phz_Bmsy > 0)
  {
    BmsyPrior = elem_div(pow ( Bmsy - mBmsy, 2 ), pow(sBmsy,2))*0.5;  
  }
  // Then Fmsy
  if (phz_Umsy > 0)
  {
    UmsyPrior = elem_div( pow ( Umsy - mUmsy, 2 ), pow(sUmsy,2)) * 0.5;  
  }
  // lnq prior (shared)
  if (phz_mlnq > 0)
  {
    if (nS == 1) lnqPrior = pow ( lnqhatSpec - mlnq, 2 ) / s2lnq /2;
    else {
      lnqPrior = 0.5*log(mfexp(lnTau2_qs)) + 0.5*norm2(lnqhatSpec - lnqbar)/mfexp(lnTau2_qs);
      if ( value(lnqbar) != value(mlnq) ) 
      {
        lnqPrior += 0.5*pow(lnqbar - mlnq,2.)/s2lnq/nS;
      }
    }
  }
  // total the schaefer model priors
  totalPrior = sum(BmsyPrior) + sum(UmsyPrior) + sum ( lnqPrior );
  if ( verbose ) cout << "Added total Schaefer Priors" << endl;
  // Now the shared effects variance (IG(alpha,beta)) if lnkappa2 is active
  if (active(lnkappa2))
  {
    kappaPrior = (alpha_kappa+1)*lnkappa2+beta_kappa/kappa2;
    totalPrior += kappaPrior;
  }
  // species specific effects variance
  if (active(lnSigma2))
  {
    SigmaPrior = alpha_Sigma+1*lnSigma2 + beta_Sigma/Sigma2;
    totalPrior += SigmaPrior;
  }
  // tau2 prior if estimated
  if (active(lnTau2))
  {
    tauPrior = (alpha_tau+1)*log(tau2)+beta_tau/tau2;
    totalPrior += tauPrior;
  }
  // Apply a Jeffries prior to msy
  //totalPrior += sum(1/msy);
  // If estimating IG parameters, apply a jeffreys prior to the
  // beta pars
  if ( active(beta_kappa) )
  {
    totalPrior += beta_kappa;
  }
  if (active(beta_Sigma))
  {
   totalPrior += beta_Sigma; 
  }
  if (active(beta_tau) )
  {
    totalPrior += beta_tau;
  }
}

void model_parameters::mcDumpOut(void)
{
  if( mcHeader == 0 )
  {
    // Output the parameters needed for performance testing of SA method.
    mcoutpar << "Bmsy msy Umsy q Sigma2 kappa2 tau2 UnT_bar BnT dep_bar " << endl;
    for ( int s=1;s<=nS;s++)
    {
      // Output the parameter values for this replicate
      mcoutpar << Bmsy(s) << " " << msy(s) <<" "<< Umsy(s) <<" "<< mfexp(lnqhatSpec(s));
      mcoutpar <<" "<< Sigma2 <<" "<< kappa2 <<" " << tau2 <<" " << UnT_bar(s) << " ";
      mcoutpar << Bt(s,nT) <<" "<< dep_bar(s) << endl;
      // Output the biomass estimates for this MC evaluation
      mcoutbio << Bt(s) << endl;
    }
    mcHeader = 1;
  }
  else
  {
    // loop over species
    for ( int s=1;s<=nS;s++)
    {
      // Output the parameter values for this replicate
      mcoutpar << Bmsy(s) << " " << msy(s) <<" "<< Umsy(s) <<" "<< mfexp(lnqhatSpec(s));
      mcoutpar <<" "<< Sigma2 <<" "<< kappa2 <<" " << tau2 <<" " << UnT_bar(s) << " ";
      mcoutpar << Bt(s,nT) <<" "<< dep_bar(s) << endl;
      // Output the biomass estimates for this MC evaluation
      mcoutbio << Bt(s) << endl;
    }
  }
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  // Output estimates and derived variables for comparison to true
  // values in sim-est experiments
  report << "## msProd Results" << endl;
  report << "## Parameter estimates (leading and concentrated) " << endl;
  report << "# Bmsy" << endl;
  report << Bmsy << endl;
  report << "# Umsy" << endl;
  report << Umsy << endl;
  report << "# kappa2" << endl;
  report << kappa2 << endl;
  report << "# Sigma2" << endl;
  report << Sigma2 << endl;
  report << "# tau2" << endl;
  report << tau2 << endl;
  report << endl;
  report << "## Scaled and correlated REs" << endl;
  report << "# epst" << endl;
  report << epst << endl;
  report << "# zetat" << endl;
  report << zetatc << endl;
  report << endl;
  report << "## Unscaled REs and correlation factors" << endl;
  report << "# omegat" << endl;
  report << omegat << endl;
  report << "# zetat" << endl;
  report << zetat << endl;
  report << "# gamma" << endl;
  report << gamma << endl;
  if (nS > 1)
  {
    report << "# c" << endl;
    report << c << endl;
    report << "# chol" << endl;
    report << chol << endl;
  }
  report << endl;
  report << "## Derived variables" << endl;
  report << "# msy" << endl;
  report << msy << endl;
  report <<"# q" << endl;
  report << mfexp(lnqhatSpec) <<endl;
  if (nS > 1)
  {
    report <<"# qbar" << endl;
    report << mfexp(lnqbar) <<endl;  
  }
  report << "# D" << endl;
  report << dep_bar << endl;
  report << "# UnT" << endl;
  report << UnT_bar << endl;
  report << "# lastCt" << endl;
  report << column(katch,nT) << endl;
  report << "# Bt" << endl;
  report << Bt << endl;
  report << "# Ut" << endl;
  report << elem_div(katch,Bt) << endl;
  report << endl;
  report << "## Prior Hyperparameters" << endl;
  report << "# mBmsy" << endl;  
  report << mBmsy << endl;
  report << "# sBmsy" << endl;  
  report << sBmsy << endl;
  report << "# mUmsy" << endl;  
  report << mUmsy << endl;
  report << "# sUmsy" << endl;  
  report << sUmsy << endl;
  report << "# alpha_kappa" << endl;
  report << alpha_kappa << endl;
  report << "# beta_kappa" << endl;
  report << beta_kappa << endl;
  report << "# alpha_tau" << endl;
  report << alpha_tau << endl;
  report << "# beta_tau" << endl;
  report << beta_tau << endl;
  report << "# mlnq" << endl;  
  report << mlnq << endl;
  report << "# slnq" << endl;  
  report << s2lnq << endl;
  report << endl;
  report << "## Data" << endl;
  report << "# nT" << endl;
  report << nT << endl;
  report << "# nS" << endl;
  report << nS << endl;
  report << "# Ct" << endl;
  report << katch << endl;
  report << "# It" << endl;
  report << It << endl;
  report << endl;
  report << "## Minimization properties" << endl;
  report << "# total_likelihood" << endl;
  report << totalLike << endl;
  report << "# obsLike" << endl;
  report << obsLike << endl;
  report << "# validObs" << endl;
  report << validObs << endl;
  report << "# ss" << endl;
  report << ss << endl;
  //report << "# omegaLike" << endl;
  //report << omegaLike << endl;
  report << "# procLike" << endl;
  report << procLike << endl;
  report << "# totalPrior" << endl;
  report << totalPrior << endl;
  report << "# BmsyPrior" << endl;
  report << BmsyPrior << endl;
  report << "# UmsyPrior" << endl;
  report << UmsyPrior << endl;
  report << "# lnqPrior" << endl;
  report << lnqPrior << endl;
  report << "# kappaPrior" << endl;
  report << kappaPrior << endl;
  report << "# SigmaPrior" << endl;
  report << SigmaPrior << endl;
  report << "# tauPrior" << endl;
  report << tauPrior << endl;
  report << "# objFun" << endl;
  report << *objective_function_value::pobjfun << endl;
  report << "# maxGrad" << endl;
  report << objective_function_value::gmax << endl;
  report << "# fpen" << endl;
  report << fpen << endl;
  report << "# iExit" << endl;
  report << iexit << endl;
  report << endl;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 50000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
  //gradient_structure::set_MAX_NVAR_OFFSET(5000);
  //gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}

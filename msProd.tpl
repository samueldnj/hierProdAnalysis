// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// msProd.tpl
// 
// Multi species surplus production model with logistic stock-recruit function,
// detects covariance among multiple species.
// 
// Author: Samuel Johnson
// Date: August 25, 2014                    
// Purpose: To fit a multi-species RH production model to commercial catch
// and survey CPUE data.
// 
// ToDo:  1. Missing data handling
//       
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

DATA_SECTION
  // Control parameters for data
  init_int nT;                        // number of time iterations
  init_int nS;                        // number of species
  init_matrix katch(1,nS,1,nT);       // catch in kg
  init_matrix It(1,nS,1,nT);          // indices of abundance

  // Estimation phases
  init_int phz_Bmsy;            // maximum sustainable yield
  init_int phz_Umsy;            // explotation rate
  init_int phz_tau;             // tau - obs error variance
  init_int phz_kappa;           // shared RE variance
  init_int phz_Sigma;           // species specific RE variance
  init_int phz_lnq;             // q - survey catchability - concentrated
  init_int phz_tau2q;           // shared q prior variance
  init_int phz_mlnq;            // lnqbar prior mean
  init_int phz_s2lnq;           // lnqbar prior variance
  init_int phz_varPriors;       // IG priors on variance terms
  init_int phz_omegat;          // shared envrionmental REs
  init_int phz_zetat;           // species spec REs
  init_int phz_AR;              // AR (gamma) for shared REs (omega)
  init_int phz_chol;            // cross correlation for zeta

  // dummy variable to check data read
  init_int dumm;                      // dummy variable
  
  // Constant for number of free pars in chol matx
  int cEntries;
  int verbose;

  // Procedure to exit if not all data is read correctly
  LOC_CALCS
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
      phz_tau2q  = -1;
    }

PARAMETER_SECTION
  // leading biological pars
  init_vector lnBmsy(1,nS,phz_Bmsy);      // Bmsy - log scale
  init_vector lnUmsy(1,nS,phz_Umsy);      // Umsy - log scale
  vector Bmsy(1,nS);        
  vector Umsy(1,nS);

  // log-variance parameters
  init_number lntau2(phz_tau);          // shared osb error var
  init_number lnkappa2(phz_kappa);      // year effect variance
  init_number lnSigma2(phz_Sigma);      // species effect var

  // normal scale variance
  number tau2;
  number kappa2;
  number Sigma2;
  
  // sd derived pars (transformed from variances)
  number tau;
  number kappa;
  number Sigma;

  // Covariance parameters
  init_bounded_number gamma(0,0.99,phz_AR);         // auto-regression factor
  init_bounded_vector c(1,cEntries,-1,1,phz_chol);  // cholesky factor off-diag entries
  matrix chol(1,nS,1,nS);                           // Matrix to hold cholesky factor

  // Prior hyperparameters
  // Biomass
  init_vector mBmsy(1,nS,-1);       // msy prior mean (species spec)
  init_vector sBmsy(1,nS,-1);       // msy prior SD (species spec)
  // Productivity
  init_vector mUmsy(1,nS,-1);       // lnUmsy prior mean (species spec)
  init_vector sUmsy(1,nS,-1);       // lnUmsy priod sd (species spec)
  // Variance
  init_bounded_number alpha_kappa(1,10,phz_varPriors);    // kappa shape (shared)
  init_bounded_number beta_kappa(0.01,4,phz_varPriors);    // kappa scale (shared)
  init_bounded_number alpha_Sigma(1,10,phz_varPriors); // Sigma prior shape (shared)
  init_bounded_number beta_Sigma(0.01,4,phz_varPriors); // Sigma prior scale (shared)
  init_bounded_number alpha_tau(1,10,phz_varPriors); // tau prior shape (shared)
  init_bounded_number beta_tau(0.01,4,phz_varPriors); // tau prior scale (shared)
  // Catchability
  init_bounded_number mlnq(-1.,1.,phz_mlnq);      // logq prior mean (shared)
  init_number s2lnq(phz_s2lnq);                          // logq prior sd (shared)
  init_bounded_number lnqbar(-5,2,phz_lnq);       // mean lnq across species
  init_bounded_number lntau2q(-5,2,phz_tau2q); // logqs prior var

  // process error deviations
  init_bounded_vector omegat(1,nT,-3.,3.,phz_omegat);   // year effect
  init_matrix zetat(1,nS,1,nT,phz_zetat);    // proc error
  
  // might be able to use dvector later for the following
  vector epst(1,nT);        // auto-correlated year effect
  matrix zetatc(1,nS,1,nT); // correlated zetat values (using chol)

  //objective function value
  objective_function_value f;

  // variables for conditional estimates of lnq(s)
  vector lnq_s(1,nS);
  sdreport_vector q(1,nS);

  // variables to hold derived values
  matrix Bt(1,nS,1,nT);       // predicted biomass
  vector msy(1,nS);           // msy
  vector UnT_bar(1,nS);       // estimated exploitation rate at nT 
  vector dep_bar(1,nS);       // depletion estimate

  // likelihood quantities
  matrix z(1,nS,1,nT);      // residues for each observation
  vector zSum(1,nS);        // sum of residues for each species
  vector ss(1,nS);          // sum of squared residuals for each species
  vector validObs(1,nS);    // counter for valid observations for each species
  vector obsLike(1,nS);     // NLL for observation error
  vector procLike(1,nS);    // NLL for REs
  number totalLike;         // total nll (sum of 3 likelihoods)
  matrix vCov(1,nS,1,nS);   // covariance matrix for combined effects
  
  // penalisers
  number fpen;               // penalty to obj function for bad dynamics
  number pospen;             // penalty to add to fpen if biomass dips

  // Prior quantities
  number totalPrior;        // total sum of priors
  vector BmsyPrior(1,nS);   // Bmsy normal prior
  vector UmsyPrior(1,nS);   // Umsy normal prior
  number kappaPrior;        // total RE variance prior
  number tauPrior;    // observation error variance prior
  number SigmaPrior;        //
  vector lnqPrior(1,nS);    // shared lnq prior density


PROCEDURE_SECTION
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

// State dynamics subroutine
FUNCTION stateDynamics
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
  tau2    = mfexp(lntau2);
  tau     = mfexp(lntau2*0.5);
  //cout << "lntau2 = " << lntau2 << endl;
  
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
      fpen += 100. * pospen;
    }
  }
  //cout << "Bt = " <<  Bt <<  endl;
  //cout << "epst = " <<  epst <<  endl;
  
  // compute derived management quantities
  dep_bar = elem_div(column(Bt,nT),Bmsy) / 2;      // depletion estimate
  UnT_bar = elem_div(column(katch,nT),column(Bt,nT)); // Most recent exp rate
  
// function to compute predicted observations and residuals
FUNCTION obsModel
  zSum=0.0;
  validObs.initialize();
  lnq_s.initialize();
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
    if (s == 1) lnq_s(s) = (mlnq/s2lnq + zSum(s)/tau2)/(1/s2lnq + validObs(s)/tau2);
    if (s > 1)
    {
      lnq_s(s) = (zSum(s)/tau2 + lnqbar/mfexp(lntau2q))/(validObs(s)/tau2 + 1/mfexp(lntau2q));
    }
    //lnqhat(s) = zSum(s)/validObs(s); // mean residual for non-hierarchical model
    //cout << "lnqhat = " << lnqhat(s) << endl;
    for (int t=1;t<=nT;t++)
    {
      //cout << "Computing squared residual for time " << t << endl;
      if(It(s,t)>0.0)
      {
        ss(s) += pow (z(s,t) - lnq_s(s),2.0);
      }
    }
    //cout << "Species " << s << " obs model complete" << endl;
  }

// Subroutine to compute negative log likelihoods for obs and proc errors
FUNCTION calcLikelihoods
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

// Subroutine to compute prior densities to add to posterior
FUNCTION calcPriors

  // Initialise prior vars
  BmsyPrior.initialize();   // optimal B
  UmsyPrior.initialize();   // optimal U
  lnqPrior.initialize();    // shared catchability prior (vector valued)
  kappaPrior.initialize();  // shared RE variance prior
  SigmaPrior.initialize();  // species specific RE variance prior (vector valued)
  tauPrior.initialize();    // obs error variance prior (vector valued)
  totalPrior.initialize();  // sum of priors

  // First, the Schaefer mode priors
  // msy
  if ( phz_Bmsy > 0)
  {
    BmsyPrior = 0.5*elem_div(pow ( Bmsy - mBmsy, 2 ), pow(sBmsy,2));  
  }
  
  // Then Umsy
  if (phz_Umsy > 0)
  {
    UmsyPrior = 0.5*elem_div( pow ( Umsy - mUmsy, 2 ), pow(sUmsy,2));  
  }
  
  // lnq prior (shared)
  if (phz_lnq > 0)
  {
    if (nS == 1) lnqPrior = pow ( lnq_s - mlnq, 2 ) / s2lnq /2;
    else {
      lnqPrior = 0.5*log(mfexp(lntau2q)) + 0.5*norm2(lnq_s - lnqbar)/mfexp(lntau2q);
      lnqPrior += lntau2q;
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
  if (active(lntau2))
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

FUNCTION mcDumpOut
  if( mcHeader == 0 )
  {
    // Output the parameters needed for performance testing of SA method.
    mcoutpar << "Bmsy msy Umsy q Sigma2 kappa2 tau2 UnT_bar BnT dep_bar " << endl;
  
    for ( int s=1;s<=nS;s++)
    {
      // Output the parameter values for this replicate
      mcoutpar << Bmsy(s) << " " << msy(s) <<" "<< Umsy(s) <<" "<< mfexp(lnq_s(s));
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
      mcoutpar << Bmsy(s) << " " << msy(s) <<" "<< Umsy(s) <<" "<< mfexp(lnq_s(s));
      mcoutpar <<" "<< Sigma2 <<" "<< kappa2 <<" " << tau2 <<" " << UnT_bar(s) << " ";
      mcoutpar << Bt(s,nT) <<" "<< dep_bar(s) << endl;

      // Output the biomass estimates for this MC evaluation
      mcoutbio << Bt(s) << endl;
    }
  }

REPORT_SECTION
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
  report << mfexp(lnq_s) <<endl;
  if (nS > 1)
  {
    report << "# qbar" << endl;
    report << mfexp(lnqbar) <<endl;  
    report << "tau2q" << endl;
    report << mfexp(lntau2q) << endl;
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


TOP_OF_MAIN_SECTION
  arrmblsize = 50000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
  //gradient_structure::set_MAX_NVAR_OFFSET(5000);
  //gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);

GLOBALS_SECTION
  #include <admodel.h>
  
  int mcHeader=0;
  int mcTrials=0;
  ofstream mcoutpar("msProdParMC.dat");
  ofstream mcoutbio("msProdBioMC.dat");



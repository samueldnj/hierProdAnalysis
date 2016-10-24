// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// msProdCV.tpl
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
  init_int phz_Bmsy;
  init_int phz_Fmsy;
  init_int phz_tau;             // tau - obs error sd
  init_int phz_sigma;           // sigma - proc error sd
  init_int phz_q;               // q - survey catchability
  init_int phz_mlnq;
  init_int phz_slnq;
  init_int phz_dev;
  init_int phz_AR;
  init_int phz_chol;
 
  // dummy variable to check data read
  init_int dumm;                      // dummy variable
  
  // Constant for number of free pars in chol matx
  int cEntries;

  // Procedure to exit if not all data is read correctly
  LOC_CALCS
    if(dumm!=999)
    {
      cout<<"Error reading data.\n Fix it!! \n dumm = " << dumm << endl;
      ad_exit(1);
    }
    cEntries=nS*(nS-1)/2;

PARAMETER_SECTION
  // parameters to estimate (mostly on log scale - found in pin file)
  init_bounded_vector lnBmsy(1,nS,2,10,phz_Bmsy);          // Bmsy - log scale
  init_bounded_vector lnFmsy(1,nS,-5,0,phz_Fmsy);          // Umsy - log scale
  init_bounded_vector lnTau2(1,nS,-5,1,phz_tau);            // tau - osb error sd
  init_bounded_vector lnSigma2(1,nS,-5,1,phz_sigma);        // Sigma - MS proc error sd
  init_bounded_number lnsigma2(-5,1,phz_sigma);             // sigma - env proc error sd
  //init_bounded_vector lnq(1,nS,-5,1,phz_q);               // survey catchability

  // Prior hyperparameters
  init_vector mBmsy(1,nS,-1);       // msy prior mean (species spec)
  init_vector sBmsy(1,nS,-1);       // msy prior SD (species spec)
  init_vector mFmsy(1,nS,-1);       // lnFmsy prior mean (species spec)
  init_vector sFmsy(1,nS,-1);       // lnFmsy priod sd (species spec)
  init_number alpha_sigmaMS(-1);    // sigma prior shape (shared)
  init_number beta_sigmaMS(-1);     // sigma prior scale (shared)
  init_vector alpha_Sigma(1,nS,-1);  // sigma prior shape (shared)
  init_vector beta_Sigma(1,nS,-1);   // sigma prior scale (shared)
  init_vector alpha_tau(1,nS,-1);    // tau prior shape (shared)
  init_vector beta_tau(1,nS,-1);     // tau prior scale (shared)
  init_number mlnq(phz_mlnq);       // logq prior mean (shared)
  init_number slnq(phz_slnq);       // logq prior sd (shared)

  // Covariance parameters
  init_bounded_number rho(0,1,phz_AR);              // auto-regression factor
  init_bounded_vector c(1,cEntries,-1,1,phz_chol);
  matrix chol(1,nS,1,nS);                           // Matrix to hold cholesky factor

  // process error deviations
  init_bounded_vector epst(1,nT,-5.,5.,phz_dev);  // environmental proc error
  init_bounded_matrix zetat(1,nS,1,nT,-5.,5.,phz_dev); // species-specific proc error

  //objective function value
  objective_function_value f;

  // back-transformed parameters
  vector Bmsy(1,nS);        
  vector Fmsy(1,nS);
  vector tau(1,nS);
  vector Sigma(1,nS);
  number sigma;
  
  // variance parameters
  vector tau2(1,nS);
  vector Sigma2(1,nS);
  number sigma2;

  // variables to hold concentrated parameters
  sdreport_vector lnqhat(1,nS);
  sdreport_vector q(1,nS);

  //penalizer
  number fpen;               // penalty to obj function for bad dynamics
  number pospen;             // penalty to add to fpen if biomass dips
  
  // vectors to hold predicted biomass and index values
  matrix Bt_bar(1,nS,1,nT);        // predicted biomass
  matrix It_bar(1,nS,1,nT);        // predicted IoA

  // variables to hold derived values
  vector msy(1,nS);           // Biomass at msy
  vector UnT_bar(1,nS);       // estimated fishing mortality 
  vector dep_bar(1,nS);       // depletion estimate
  vector epstCorr(1,nT);      // Autocorrelated epst values
  matrix zetatCorr(1,nS,1,nT);// Correlated zetat values

  // likelihood quantities
  number envNLL;            // NLL for environmental proc error
  vector procNLL(1,nS);     // NLL for ms interaction error
  vector obsNLL(1,nS);      // NLL for observation error

  vector nll(1,nS);         // total nll (sum of above)

  // Prior quantities
  vector prior(1,nS);         // total prior log density
  
  // Variables to hold quantities used for concentration of parameters
  vector SSRobs(1,nS);    // sum of squared resids of obs error


PROCEDURE_SECTION
  // Initialise obj function var
  f = 0.;

  // initialise derived variables
  lnqhat=0.;
  q = 0.;
  chol = 0.;
  dep_bar = 0.;
  UnT_bar = 0.;

  // Run state dynamics function
  stateDynamics();
  // cout << "State Dynamic complete" << endl;
  // apply observation model, compute conditional value for lnq
  obsModel();
  // cout << "Observation model complete" << endl;
  // Compute likelihood for dynamics
  calcLikelihoods();
  // cout << "Likelihood computation complete" << endl;
  // cout << "nll = " << nll << endl;
  // Compute priors
  calcPriors();
  // cout << "Prior computation complete" << endl;
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
  f = sum(nll) + envNLL;  // Likelihood and proc error pen
  f += sum(prior);        // sum species spec priors
  f += fpen;              // positive function value penalties

// State dynamics subroutine
FUNCTION stateDynamics
  // exponentiate leading parameters
  Fmsy = mfexp ( lnFmsy );         // optimal fishing mortality
  Bmsy = mfexp ( lnBmsy );           // msy
  msy = elem_prod( Bmsy, Fmsy);     // optimal biomass

  // initialise variables
  Bt_bar = 0.; chol = 0.; epstCorr = 0.; zetatCorr = 0.;

  // Create cholesky factor of correlation mtx
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

  // Compute auto-correlated epst values
  epstCorr(1) = epst(1);
  for(int t=1;t<nT;t++)
  {
    epstCorr(t+1) = rho*epstCorr(t) + epst(t+1);
  }

  // Compute correlated zetat values
  zetatCorr = chol * zetat;
  
  // reinitialise penalisers
  fpen = 0.; pospen = 0.;
  
  // Run a loop for population dynamics of each stock
  // Assume stocks are at equilibrium initially ( B1 = B0e^delta1 )
  for (int s=1;s<=nS;s++)
  {
    Bt_bar(s,1) = ( 2. * Bmsy(s) ) * mfexp(epstCorr(1)+zetatCorr(s,1));
    // loop over time periods to build dynamic model for t > 1
    for(int t=1; t<nT; t++)
    {
      pospen = 0.;
      // Run state dynamics equation, add process error
      Bt_bar(s,t+1) = ( Bt_bar(s,t) + 
                      2. * Fmsy(s) * Bt_bar(s,t) * (1 - Bt_bar(s,t)/Bmsy(s)/2.0 ) - 
                      katch(s,t) ) * mfexp ( epstCorr(t+1)+zetatCorr(s,t+1) );
      Bt_bar(s,t+1) = posfun ( Bt_bar(s,t+1), 10e-1, pospen );
      
      // Increment function penaliser variable
      fpen += 1000. * pospen;
    }
  }
  
  // compute derived management quantities
  dep_bar = elem_div(column(Bt_bar, nT ) ,Bmsy) / 2;      // depletion estimate
  UnT_bar = elem_div(column(katch,nT),column(Bt_bar,nT)); // Most recent exp rate
  
// function to compute predicted observations and residuals
FUNCTION obsModel
  for ( int s=1;s<=nS;s++)
  {
    // Compute conditional estimate of lnq_hat
    lnqhat(s) = sum ( log ( It(s) ) - log ( Bt_bar(s) ) ) / nT ;
    q(s) = mfexp ( lnqhat(s));
    // Compute predicted index of abundance
    It_bar(s) = q(s) * Bt_bar(s);

    // Compute the sum of squared residuals
    SSRobs(s) = norm2 ( log ( It(s) ) - log ( It_bar(s) ) );
  }
// Subroutine to compute negative log likelihoods
FUNCTION calcLikelihoods
  // Initialise NLL variables
  obsNLL = 0.; procNLL = 0.; nll = 0.; envNLL=0.;
  
  // back-calculate leading pars
  sigma2   = mfexp(lnsigma2);
  Sigma2   = mfexp(lnSigma2); 
  tau2     = mfexp(lnTau2);

  // Compute variances
  sigma   = sqrt(sigma2);
  Sigma   = sqrt(Sigma2);
  tau     = sqrt(tau2);
  
  // compute observation model likelihood
  obsNLL = 0.5*nT*log(tau2) + 0.5*elem_div(SSRobs,tau2);
  // Compute proc error conditional variance, and proc error likelihood
  for ( int s=1;s<=nS;s++)
  {
    procNLL(s) = nT*0.5*log(Sigma2(s)) + 0.5*norm2(zetat(s))/Sigma2(s);
  }
  // Compute environmental error likelihood - this is added at line #157 ish
  envNLL = 0.5*nT*log(sigma2) + 0.2*norm2(epst)/sigma2;

  // sum likelihoods
  nll = obsNLL + procNLL;

// Subroutine to compute prior densities to add to posterior
FUNCTION calcPriors
  // Initialise prior var
  prior = 0.;
  // First, msy
  prior = elem_div(pow ( Bmsy - mBmsy, 2 ), pow(sBmsy,2)) / 2.;
  // Then Fmsy
  prior += elem_div(pow ( Fmsy - mFmsy, 2 ), pow(sFmsy,2)) / 2.;
  // Now Sigma2hat prior
  prior += elem_prod((alpha_Sigma+1),log(Sigma2))+elem_div(beta_Sigma,Sigma2);
  // sigma2hat prior
  prior += (alpha_sigmaMS+1)*log(sigma2) + beta_sigmaMS/sigma2;
  // tau2hat prior
  prior += elem_prod((alpha_tau+1),log(tau2))+elem_div(beta_tau,tau2);
  // lnq prior
  prior += pow( lnqhat - mlnq, 2 ) / slnq / slnq / 2;


FUNCTION mcDumpOut
  if( mcHeader == 0 )
  {
    // Output the parameters needed for performance testing of SA method.
    mcoutpar << "Bmsy msy Fmsy q sigma2 tau2 UnT_bar BnT dep_bar Sigma2" << endl;
  
    for ( int s=1;s<=nS;s++)
    {
      // Output the parameter values for this replicate
      mcoutpar << Bmsy(s) << " " << msy(s) <<" "<< Fmsy(s) <<" "<< q(s);
      mcoutpar <<" "<< sigma2 <<" "<< tau2(s) <<" " << UnT_bar(s) <<" ";
      mcoutpar << Bt_bar(s,nT) <<" "<< dep_bar(s) << " " << Sigma2(s) << endl;

      // Output the biomass estimates for this MC evaluation
      mcoutbio << Bt_bar(s) << endl;
    }

    mcHeader = 1;
  }
  else
  {
    // loop over species
    for ( int s=1;s<=nS;s++)
    {
      // Output the parameter values for this replicate
      mcoutpar << Bmsy(s) << " " << msy(s) <<" "<< Fmsy(s) <<" "<< q(s);
      mcoutpar <<" "<< sigma2 <<" "<< tau2(s) <<" " << UnT_bar(s) <<" ";
      mcoutpar << Bt_bar(s,nT) <<" "<< dep_bar(s) << " " <<  Sigma2(s) << endl;

      // Output the biomass estimates for this MC evaluation
      mcoutbio << Bt_bar(s) << endl;
    }
  }

REPORT_SECTION
  // Output estimates and derived variables for comparison to true
  // values in sim-est experiments
  report << "## Single Species Production Model Results" << endl;
  report << "## Parameter estimates " << endl;
  report << "# Bmsy" << endl;
  report << Bmsy << endl;
  report << "# Fmsy" << endl;
  report << Fmsy << endl;
  report << "# sigma" << endl;
  report << sigma << endl;
  report << "# Sigma" << endl;
  report << Sigma << endl;
  report << "# tau" << endl;
  report << tau << endl;
  report << "# epst" << endl;
  report << epst << endl;
  report << "# zetat" << endl;
  report << zetat << endl;
  report << "# c" << endl;
  report << c << endl;
  report << "# rho" << endl;
  report << rho << endl;
  report << endl;
  
  report << "## Derived variables" << endl;
  report << "# msy" << endl;
  report << msy << endl;
  report <<"# q" << endl;
  report << q <<endl;
  report << "# Sigma2" << endl;
  report << Sigma2 << endl;
  report << "# sigma2" << endl;
  report << sigma2 << endl;
  report << "# tau2" << endl;
  report << tau2 << endl;
  report << "# D" << endl;
  report << dep_bar << endl;
  report << "# lastCt" << endl;
  report << column(katch,nT) << endl;
  report << "# Bt" << endl;
  report << Bt_bar << endl;
  report << "# Ut" << endl;
  report << elem_div(katch,Bt_bar) << endl;
  report << "# epstCorr" << endl;
  report << epstCorr << endl;
  report << "# zetatCorr" << endl;
  report << zetatCorr << endl;
  report << "# chol" << endl;
  report << chol << endl;
  report << endl;


  report << "## Priors" << endl;
  report << "# mBmsy" << endl;  
  report << mBmsy << endl;
  report << "# sBmsy" << endl;  
  report << sBmsy << endl;
  report << "# mFmsy" << endl;  
  report << mFmsy << endl;
  report << "# sdFmsy" << endl;  
  report << sFmsy << endl;
  report << "# alpha_sigma" << endl;
  report << alpha_sigmaMS << endl;
  report << "# beta_sigma" << endl;
  report << beta_sigmaMS << endl;
  report << "# alpha_Sigma" << endl;
  report << alpha_Sigma << endl;
  report << "# beta_Sigma" << endl;
  report << beta_Sigma << endl;
  report << "# alpha_tau" << endl;
  report << alpha_tau << endl;
  report << "# beta_tau" << endl;
  report << beta_tau << endl;
  report << "# mlnq" << endl;  
  report << mlnq << endl;
  report << "# slnq" << endl;  
  report << slnq << endl;
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
  report << nll << endl;
  report << "# Env_err_lik" << endl;
  report << envNLL << endl;
  report << "# total_priors" << endl;
  report << prior << endl;
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
  ofstream mcoutpar("msProdCVParMC.dat");
  ofstream mcoutbio("msProdCVBioMC.dat");



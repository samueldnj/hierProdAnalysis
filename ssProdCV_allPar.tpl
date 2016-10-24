// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// ssProd.tpl
// 
// Single species surplus production model with logistic stock-recruit function
// 
// Author: Samuel Johnson
// Date: August 25, 2014                    
// Purpose: To fit a single species production model to commercial catch
// and survey CPUE data.
// 
// ToDo:  1. logit transformations for bounding pars?
//        2. Hope the motherfucker fits.
//        3. Input estimation phases for each par group
//       
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

DATA_SECTION
  // Control parameters for data
  init_int nT;                        // number of time iterations
  init_vector katch(1,nT);            // catch in kg
  init_vector It(1,nT);               // indices of abundance

  // Estimation phases
  init_int phz_Bmsy;            // Bmsy - scale
  init_int phz_Fmsy;            // Fmsy - rate
  init_int phz_tau;             // tau - obs error sd
  init_int phz_sigma;           // sigma - proc error sd
  init_int phz_q;               // q - survey catchability
  init_int phz_mlnq;            // q prior mean
  init_int phz_slnq;            // q prior sd
  init_int phz_dev;             // proc error deviations
  init_int phz_AR;              // AR(1) factor
 
  // dummy variable to check data read
  init_int dumm;                      // dummy variable
  
  // Procedure to exit if not all data is read correctly
  LOC_CALCS
    if(dumm!=999)
    {
      cout<<"Error reading data.\n Fix it!! \n dumm = " << dumm << endl;
      ad_exit(1);
    }

PARAMETER_SECTION
  //parameters to estimate (mostly on log scale - found in pin file)
  init_bounded_number lnBmsy(2,10,phz_Bmsy);         // carrying capacity - log scale
  init_bounded_number lnFmsy(-5,0,phz_Fmsy);         // intrinsic rate of growth - log scale
  init_bounded_number lnTau2(-5,1,phz_tau);           // obs error sd
  init_bounded_number lnsigma2(-5,1,phz_sigma);       // proc error sd
  init_bounded_number lnq(-5,1,phz_q);               // survey catchability

  // Fixed parameters and hyperparameters
  init_number mBmsy(-1);              // msy prior mean
  init_number sBmsy(-1);              // msy prior SD
  init_number mFmsy(-1);              // lnFmsy prior mean
  init_number sFmsy(-1);              // lnFmsy priod sd
  init_number alpha_sigma(-1);         // sigma prior mean
  init_number beta_sigma(-1);          // sigma prior sd
  init_number alpha_tau(-1);           // sigma prior mean
  init_number beta_tau(-1);            // sigma prior sd
  init_number mlnq(phz_mlnq);         // log q prior mean
  init_number slnq(phz_slnq);         // log q prior sd

  // autocorrelation factor
  init_bounded_number rho(0,0.9,phz_AR);         // AR(1) factor

  // process error deviations
  init_bounded_vector epst(1,nT,-2.,2.,phz_dev);

  //objective function value
  objective_function_value f;

  // back-transformed leading parameters
  number Bmsy;         
  number Fmsy;
  number sigma;
  number tau;
  sdreport_number q;
  
  // variances
  number sigma2;
  number sigma2corr;
  number tau2;

  // variables to hold concentrated parameters
  sdreport_number lnqhat;

  // Autocorrelated proc error deviations
  vector epstCorr(1,nT);

  //penalizer
  number fpen;               // penalty to obj function for bad dynamics
  number pospen;             // penalty to add to fpen if biomass dips
  
  // vectors to hold predicted biomass and index values
  vector Bt_bar(1,nT);        // predicted biomass
  vector It_bar(1,nT);        // predicted IoA

  // variables to hold derived values
  number UnT_bar;     // estimated fishing mortality 
  number dep_bar;     // depletion estimate

  // likelihood quantities
  number obsNLL;      // NLL for observation error
  number procNLL;     // NLL for process error
  number nll;         // total nll (sum of above)   

  // Prior quantities
  number prior;         // total prior log density
  
  // Variables to hold concentration quantities
  number SSRobs;    // sum of squared resids of obs error
  !! cout << "OK to end of Parameter Section" << endl;

PROCEDURE_SECTION
  // Initialise obj function var
  f = 0.;

  //cout << "Starting Procedure Section" << endl;

  // initialise derived variables
  lnqhat = 0.;
  dep_bar=0.; UnT_bar=0.;

  // Run state dynamics function
  stateDynamics();
  //cout << "stateDynamics ok" << endl;

  // apply observation model, compute conditional value for lnq
  obsModel();
  //cout << "obsModel ok" << endl;

  // Compute likelihood for dynamics
  calcLikelihoods();
  //cout << "calcLikelihoods ok" << endl;

  // Compute priors
  calcPriors();
  //cout << "calcPriors ok" << endl;

  // Output MCMC data if running mcmc trials
  if(mceval_phase())
  { 
    mcDumpOut();
  }

  // compute objective function value
  f = nll;      // Likelihood and proc error pen
  f += prior;   // priors
  f += fpen;    // function value penalties

// State dynamics subroutine
FUNCTION stateDynamics
  // exponentiate leading parameters
  Fmsy = mfexp ( lnFmsy );            // optimal exploitation rate
  //cout << "Fmsy" << " = " << Fmsy << endl;
  Bmsy = mfexp ( lnBmsy );            // optimal Biomass
  //cout << "Bmsy" << " = " << Bmsy << endl;
  
  // reinitialise penalisers
  fpen = 0.; pospen = 0.;

  // Start populating AR(1) proc error deviations
  epstCorr(1) = epst(1);
  
  // Run a loop for population dynamics of each stock
  // Assume stocks are at equilibrium initially ( B1 = B0e^delta1 )
  Bt_bar(1) = ( 2. * Bmsy ) * mfexp(epstCorr(1));
  //cout << "B1" << " = " << Bt_bar(1) << endl;
  
  // loop over time periods to build dynamic model for t > 1
  for(int t=1; t<nT; t++)
  { 
    //cout << "epstCorr" << t << " = " << epstCorr(t) << endl;
    epstCorr(t+1) = rho * epstCorr(t) + epst(t+1);

    pospen = 0.;
    // Run state dynamics equation, add process error
    Bt_bar(t+1) = ( Bt_bar(t) + 
                    2. * Fmsy * Bt_bar(t) * (1 - Bt_bar(t)/Bmsy/2. ) - 
                    katch(t) ) * mfexp ( epstCorr(t+1) );
    Bt_bar(t+1) = posfun ( Bt_bar(t+1), 1., pospen );
    //cout << "Computing B" << t+1 << " = " << Bt_bar(t+1) << endl;
    // Increment function penaliser variable
    fpen += 100. * pospen;
  }

  // compute derived performance values
  UnT_bar = katch ( nT ) / Bt_bar(nT);  // comparison of F
  dep_bar = Bt_bar ( nT ) / Bmsy / 2;     // depletion estimate
  
// function to compute predicted observations and residuals
FUNCTION obsModel
  // Compute conditional estimate of lnq_hat
  //lnqhat = sum ( log ( It ) - log ( Bt_bar ) ) / nT ;
  q = mfexp(lnq);
  //cout << "q = " << q << endl;
  // Compute predicted index of abundance
  It_bar = q * Bt_bar;
  //cout << "It = " << It_bar << endl;

  // Compute the sum of squares
  SSRobs = norm2 ( log ( It ) - log ( It_bar ) );
  //cout << "SSRobs = " << SSRobs << endl;

// Subroutine to compute negative log likelihoods
FUNCTION calcLikelihoods
  // Initialise NLL variables
  obsNLL = 0.; procNLL = 0.; nll = 0.;
  sigma = 0.; tau = 0.; sigma2 = 0.; tau2 = 0.;

  // Back-transform sigma and tau
  sigma2    = mfexp(lnsigma2);
  tau2      = mfexp(lnTau2);
  sigma2corr= sigma2/(1-rho*rho);

  // compute variances
  sigma       = sqrt(sigma2);
  tau         = sqrt(tau2);

  // compute observation model likelihood
  obsNLL = 0.5*nT*log(tau2) + 0.5*SSRobs/tau2;
  // Compute proc error conditional variance
  procNLL = 0.5*nT*log(sigma2corr) + 0.5*norm2(epstCorr)/sigma2corr;

  // sum likelihoods
  nll = obsNLL + procNLL;

// Subroutine to compute prior densities to add to posterior
FUNCTION calcPriors
  // Initialise prior var
  prior = 0.;
  // First, msy
  prior = ( Bmsy - mBmsy)*(Bmsy - mBmsy) / sBmsy / sBmsy / 2.;
  // Then Fmsy
  prior += (Fmsy - mFmsy)*(Fmsy - mFmsy) / sFmsy / sFmsy / 2.;
  // Now sigma2hat prior
  prior += (alpha_sigma + 1) * log(sigma2) + beta_sigma / sigma2;
  // tau2hat prior
  prior += (alpha_tau + 1) * log(tau2) + beta_tau / tau2;
  // lnq prior
  prior += (lnq - mlnq)*(lnq-mlnq) / slnq / slnq / 2;


FUNCTION mcDumpOut
  if( mcHeader == 0 )
  {
    // Output the parameters needed for performance testing of SA method.
    mcoutpar << "Bmsy msy Fmsy q sigma2 tau2 UnT_bar BnT dep_bar" << endl;
  
    // Output the parameter values for this replicate
    mcoutpar << Bmsy << " " << Bmsy*Fmsy <<" "<< Fmsy <<" "<< q;
    mcoutpar <<" "<< sigma2 <<" "<< tau2 <<" " << UnT_bar <<" ";
    mcoutpar << Bt_bar(nT) <<" "<< dep_bar << endl;

    // Output the biomass estimates for this MC evaluation
    mcoutbio << Bt_bar << endl;

    mcHeader = 1;
  }
  else
  {
    // Output the parameter values for this replicate
    mcoutpar << Bmsy << " " << Bmsy*Fmsy <<" "<< Fmsy <<" "<< q; 
    mcoutpar <<" "<< sigma2 <<" "<< tau2 <<" " << UnT_bar <<" ";
    mcoutpar << Bt_bar(nT) <<" "<< dep_bar << endl;

    // Output the biomass estimates for this MC evaluation
    mcoutbio << Bt_bar << endl;
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
  report << "# tau" << endl;
  report << tau << endl;
  report << "# sigma" << endl;
  report << sigma << endl;
  report << "# epst" << endl;
  report << epst << endl;
  report << "# rho" << endl;
  report << rho << endl;
  report << endl;
  
  report << "## Derived variables" << endl;
  report << "# msy" << endl;
  report << Bmsy * Fmsy << endl;
  report <<"# q" << endl;
  report << q <<endl;
  report << "# sigma2" << endl;
  report << sigma*sigma << endl;
  report << "# tau2" << endl;
  report << tau*tau << endl;
  report << "# D" << endl;
  report << dep_bar << endl;
  report << "# lastCt" << endl;
  report << katch(nT) << endl;
  report << "# Bt" << endl;
  report << Bt_bar << endl;
  report << "# Ut" << endl;
  report << elem_div(katch,Bt_bar) << endl;
  report << "# epstCorr" << endl;
  report << epstCorr << endl;
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
  report << alpha_sigma << endl;
  report << "# beta_sigma" << endl;
  report << beta_sigma << endl;
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
  report << "# Ct" << endl;
  report << katch << endl;
  report << "# It" << endl;
  report << It << endl;
  report << endl;
  
  report << "## Minimization properties" << endl;
  report << "# total_likelihood" << endl;
  report << nll << endl;
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
  ofstream mcoutpar("ssProdCVParMC.dat");
  ofstream mcoutbio("ssProdCVBioMC.dat");



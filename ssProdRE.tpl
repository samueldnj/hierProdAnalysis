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
//       
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

DATA_SECTION
  // Control parameters for data
  init_int nT;                        // number of time iterations
  init_vector katch(1,nT);            // catch in kg
  init_vector It(1,nT);               // indices of abundance
 
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
  // Leading pop dynamics pars (mostly on log scale - found in pin file)
  init_bounded_number lnBMSY(0,20,1);            // carrying capacity - log scale
  init_bounded_number lnFMSY(-5,0,-2);            // intrinsic rate of growth - log scale
  init_number lnq(-3);              // catchability

  // Hyperparameters for prior distributions
  init_number mBMSY(-1);           // MSY prior mean
  init_number sBMSY(-1);           // MSY prior SD
  init_number mFMSY(-1);          // lnFMSY prior mean
  init_number sFMSY(-1);          // lnFMSY priod sd
  init_number alphaSigma(-1);     // sigma prior mean
  init_number betaSigma(-1);      // sigma prior sd
  init_number alphaTau(-1);       // sigma prior mean
  init_number betaTau(-1);        // sigma prior sd
  init_number mlnq(-1);           // log q prior mean
  init_number slnq(-1);           // log q prior sd

  // error stuff
  init_number lntau2(-2);
  init_number lnsigma2(-2);                          // process error sd log scale
  init_bounded_number rho(0,0.99,-3);                // AR(1) factor
  init_bounded_vector u(1,nT,-3.,3.,-2);              // standard normal deviations


  //objective function value
  objective_function_value f;

  // back-transformed parameters
  number q;
  number sigma2; number sigma;
  number tau2; number tau;

  // vector for proc error deviations
  vector epst(1,nT);        // uncorrelated
  vector epstCorr(1,nT);    // AR(1)
  
  // variables to hold concentrated parameters
  sdreport_number lnqhat;
  number tau2hat;

  //penalizer
  number fpen;               // penalty to obj function for bad dynamics
  number pospen;             // penalty to add to fpen if biomass dips
  
  // vectors to hold predicted biomass and index values
  vector Bt_bar(1,nT);        // predicted biomass
  vector It_bar(1,nT);        // predicted IoA

  // variables to hold derived values
  number BMSY;        // Biomass at MSY
  number FMSY;
  number FnT_bar;     // estimated fishing mortality 
  number dep_bar;     // depletion estimate

  // likelihood quantities
  number obsNLL;      // NLL for observation error
  number procNLL;     // NLL for process error
  number nll;         // total nll (sum of above)   

  // Prior quantities
  number prior;         // total prior log density
  
  // Variables to hold concentration quantities
  number SSRobs;    // sum of squared resids of obs error


PROCEDURE_SECTION
  // Initialise obj function var
  f = 0.;

  // initialise derived variables
  //lnqhat = 0.; tau2hat = 0.;
  q = mfexp(lnq);
  tau2 = mfexp(lntau2); tau=sqrt(tau2);
  sigma2 = mfexp(lnsigma2); sigma=sqrt(sigma2);

  // Run state dynamics function
  stateDynamics();
  cout << "Pop Dynamics done" << endl;
  // apply observation model, compute conditional value for lnq
  obsModel();
  cout << "obs model done" << endl;
  // Compute likelihood for dynamics
  calcLikelihoods();
  cout << "Like done" << endl;
  // Compute priors
  calcPriors();
  cout << "Priors done" << endl;

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
  cout << "lnFMSY = " << lnFMSY << endl;
  FMSY = mfexp ( lnFMSY );     // optimal fishing mortality
  cout << "lnBMSY = " << lnBMSY << endl;
  BMSY = mfexp ( lnBMSY );            // MSY
  

  // reinitialise penalisers
  fpen = 0.; pospen = 0.;

  // start filling in epstCorr
  epst = sigma*u;
  cout << "u = " << u << endl;
  cout << "sigma = " << sigma << endl;
  cout << "epst = " << epst << endl;
  epstCorr(1) = epst(1);
  
  // Run a loop for population dynamics of each stock
  // Assume stocks are at equilibrium initially ( B1 = B0e^delta1 )
  Bt_bar(1) = ( 2. * BMSY ) * mfexp(epstCorr(1));
  
  // loop over time periods to build dynamic model for t > 1
  for(int t=1; t<nT; t++)
  {
    pospen = 0.;
    epstCorr(t+1) = rho*epstCorr(t) + epst(t+1);
    // Run state dynamics equation, add process error
    Bt_bar(t+1) = ( Bt_bar(t) + 
                    2. * FMSY * Bt_bar(t) * (1 - Bt_bar(t)/BMSY/2.0 ) - 
                    katch(t) ) * mfexp ( epstCorr(t+1) );
    Bt_bar(t+1) = posfun ( Bt_bar(t+1), 10e-1, pospen );
    
    // Increment function penaliser variable
    fpen += 1000. * pospen;
  }

  cout << "Bt_bar = " << Bt_bar << endl;

  // compute derived performance values
  //FnT_bar = katch ( nT - 1 ) / FMSY;   // comparison of F
  //dep_bar = Bt_bar ( nT ) / BMSY / 2;     // depletion estimate
  
// function to compute predicted observations and residuals
FUNCTION obsModel

  // Compute conditional estimate of lnq_hat
  //lnqhat = sum ( log ( It ) - log ( Bt_bar ) ) / nT ;

  // Compute predicted index of abundance
  It_bar = q * Bt_bar;
  //cout << "lnqhat = " << lnqhat << endl;

  // Compute the sum of squares
  SSRobs = norm2 ( log ( It ) - log ( It_bar ) );

// Subroutine to compute negative log likelihoods
FUNCTION calcLikelihoods
  // Initialise NLL variables
  obsNLL = 0.; procNLL = 0.; nll = 0.;
  // Compute observation error conditional variance
  //tau2hat = SSRobs/nT;
  //cout << "tau2hat = " << tau2hat << endl;
  // compute observation model likelihood
  obsNLL = 0.5*nT*log(tau2) + 0.5*SSRobs/tau2;
  // Compute proc error conditional variance
  procNLL =  0.5*norm2(epst)/sigma2 + 0.5*nT*log(sigma2);

  // sum likelihoods
  nll = obsNLL + procNLL;

// Subroutine to compute prior densities to add to posterior
FUNCTION calcPriors
  // Initialise prior var
  prior = 0.;
  // First, MSY
  prior = pow ( mfexp(lnBMSY) - mBMSY, 2 ) / sBMSY / sBMSY / 2.;
  // Then FMSY
  prior += pow ( mfexp(lnFMSY) - mFMSY, 2 ) / sFMSY / sFMSY / 2.;
  // Now sigma2 prior
  prior += (alphaSigma + 1) * log(sigma*sigma) + betaSigma/sigma/sigma;
  // tau2hat prior
  prior += (alphaTau + 1) * log(tau2hat) + betaTau / tau2hat;
  // lnq prior
  prior += pow( lnqhat - mlnq, 2 ) / slnq / slnq / 2;


FUNCTION mcDumpOut
  if( mcHeader == 0 )
  {
    // Output the parameters needed for performance testing of SA method.
    mcoutpar << "BMSY FMSY q sigma2 tau2 FnT_bar BnT dep_bar" << endl;
  
    // Output the parameter values for this replicate
    mcoutpar << mfexp(lnBMSY) << " " << mfexp(lnFMSY) <<" "<< mfexp(lnqhat);
    mcoutpar <<" "<< sigma*sigma <<" "<< tau2hat <<" " << FnT_bar <<" ";
    mcoutpar << Bt_bar(nT) <<" "<< dep_bar << endl;

    // Output the biomass estimates for this MC evaluation
    mcoutbio << Bt_bar << endl;

    mcHeader = 1;
  }
  else
  {
    // Condition on "good" starting values
    if( value(pospen)==0 )
    {
      // Output the parameter values for this replicate
      mcoutpar << mfexp(lnBMSY) <<" "<< mfexp(lnFMSY) <<" "<< mfexp(lnqhat); 
      mcoutpar <<" "<< sigma*sigma <<" "<< tau2hat <<" " << FnT_bar <<" ";
      mcoutpar << Bt_bar(nT) <<" "<< dep_bar << endl;

      // Output the biomass estimates for this MC evaluation
      mcoutbio << Bt_bar << endl;
    }
  }

REPORT_SECTION
  // Output estimates and derived variables for comparison to true
  // values in sim-est experiments
  report << "## Single Species Production Model Results" << endl;
  report << "## Parameter estimates " << endl;
  report << "# BMSY" << endl;
  report << mfexp(lnBMSY) << endl;
  report << "# FMSY" << endl;
  report << mfexp(lnFMSY) << endl;
  report << "# u" << endl;
  report << u << endl;
  report << "# sigma2" << endl;
  report << sigma*sigma << endl;
  report << "# rho" << endl;
  report << rho << endl;
  report << endl;
  
  report << "## Derived variables" << endl;
  report << "# epst" << endl;
  report << epst << endl;
  report << "# epstCorr" << endl;
  report << epstCorr << endl;
  report <<"# q" << endl;
  report << mfexp(lnqhat) <<endl;
  report << "# tau2" << endl;
  report << tau2hat << endl;
  report << "# D" << endl;
  report << dep_bar << endl;
  report << "# lastCt" << endl;
  report << katch(nT) << endl;
  report << "# Bt" << endl;
  report << Bt_bar << endl;
  report << "# Ut" << endl;
  report << elem_div(katch,Bt_bar) << endl;
  report << "# Ft" << endl;
  report << -1. * log(1 - elem_div(katch,Bt_bar) ) << endl;
  report << endl;

  report << "## Priors" << endl;
  report << "# mBMSY" << endl;  
  report << mBMSY << endl;
  report << "# sBMSY" << endl;  
  report << sBMSY << endl;
  report << "# mFMSY" << endl;  
  report << mFMSY << endl;
  report << "# sdFMSY" << endl;  
  report << sFMSY << endl;
  report << "# alphaSigma" << endl;
  report << alphaSigma << endl;
  report << "# betaSigma" << endl;
  report << betaSigma << endl;
  report << "# alphaTau" << endl;
  report << alphaTau << endl;
  report << "# betaTau" << endl;
  report << betaTau << endl;
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
  ofstream mcoutpar("ssProdParMC.dat");
  ofstream mcoutbio("ssProdBioMC.dat");



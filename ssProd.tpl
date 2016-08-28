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
  //parameters to estimate (mostly on log scale - found in pin file)
  init_number lnMSY(1);           // carrying capacity - log scale
  init_number lnFMSY(2);          // intrinsic rate of growth - log scale
  
  // process error deviations
  init_bounded_dev_vector epst(1,nT,-5.,5.,3);
  
  // Fixed parameters and hyperparameters
  init_number mMSY(-1);           // MSY prior mean
  init_number sMSY(-1);           // MSY prior SD
  init_number mFMSY(-1);          // lnFMSY prior mean
  init_number sFMSY(-1);          // lnFMSY priod sd
  init_number alphaSigma(-1);     // sigma prior mean
  init_number betaSigma(-1);      // sigma prior sd
  init_number alphaTau(-1);       // sigma prior mean
  init_number betaTau(-1);        // sigma prior sd
  init_number mlnq(-1);           // log q prior mean
  init_number slnq(-1);           // log q prior sd

  //objective function value
  objective_function_value f;

  // back-transformed parameters
  number MSY;         
  number FMSY;
  
  // variables to hold concentrated parameters
  sdreport_number lnqhat;
  number sigma2hat;
  number tau2hat;

  //penalizer
  number fpen;               // penalty to obj function for bad dynamics
  number pospen;             // penalty to add to fpen if biomass dips
  
  // vectors to hold predicted biomass and index values
  vector Bt_bar(1,nT);        // predicted biomass
  vector It_bar(1,nT);        // predicted IoA

  // variables to hold derived values
  number BMSY;        // Biomass at MSY
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
  lnqhat = 0.; sigma2hat = 0.; tau2hat = 0.;

  // Run state dynamics function
  stateDynamics();
  // apply observation model, compute conditional value for lnq
  obsModel();
  // Compute likelihood for dynamics
  calcLikelihoods();
  // Compute priors
  calcPriors();

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
  FMSY = mfexp ( lnFMSY );         // optimal fishing mortality
  MSY = mfexp ( lnMSY );           // MSY
  BMSY = MSY / FMSY;               // optimal biomass

  // reinitialise penalisers
  fpen = 0.; pospen = 0.;
  
  // Run a loop for population dynamics of each stock
  // Assume stocks are at equilibrium initially ( B1 = B0e^delta1 )
  Bt_bar(1) = ( 2. * BMSY ) * mfexp(epst(1));
  
  // loop over time periods to build dynamic model for t > 1
  for(int t=1; t<nT; t++)
  {
    pospen = 0.;
    // Run state dynamics equation, add process error
    Bt_bar(t+1) = ( Bt_bar(t) + 
                    2. * FMSY * Bt_bar(t) * (1 - Bt_bar(t)/BMSY/2.0 ) - 
                    katch(t) ) * exp ( epst(t+1) );
    Bt_bar(t+1) = posfun ( Bt_bar(t+1), 10e-3, pospen );
    
    // Increment function penaliser variable
    fpen += 100. * pospen;
  }

  // compute derived performance values
  //FnT_bar = katch ( nT - 1 ) / FMSY;   // comparison of F
  //dep_bar = Bt_bar ( nT ) / BMSY / 2;     // depletion estimate
  
// function to compute predicted observations and residuals
FUNCTION obsModel

  // Compute conditional estimate of lnq_hat
  lnqhat = sum ( log ( It ) - log ( Bt_bar ) ) / nT ;

  // Compute predicted index of abundance
  It_bar = mfexp(lnqhat) * Bt_bar;

  // Compute the sum of squares
  SSRobs = norm2 ( log ( It ) - log ( It_bar ) );

// Subroutine to compute negative log likelihoods
FUNCTION calcLikelihoods
  // Initialise NLL variables
  obsNLL = 0.; procNLL = 0.; nll = 0.;
  // Compute observation error conditional variance
  tau2hat = SSRobs/nT;
  // compute observation model likelihood
  obsNLL = 0.5*nT*log(SSRobs);
  // Compute proc error conditional variance
  sigma2hat = norm2(epst)/nT;
  procNLL = 0.5*nT*log(norm2(epst));

  // sum likelihoods
  nll = obsNLL + procNLL;

// Subroutine to compute prior densities to add to posterior
FUNCTION calcPriors
  // Initialise prior var
  prior = 0.;
  // First, MSY
  prior = pow ( MSY - mMSY, 2 ) / sMSY / sMSY / 2.;
  // Then FMSY
  prior += pow ( FMSY - mFMSY, 2 ) / sFMSY / sFMSY / 2.;
  // Now sigma2hat prior
  prior += (alphaSigma + 1) * log(sigma2hat) + betaSigma / sigma2hat;
  // tau2hat prior
  prior += (alphaTau + 1) * log(tau2hat) + betaTau / tau2hat;
  // lnq prior
  prior += pow( lnqhat - mlnq, 2 ) / slnq / slnq / 2;


FUNCTION mcDumpOut
  if( mcHeader == 0 )
  {
    // Output the parameters needed for performance testing of SA method.
    mcoutpar << "BMSY MSY FMSY q sigma2 tau2 FnT_bar BnT dep_bar" << endl;
  
    // Output the parameter values for this replicate
    mcoutpar << BMSY << " " << MSY <<" "<< FMSY <<" "<< mfexp(lnqhat);
    mcoutpar <<" "<< sigma2hat <<" "<< tau2hat <<" " << FnT_bar <<" ";
    mcoutpar << Bt_bar(nT) <<" "<< dep_bar << endl;

    // Output the biomass estimates for this MC evaluation
    mcoutbio << Bt_bar << endl;

    mcHeader = 1;
  }
  else
  {
    // Condition on "good" starting values
    if( pospen==0 )
    {
      // Output the parameter values for this replicate
      mcoutpar << BMSY << " " << MSY <<" "<< FMSY <<" "<< mfexp(lnqhat); 
      mcoutpar <<" "<< sigma2hat <<" "<< tau2hat <<" " << FnT_bar <<" ";
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
  report << "# MSY" << endl;
  report << MSY << endl;
  report << "# FMSY" << endl;
  report << FMSY << endl;
  report << "# epst" << endl;
  report << epst << endl;
  report << endl;
  
  report << "## Derived variables" << endl;
  report << "# BMSY" << endl;
  report << BMSY << endl;
  report <<"# q" << endl;
  report << mfexp(lnqhat) <<endl;
  report << "# sigma2" << endl;
  report << sigma2hat << endl;
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
  report << "# mMSY" << endl;  
  report << mMSY << endl;
  report << "# sMSY" << endl;  
  report << sMSY << endl;
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
  ofstream mcoutpar("mcoutSSpar.dat");
  ofstream mcoutbio("mcoutSSbio.dat");



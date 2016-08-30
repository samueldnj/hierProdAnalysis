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
  init_int nS;                        // number of species
  init_matrix katch(1,nS,1,nT);       // catch in kg
  init_matrix It(1,nS,1,nT);          // indices of abundance
 
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
  init_vector lnMSY(1,nS,1);         // msy - log scale
  init_vector lnFMSY(1,nS,2);        // Umsy - log scale
  
  // process error deviations
  init_matrix epst(1,nS,1,nT,3);     // uncorrelated - cholesky factor to be added
  
  // Fixed parameters and hyperparameters
  init_vector mMSY(1,nS,-1);      // MSY prior mean (species spec)
  init_vector sMSY(1,nS,-1);      // MSY prior SD (species spec)
  init_vector mFMSY(1,nS,-1);     // lnFMSY prior mean (species spec)
  init_vector sFMSY(1,nS,-1);     // lnFMSY priod sd (species spec)
  init_number alphaSigma(-1);     // sigma prior shape (shared)
  init_number betaSigma(-1);      // sigma prior scale (shared)
  init_number alphaTau(-1);       // tau prior shape (shared)
  init_number betaTau(-1);        // tau prior scale (shared)
  init_number mlnq(-1);           // logq prior mean (shared)
  init_number slnq(-1);           // logq prior sd (shared)

  //objective function value
  objective_function_value f;

  // back-transformed parameters
  vector MSY(1,nS);        
  vector FMSY(1,nS);
  
  // variables to hold concentrated parameters
  sdreport_vector lnqhat(1,nS);
  vector sigma2hat(1,nS);
  vector tau2hat(1,nS);

  //penalizer
  number fpen;               // penalty to obj function for bad dynamics
  number pospen;             // penalty to add to fpen if biomass dips
  
  // vectors to hold predicted biomass and index values
  matrix Bt_bar(1,nS,1,nT);        // predicted biomass
  matrix It_bar(1,nS,1,nT);        // predicted IoA

  // variables to hold derived values
  vector BMSY(1,nS);        // Biomass at MSY
  vector FnT_bar(1,nS);     // estimated fishing mortality 
  vector dep_bar(1,nS);     // depletion estimate

  // likelihood quantities
  vector obsNLL(1,nS);      // NLL for observation error
  vector procNLL(1,nS);     // NLL for process error
  vector nll(1,nS);         // total nll (sum of above)   

  // Prior quantities
  vector prior(1,nS);         // total prior log density
  
  // Variables to hold quantities used for concentration of parameters
  vector SSRobs(1,nS);    // sum of squared resids of obs error


PROCEDURE_SECTION
  // Initialise obj function var
  f = 0.;

  // initialise derived variables
  lnqhat = 0.; sigma2hat = 0.; tau2hat = 0.;

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

  // cout << "lnMSY = " << lnMSY << endl;
  // cout << "lnFMSY = " << lnFMSY << endl;
  // cout << "epst = " << epst << endl;

  // Output MCMC data if running mcmc trials
  if(mceval_phase())
  { 
    mcDumpOut();
  }

  // compute objective function value
  f = sum(nll);       // Likelihood and proc error pen
  f += sum(prior);    // sum species spec priors
  f += fpen;          // positive function value penalties

// State dynamics subroutine
FUNCTION stateDynamics
  // exponentiate leading parameters
  FMSY = mfexp ( lnFMSY );         // optimal fishing mortality
  MSY = mfexp ( lnMSY );           // MSY
  BMSY = elem_div( MSY, FMSY);     // optimal biomass

  // reinitialise penalisers
  fpen = 0.; pospen = 0.;
  
  // Run a loop for population dynamics of each stock
  // Assume stocks are at equilibrium initially ( B1 = B0e^delta1 )
  for (int s=1;s<=nS;s++)
  {
    // cout << "In State Dynamics looping over species" << endl;
    // cout << "Species " << s << endl;
    Bt_bar(s,1) = ( 2. * BMSY(s) ) * mfexp(epst(s,1));
    // loop over time periods to build dynamic model for t > 1
    for(int t=1; t<nT; t++)
    {
      pospen = 0.;
      // Run state dynamics equation, add process error
      Bt_bar(s,t+1) = ( Bt_bar(s,t) + 
                      2. * FMSY(s) * Bt_bar(s,t) * (1 - Bt_bar(s,t)/BMSY(s)/2.0 ) - 
                      katch(s,t) ) * exp ( epst(s,t+1) );
      Bt_bar(s,t+1) = posfun ( Bt_bar(s,t+1), 10e-3, pospen );
      
      // Increment function penaliser variable
      fpen += 100. * pospen;
    }
    // cout << "Bt_bar = " << Bt_bar << endl;
  }
  
  
  

  // compute derived performance values
  //FnT_bar = katch ( nT - 1 ) / FMSY;   // comparison of F
  //dep_bar = Bt_bar ( nT ) / BMSY / 2;     // depletion estimate
  
// function to compute predicted observations and residuals
FUNCTION obsModel
  for ( int s=1;s<=nS;s++)
  {
    // Compute conditional estimate of lnq_hat
    lnqhat(s) = sum ( log ( It(s) ) - log ( Bt_bar(s) ) ) / nT ;

    // Compute predicted index of abundance
    It_bar(s) = mfexp(lnqhat(s)) * Bt_bar(s);

    // Compute the sum of squares
    SSRobs(s) = norm2 ( log ( It(s) ) - log ( It_bar(s) ) );
  }
// Subroutine to compute negative log likelihoods
FUNCTION calcLikelihoods
  // Initialise NLL variables
  obsNLL = 0.; procNLL = 0.; nll = 0.;
  // Compute observation error conditional variance
  tau2hat = SSRobs/nT;
  // compute observation model likelihood
  obsNLL = 0.5*nT*log(SSRobs);
  // Compute proc error conditional variance
  for ( int s=1;s<=nS;s++)
  {
    sigma2hat(s) = norm2(epst(s))/nT;
    procNLL(s) = 0.5*nT*log(norm2(epst(s)));
  
  }
  // sum likelihoods
  nll = obsNLL + procNLL;

// Subroutine to compute prior densities to add to posterior
FUNCTION calcPriors
  // Initialise prior var
  prior = 0.;
  // First, MSY
  prior = elem_div(pow ( MSY - mMSY, 2 ), pow(sMSY,2)) / 2.;
  // Then FMSY
  prior += elem_div(pow ( FMSY - mFMSY, 2 ), pow(sFMSY,2)) / 2.;
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
  
    for ( int s=1;s<=nS;s++)
    {
      // Output the parameter values for this replicate
      mcoutpar << BMSY(s) << " " << MSY(s) <<" "<< FMSY(s) <<" "<< mfexp(lnqhat(s));
      mcoutpar <<" "<< sigma2hat(s) <<" "<< tau2hat(s) <<" " << FnT_bar(s) <<" ";
      mcoutpar << Bt_bar(s,nT) <<" "<< dep_bar(s) << endl;

      // Output the biomass estimates for this MC evaluation
      mcoutbio << Bt_bar(s) << endl;
    }

    mcHeader = 1;
  }
  else
  {
    // Condition on "good" starting values
    if( pospen==0 )
    {
      for ( int s=1;s<=nS;s++)
      {
        // Output the parameter values for this replicate
        mcoutpar << BMSY(s) << " " << MSY(s) <<" "<< FMSY(s) <<" "<< mfexp(lnqhat(s));
        mcoutpar <<" "<< sigma2hat(s) <<" "<< tau2hat(s) <<" " << FnT_bar(s) <<" ";
        mcoutpar << Bt_bar(s,nT) <<" "<< dep_bar(s) << endl;

        // Output the biomass estimates for this MC evaluation
        mcoutbio << Bt_bar(s) << endl;
      }
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
  report << column(katch,nT) << endl;
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
  ofstream mcoutpar("msProdParMC.dat");
  ofstream mcoutbio("msProdBioMC.dat");



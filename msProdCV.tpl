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
  //parameters to estimate (mostly on log scale - found in pin file)
  init_vector lnBMSY(1,nS,1);         // Bmsy - log scale
  init_vector lnFMSY(1,nS,2);         // Umsy - log scale

  // Fixed parameters and hyperparameters
  init_vector mBMSY(1,nS,-1);       // MSY prior mean (species spec)
  init_vector sBMSY(1,nS,-1);       // MSY prior SD (species spec)
  init_vector mFMSY(1,nS,-1);       // lnFMSY prior mean (species spec)
  init_vector sFMSY(1,nS,-1);       // lnFMSY priod sd (species spec)
  init_number alpha_sigmaMS(-1);    // sigma prior shape (shared)
  init_number beta_sigmaMS(-1);     // sigma prior scale (shared)
  init_vector alphaSigma(1,nS,-1);  // sigma prior shape (shared)
  init_vector betaSigma(1,nS,-1);   // sigma prior scale (shared)
  init_vector alphaTau(1,nS,-1);    // tau prior shape (shared)
  init_vector betaTau(1,nS,-1);     // tau prior scale (shared)
  init_number mlnq(3);             // logq prior mean (shared)
  init_number slnq(-3);             // logq prior sd (shared)

  // process error deviations
  init_bounded_dev_vector epst(1,nT,-5.,5.,2); // auto-regressive proc error
  init_matrix zetat(1,nS,1,nT,2);       // correlated across species

  // Covariance parameters
  init_bounded_number rho(0,1,2);              // auto-regression factor
  init_bounded_vector c(1,cEntries,-1,1,2);
  //init_vector logit_c(1,cEntries,3);     // Initialise cholesky entry vector
  matrix chol(1,nS,1,nS);                 // Matrix to hold cholesky factor

  //objective function value
  objective_function_value f;

  // back-transformed parameters
  vector BMSY(1,nS);        
  vector FMSY(1,nS);
  //number rho;
  
  // variables to hold concentrated parameters
  sdreport_vector lnqhat(1,nS);
  vector Sigma2hat(1,nS);
  number sigma2hat;
  vector tau2hat(1,nS);

  //penalizer
  number fpen;               // penalty to obj function for bad dynamics
  number pospen;             // penalty to add to fpen if biomass dips
  
  // vectors to hold predicted biomass and index values
  matrix Bt_bar(1,nS,1,nT);        // predicted biomass
  matrix It_bar(1,nS,1,nT);        // predicted IoA

  // variables to hold derived values
  vector MSY(1,nS);           // Biomass at MSY
  vector FnT_bar(1,nS);       // estimated fishing mortality 
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
  lnqhat = 0.; sigma2hat = 0.; tau2hat = 0.; Sigma2hat = 0.;
  chol = 0.; 

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
  f = sum(nll) + envNLL;  // Likelihood and proc error pen
  f += sum(prior);        // sum species spec priors
  f += fpen;              // positive function value penalties

// State dynamics subroutine
FUNCTION stateDynamics
  // exponentiate leading parameters
  FMSY = mfexp ( lnFMSY );         // optimal fishing mortality
  BMSY = mfexp ( lnBMSY );           // MSY
  MSY = elem_prod( BMSY, FMSY);     // optimal biomass

  // initialise variables
  Bt_bar = 0.; chol = 0.; epstCorr = 0.; zetatCorr = 0.;

  // Compute cholesky factor of correlation matrix
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
    Bt_bar(s,1) = ( 2. * BMSY(s) ) * mfexp(epstCorr(1)+zetatCorr(s,1));
    // loop over time periods to build dynamic model for t > 1
    for(int t=1; t<nT; t++)
    {
      pospen = 0.;
      // Run state dynamics equation, add process error
      Bt_bar(s,t+1) = ( Bt_bar(s,t) + 
                      2. * FMSY(s) * Bt_bar(s,t) * (1 - Bt_bar(s,t)/BMSY(s)/2.0 ) - 
                      katch(s,t) ) * mfexp ( epstCorr(t+1)+zetatCorr(s,t+1) );
      Bt_bar(s,t+1) = posfun ( Bt_bar(s,t+1), 10e-1, pospen );
      
      // Increment function penaliser variable
      fpen += 1000. * pospen;
    }
  }
  
  // compute derived management quantities
  dep_bar = elem_div(column(Bt_bar, nT ) ,BMSY) / 2;     // depletion estimate
  
// function to compute predicted observations and residuals
FUNCTION obsModel
  for ( int s=1;s<=nS;s++)
  {
    // Compute conditional estimate of lnq_hat
    lnqhat(s) = sum ( log ( It(s) ) - log ( Bt_bar(s) ) ) / nT ;

    // Compute predicted index of abundance
    It_bar(s) = mfexp(lnqhat(s)) * Bt_bar(s);

    // Compute the sum of squared residuals
    SSRobs(s) = norm2 ( log ( It(s) ) - log ( It_bar(s) ) );
  }
// Subroutine to compute negative log likelihoods
FUNCTION calcLikelihoods
  // Initialise NLL variables
  obsNLL = 0.; procNLL = 0.; nll = 0.; envNLL=0.;
  // Compute observation error conditional variance
  tau2hat = SSRobs/nT;
  // compute observation model likelihood
  obsNLL = 0.5*nT*log(SSRobs);
  // Compute proc error conditional variance, and proc error likelihood
  for ( int s=1;s<=nS;s++)
  {
    Sigma2hat(s) = norm2(zetatCorr(s))/nT;
    procNLL(s) = (nT)*log(norm2(zetatCorr(s)))/2.;
  }
  // Compute environmental error conditional variance and likelihood
  sigma2hat = norm2(epst)/nT;
  envNLL = 0.5*nT*log(norm2(epstCorr));
  // sum likelihoods
  nll = obsNLL + procNLL;

// Subroutine to compute prior densities to add to posterior
FUNCTION calcPriors
  // Initialise prior var
  prior = 0.;
  // First, MSY
  prior = elem_div(pow ( BMSY - mBMSY, 2 ), pow(sBMSY,2)) / 2.;
  // Then FMSY
  prior += elem_div(pow ( FMSY - mFMSY, 2 ), pow(sFMSY,2)) / 2.;
  // Now Sigma2hat prior
  prior += elem_prod((alphaSigma+1),log(Sigma2hat))+elem_div(betaSigma,Sigma2hat);
  // sigma2hat prior
  prior += (alpha_sigmaMS+1)*log(sigma2hat) + beta_sigmaMS/sigma2hat;
  // tau2hat prior
  prior += elem_prod((alphaTau+1),log(tau2hat))+elem_div(betaTau,tau2hat);
  // lnq prior
  prior += pow( lnqhat - mlnq, 2 ) / slnq / slnq / 2;


FUNCTION mcDumpOut
  if( mcHeader == 0 )
  {
    // Output the parameters needed for performance testing of SA method.
    mcoutpar << "BMSY MSY FMSY q Sigma2 tau2 FnT_bar BnT dep_bar sigma2" << endl;
  
    for ( int s=1;s<=nS;s++)
    {
      // Output the parameter values for this replicate
      mcoutpar << BMSY(s) << " " << MSY(s) <<" "<< FMSY(s) <<" "<< mfexp(lnqhat(s));
      mcoutpar <<" "<< Sigma2hat(s) <<" "<< tau2hat(s) <<" " << FnT_bar(s) <<" ";
      mcoutpar << Bt_bar(s,nT) <<" "<< dep_bar(s) << " " << sigma2hat << endl;

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
        mcoutpar <<" "<< Sigma2hat(s) <<" "<< tau2hat(s) <<" " << FnT_bar(s) <<" ";
        mcoutpar << Bt_bar(s,nT) <<" "<< dep_bar(s) << " " <<  sigma2hat << endl;

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
  report << "# BMSY" << endl;
  report << BMSY << endl;
  report << "# FMSY" << endl;
  report << FMSY << endl;
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
  report << "# MSY" << endl;
  report << MSY << endl;
  report <<"# q" << endl;
  report << mfexp(lnqhat) <<endl;
  report << "# Sigma2" << endl;
  report << Sigma2hat << endl;
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
  report << "# epstCorr" << endl;
  report << epstCorr << endl;
  report << "# zetatCorr" << endl;
  report << zetatCorr << endl;
  report << "# chol" << endl;
  report << chol << endl;
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
  report << "# alpha.sigma" << endl;
  report << alpha_sigmaMS << endl;
  report << "# beta.sigma" << endl;
  report << beta_sigmaMS << endl;
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



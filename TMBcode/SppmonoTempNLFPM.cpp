
#include <TMB.hpp> 
#include <iostream>

template<class Type> //only list this in whichever cpp file is alphabetically first since thats the first cpp file read in and the others aren't include guarded meaning it will tell you that it has been re-defined
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;
  Type nll = 0.0;
  Type one = 1.0;
  
  //READ IN THE INPUT DATA
  DATA_INTEGER(n);
  DATA_INTEGER(nyrs);
  // DATA_INTEGER(ndex);
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(idex);
  DATA_IVECTOR(sppex);
  DATA_SCALAR(k);
  DATA_VECTOR(pa);
  DATA_IVECTOR(idmod);
  //TEMPERATURE MODIFIER
  DATA_VECTOR(temp);
  
  
  //READ THE PARAMETERS
  PARAMETER_VECTOR(iye);
  PARAMETER(logrw_var);  
  PARAMETER_VECTOR(lchi);
  PARAMETER_VECTOR(lbeta);
  PARAMETER_VECTOR(lmink);
  PARAMETER_VECTOR(slo);
  PARAMETER_VECTOR(lmidp);
  
  
  //TRANSFORM PARAMETERS AND CREATE DERIVED PARAMETERS
  Type rw_var = exp(logrw_var);
  vector<Type> beta = exp(lbeta);
  vector<Type> chi = exp(lchi);
  
  
  // SCALING PARAMETER FOR TEMPERATURE MODIFICATION
  vector<Type> mink(lmink.size());
  for(int i = 0; i < lmink.size(); i++){
    if(idmod(i) == 0){mink(i) = k;}
    if(idmod(i) == 5){mink(i) = 1/(1+exp(-lmink(i)));}
  }
  
  // vector<Type> slo(lslo.size());
  // for(int i = 0; i < lslo.size(); i++){
  //   if(idmod(i) == 0){slo(i) = k;}
  //   if(idmod(i) == 5){slo(i) = exp(lslo(i));}
  // }
  
  vector<Type> midp(lmidp.size());
  for(int i = 0; i < lmidp.size(); i++){
    if(idmod(i) == 0){midp(i) = k;}
    if(idmod(i) == 5){midp(i) = exp(lmidp(i));}
  }
  
  
  // SET UP ENV FOR VARIABLES USED IN FUNCTIONAL RESPONSE EQUATION
  vector<Type> mu(nyrs);
  vector<Type> new_mu(n);
  vector<Type> theta(n);
  
  // SET UP ENV FOR RESIDUALS
  vector<Type> sto_resid(n);
  vector<Type> sgn_resid(n);
  vector<Type> dev_resid(n);
  
  
  int id, iy, spp;
  
  for(int i = 0; i < n; ++i){
    iy = iyear(i);
    id = idex(i);
    spp = sppex(i);
    
    mu(iy) = one-exp(-iye(iy));  //EXPONENTIAL FOR TRAWL DATA
    
    //TRAWL PROBABILITY UNBIASED
    if(id == 0){      
      if(idmod(i) == 0){
        new_mu(i) = mu(iy);
      }  
      if(idmod(i) == 5){
        theta(i) = 0;  //TEMPERATURE HAS NO EFFECT ON TRAWL BUT NEED TO KEEP INDEXING
        
        new_mu(i) = mu(iy);
      }
    } 
    
    
    //FUNCTIONAL RESPONSE FOR EACH STOMACH CONTENT DATA SET
    if(id > 0){
      if(idmod(i) == 0){
        new_mu(i) = (slo(spp-1)*pow(mu(iy), beta(id-1)))/(((pow((chi(id-1)), beta(id-1))))+pow(mu(iy), beta(id-1))); 
      }
      if(idmod(i) == 5){
        theta(i) = 1+((mink(spp-1)-1)/(1+exp(slo(spp-1)*(temp(i)-midp(spp-1)))));  //TEMPERATURE MODIFIER SIGMOIDAL EQUATION
        
        new_mu(i) = (theta(i)*(pow(mu(iy), beta(id-1))))/(((pow((chi(id-1)), beta(id-1))))+pow(mu(iy), beta(id-1)));  //WITH TEMPERATURE MODIFIER
      } 
    }
    
    
    if(isNA(pa(i)) == false){
      nll -= dbinom(pa(i), one, new_mu(i), true); //BERNOULLI LIKELIHOOD FOR PRESENCE/ABSENCE DATA
      
      //SIMULATE RESIDUALS
      sto_resid(i) = pa(i) - new_mu(i);
      sgn_resid(i) = (sto_resid(i) > 0) - (sto_resid(i) < 0);
      dev_resid(i) = sgn_resid(i)*sqrt(-2.0*((pa(i)*log(new_mu(i)))+(1-pa(i))*log(1.0-new_mu(i))));
      
      SIMULATE{
        pa(i) = rbinom(one, new_mu(i));
      }
    }
    
  }
  
  
  //GAUSSIAN RW FOR THE ABUNDANCE INDEX
  vector<Type> del_iye = log(iye);
  nll -= dnorm(del_iye(0), Type(10.0), rw_var, true);
  
  for(int i = 1; i < nyrs; ++i){
    nll -= dnorm(del_iye(i), del_iye(i-1), rw_var, true);
  }
  
  
  REPORT(iyear);
  REPORT(idex);
  REPORT(rw_var);
  REPORT(mu);
  REPORT(idmod);
  REPORT(temp);
  REPORT(theta);
  
  REPORT(sto_resid);
  REPORT(dev_resid);
  
  REPORT(iye);  
  REPORT(new_mu);  
  REPORT(beta);
  REPORT(chi);
  REPORT(mink);
  REPORT(slo);
  REPORT(midp);
  
  
  ADREPORT(iye);
  ADREPORT(beta);
  ADREPORT(chi);
  ADREPORT(mink);
  ADREPORT(slo);
  ADREPORT(midp);
  
  return nll;
}


/*
 * logisticGrowth.cpp: TMB modeling interface, links TMB to model.hpp
 * 
 * This example is not a module of FIMS and is for learning purposes only 
 * to demonstrate the concepts of modularity and portability. 
 */
#define TMB_MODEL //defines TMB_MODEL: Common.hpp sets up TMB specific preprocessing macros 
#include "model.hpp" //contains the model
#include "DerivedQuantities.hpp" //functions to calculate derived quantities

//set-up TMB objective function type
template<class Type>
Type objective_function<Type>::operator()(){
  /*
   * create pointer, inst, that points to singleton class of logisticGrowth model
   * getinstance is defined in model.hpp
   */ 
  logisticGrowth<Type>* inst = logisticGrowth<Type>::getinstance();
  
  //Set up TMB macros for data and parameters
  DATA_VECTOR(y);
  
  PARAMETER_VECTOR(theta);
  PARAMETER(ln_sig);
  PARAMETER(ln_tau);
  PARAMETER_VECTOR(u);
  
  /*
   * access and assign members of logisticGrowth class using inst pointer
   */
  
  inst->y = y;
  inst->ln_sig = ln_sig;
  inst->ln_tau = ln_tau;
  inst->u = u;
  
  Type r = exp(theta[0]);
  Type K = exp(theta[1]);
  inst->r = r;
  inst->K = K;
  
  Type sigma = exp(ln_sig);
  Type tau = exp(ln_tau);
  
  /*
   *   create Type nll and assign value to the return of the 
   *   evaluate() function defined in model.hpp
   */   
  Type nll = inst -> evaluate();
  //Use TMB's simulate function to generate data
  SIMULATE{
    /*
     * calculateEta is a member function inside the inst object (instance of logisticGrowth class) 
       so needs to be accessed through the pointer, inst
     */
    vector<Type> eta = inst -> calculateEta(inst ->u, inst ->r, inst ->K);
    u = exp(rnorm(log(eta), sigma));
    y = exp(rnorm(log(u), tau));
    REPORT(u);
    REPORT(y);
  }
  //REPORT out quantities
  REPORT(r);
  REPORT(K);
  ADREPORT(r);
  ADREPORT(K);
  REPORT(u);

  
  REPORT(sigma);
  REPORT(tau);
  ADREPORT(sigma);
  ADREPORT(tau);

  //calculate any derived quantities and report out
  Type H = MSY(r,K);
  REPORT(H);
  ADREPORT(H);
  
  return nll;

}


GLOBALS_SECTION
  #define ADMB_MODEL
  #include "../model.hpp"

DATA_SECTION
  init_int n
  init_vector y(1,n)
PARAMETER_SECTION
  init_number ln_sig
  init_number ln_tau
  init_bounded_number ln_r(-5,-0.001);
  init_bounded_number ln_K(3.5,7);
  sdreport_number r;
  sdreport_number K;
  sdreport_number sigma;
  sdreport_number tau;
  random_effects_bounded_vector u(1,n,1,150);
  objective_function_value f
LOCAL_CALCS

  r = exp(ln_r);
  K = exp(ln_K);

  logisticGrowth<dvariable>* instance = logisticGrowth<dvariable>:: getinstance();

  for(int i =1; i <= n; i++){

    instance->y.push_back(this->y[i]);
    instance->u.push_back(this->u[i]);
    
  }
  //set parameter values
  instance->ln_sig = this->ln_sig;
  instance->ln_tau = this->ln_tau;
  instance->r = this->r;
  instance->K = this->K;

PROCEDURE_SECTION

  //get singleton
  logisticGrowth<dvariable>* instance = logisticGrowth<dvariable>:: getinstance();

  //set re values
  for(int i =1; i <= n; i++){
     instance->u[i-1] = this->u[i];
  }
  //set parameter values
  instance->ln_sig = this->ln_sig;
  instance->ln_tau = this->ln_tau;
  instance->r = this->r;
  instance->K = this->K;
  //evaluate the objective function
  this->f = instance->evaluate();
    


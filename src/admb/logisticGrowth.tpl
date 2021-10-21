DATA_SECTION
  init_int n
  init_vector y(1,n)
  
PARAMETER_SECTION
  init_bounded_number ln_r(-5,-.001)
  init_bounded_number ln_K(3.5,7)
  init_number ln_sig
  init_number ln_tau
  sdreport_number r;
  sdreport_number K;
  sdreport_number sigma;
  sdreport_number tau;

  random_effects_bounded_vector u(1,n,2,150,2);

  objective_function_value jnll


PROCEDURE_SECTION
  r = exp(ln_r); 
  K = exp(ln_K);
  sigma = exp(ln_sig);
  tau = exp(ln_tau);
  for(int i=2; i<=n; ++i){
    step(u(i-1),u(i),r,K,sigma);  }

  for(int i=1; i<=n; ++i){ 
    obs(u(i),tau,i);}

SEPARABLE_FUNCTION void step(const dvariable& x1, const dvariable& x2, const dvariable& r, const dvariable& K, const dvariable& sigma)
  dvariable var=pow(sigma,2);
  dvariable m=log(x1 + r * x1 * (1.0 - x1/K)); 
  jnll+=0.5*(log(2.0*M_PI*var)+square(log(x2)-m)/var) + log(x2);

SEPARABLE_FUNCTION void obs(const dvariable& x, const dvariable& tau, int i)
  dvariable var=pow(tau,2);
  jnll+=0.5*(log(2.0*M_PI*var)+square(log(y(i))-log(x))/var) + log(y(i));

REPORT_SECTION
  report << "r =" << r << endl;
  report << "K =" << K << endl;
  report << "sigma =" << sigma << endl;
  report << "tau =" << tau << endl;
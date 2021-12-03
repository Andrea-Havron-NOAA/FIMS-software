DATA_SECTION
  init_int n
  init_vector y(1,n)
  
PARAMETER_SECTION
  init_bounded_number ln_r(-5,-0.001)
  init_bounded_number ln_K(3.5,7)
  init_number ln_sig
  init_number ln_tau
  sdreport_number r;
  sdreport_number K;
  sdreport_number sigma;
  sdreport_number tau;

  random_effects_bounded_vector u(1,n,1,150);

  objective_function_value jnll


PROCEDURE_SECTION
  r = exp(ln_r); 
  K = exp(ln_K);
  sigma = exp(ln_sig);
  tau = exp(ln_tau);
  for(int i=2; i<=n; ++i){
    step(u(i-1),u(i),ln_r,ln_K,ln_sig);  }

  for(int i=1; i<=n; ++i){ 
    obs(u(i),ln_tau,i);}

SEPARABLE_FUNCTION void step(const dvariable& x1, const dvariable& x2, const dvariable& lnr, const dvariable& lnK, const dvariable& lnsig)
  dvariable var=pow(exp(lnsig),2);
  dvariable m=log(x1 + exp(lnr) * x1 * (1.0 - x1/exp(lnK))); 
  jnll+=0.5*(log(2.0*M_PI*var)+square(log(x2)-m)/var) + log(x2);

SEPARABLE_FUNCTION void obs(const dvariable& x, const dvariable& lntau, int _i)
  dvariable var=pow(exp(lntau),2);
  jnll+=0.5*(log(2.0*M_PI*var)+square(log(y(_i))-log(x))/var) + log(y(_i));

REPORT_SECTION
  report << "r =" << r << endl;
  report << "K =" << K << endl;
  report << "sigma =" << sigma << endl;
  report << "tau =" << tau << endl;
  
TOP_OF_MAIN_SECTION
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(4604040);
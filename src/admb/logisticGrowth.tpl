DATA_SECTION
  init_int n
  init_vector y(1,n)
  
PARAMETER_SECTION
  init_bounded_number ln_r(-5,-.001)
  init_bounded_number ln_K(3.5,7)
  init_number ln_sig
  init_number ln_tau

  random_effects_vector u(1,n);

  objective_function_value jnll

PROCEDURE_SECTION
  for(int i=2; i<=n; ++i){
    step(u(i-1),u(i),ln_r,ln_K,ln_sig);  }

  for(int i=1; i<=n; ++i){ 
    obs(u(i),ln_tau,i);}

SEPARABLE_FUNCTION void step(const dvariable& x1, const dvariable& x2, const dvariable& ln_r, const dvariable& ln_K, const dvariable& ln_sig)
  dvariable var=exp(ln_sig);
  dvariable m=log(x1 + exp(ln_r) * x1 * (1.0 - x1/exp(ln_K))); 
  jnll+=0.5*(log(2.0*M_PI*var)+square(log(x2)-m)/var) + log(x2);

SEPARABLE_FUNCTION void obs(const dvariable& x, const dvariable& ln_tau, int i)
  dvariable var=exp(ln_tau);
  jnll+=0.5*(log(2.0*M_PI*var)+square(log(y(i))-log(x))/var) + log(y(i));

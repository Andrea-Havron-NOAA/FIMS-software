DATA_SECTION
  init_int N
  init_vector y(1,N)

PARAMETER_SECTION
  init_number logr0
  init_number logpsi
  init_bounded_number logK(4.6,7.6)
  init_number logQ
  init_number logR

  random_effects_vector u(1,N);

  objective_function_value jnll
  			   
PROCEDURE_SECTION
  for(int i=2; i<=N; ++i){
    step(u(i-1),u(i),logr0,logK,logpsi,logQ);  }

  for(int i=1; i<=N; ++i){ 
    obs(u(i),logR,i);}

SEPARABLE_FUNCTION void step(const dvariable& x1, const dvariable& x2, const dvariable& logr0, const dvariable& logK, const dvariable& logpsi, const dvariable& logQ)
  dvariable var=exp(logQ);
  dvariable m=x1 + exp(logr0) * (1.0 - pow(exp(x1)/exp(logK),exp(logpsi))); 
  jnll+=0.5*(log(2.0*M_PI*var)+square(x2-m)/var);

SEPARABLE_FUNCTION void obs(const dvariable& x, const dvariable& logR, int i)
  dvariable var=exp(logR);
  jnll+=0.5*(log(2.0*M_PI*var)+square(x-y(i))/var);

DATA_SECTION
  init_int n
  init_vector y(1,n)
  
PARAMETER_SECTION
  init_bounded_number ln_r(-5,-0.001)
  init_bounded_number ln_K(3.5,7)
  init_number ln_sig(2)
  init_number ln_tau(2)
  sdreport_number r;
  sdreport_number K;
  sdreport_number sigma;
  sdreport_number tau;
  sdreport_number foo1;
  sdreport_number foo2;

  random_effects_bounded_vector u(1,n,1,150);

  objective_function_value jnll


PROCEDURE_SECTION
  r = exp(ln_r); 
  K = exp(ln_K);
  sigma = exp(ln_sig);
  tau = exp(ln_tau);
  dvariable tau2 = pow(tau,2);
  dvariable sig2 = pow(sigma,2);
  
  for(int i=2; i<=n; ++i){
    dvariable m=log(u(i-1) + r * u(i-1) * (1.0 - u(i-1)/K)); 
    jnll+=0.5*(log(2.0*M_PI*sig2)+square(log(u(i))-m)/sig2) + log(u(i));
    foo1 = i;
  }

  for(int i=1; i<=n; ++i){ 
    jnll+=0.5*(log(2.0*M_PI*tau2)+square(log(y(i))-log(u(i)))/tau2) + log(y(i));
    foo2 = i;
  }


REPORT_SECTION
  report << "r =" << r << endl;
  report << "K =" << K << endl;
  report << "sigma =" << sigma << endl;
  report << "tau =" << tau << endl;
  report << "foo1 =" << foo1 << endl;
  report << "foo2 =" << foo2 << endl;
  report << "u = " << u << endl;
  
TOP_OF_MAIN_SECTION
  arrmblsize = 80000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(460404);
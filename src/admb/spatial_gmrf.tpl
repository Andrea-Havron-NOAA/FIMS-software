// admb spatial.tpl modified to GMRF based on precision matrix, Q


DATA_SECTION
  init_int n			// Number of observations
  init_int m     //Number of random effects
  init_vector y(1,n)		// Poisson counts
  init_ivector v_i(1,n)
  init_matrix M0(1,m,1,m)
  init_matrix M1(1,m,1,m)
  init_matrix M2(1,m,1,m)
  init_number log_kappa
  init_number log_tau


PARAMETER_SECTION
  init_number b(1)
  random_effects_vector u(1,m,2)

  objective_function_value l

PRELIMINARY_CALCS_SECTION
  //cout << setprecision(4);
  
PROCEDURE_SECTION
  int i,j;
  
  gmrf(m, u);
  for (i=1;i<=n;i++)
    pois_loglik(i,u(v_i(i)),b);			// Likelilhood contribution from i'th point

SEPARABLE_FUNCTION void pois_loglik(int i, const dvariable& ui, const dvariable& _b)
    dvariable eta = _b + ui/exp(log_tau);	// Linear predictor
    dvariable lambda = mfexp(eta);			// Mean in Poisson distribution
    l += lambda-y(i)*eta;
    
SEPARABLE_FUNCTION void gmrf(int _m, const dvar_vector& _u)  
  dvar_matrix Ucol(1,m,1,1);
  dvar_matrix Urow(1,1,1,m);
  for(int j=1; j<=m; j++){
     Ucol(j,1) = _u(j);
     Urow(1,j) = _u(j);
  }
   dvar_matrix Q(1,_m,1,_m);
   dvariable ln_detQ;

   Q = pow(exp(log_kappa),4) * M0 + 2 * pow(exp(log_kappa),2) * M1 + M2;
   ln_detQ = ln_det(Q);
   l +=  0.5 * sum(-ln_detQ + Urow*Q*Ucol ); 


REPORT_SECTION
  report << "b =" << b << endl;
  //report << "tau =" << exp(log_tau) << endl;
  //report << "kappa =" << exp(log_kappa) << endl;
  //report << "spatial var =" << pow(1/(sqrt(2*M_PI)*exp(log_kappa)*exp(log_tau)),2) << endl;
  //report << "range =" << sqrt(8)/exp(log_kappa) << endl;
  report << "u =" << u << endl;
  

TOP_OF_MAIN_SECTION
  arrmblsize=40000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  gradient_structure::set_MAX_NVAR_OFFSET(2000000);

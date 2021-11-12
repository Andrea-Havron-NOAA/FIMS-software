// Poisson regression with spatially smoothed mean.

DATA_SECTION
  init_int n			// Number of observations
  init_int m     //Number of random effects
  init_vector y(1,n)		// Poisson counts
  init_ivector v_i(1,n)
  init_matrix M0(1,m,1,m)
  init_matrix M1(1,m,1,m)
  init_matrix M2(1,m,1,m)
  
  //init_int p			// Number of fixed effects (b's)
  //init_matrix X(1,n,1,p)	// Covariate matrix
  //init_matrix Z(1,n,1,2)	// Locations (xy-coordinates)
  //matrix dd(1,n,1,n);		// Distance matrix


PARAMETER_SECTION
  init_bounded_number log_kappa(-3, 10)
  init_bounded_number log_tau(-3,10)
  init_number b
  random_effects_vector u(1,m)
  //not sure what this does
  normal_prior M(u);
  
  
  //init_bounded_vector b(1,p,-100,100)  
  //init_bounded_number a(0.01,3.0)
  //init_bounded_number log_sigma(-3.0,3.0)

  objective_function_value l

PRELIMINARY_CALCS_SECTION
  //cout << setprecision(4);
  
 
  
PROCEDURE_SECTION

  int i;
  for (i=1;i<=n;i++)
    pois_loglik(i,u(v_i(i)),b,log_tau);			// Likelilhood contribution from i'th point

SEPARABLE_FUNCTION void pois_loglik(int i,const dvariable& ui,const dvar_vector& _b,const dvariable& _log_tau)
    dvariable eta = b + ui/exp(_log_tau);	// Linear predictor
    dvariable lambda = mfexp(eta);			// Mean in Poisson distribution
    l += lambda-y(i)*eta;
    
   

   ///////////////rewrite
              NORMAL_PRIOR_FUNCTION void get_M(const dvariable& log_kappa, const dvariable& log_tau )
                int i,j;
                dvar_matrix Q(1,m,1,m);
                //do I need to write this as a loop?
                Q = pow(exp(log_kappa),4) * M0 + 2 * pow(exp(log_kappa),2) * M1 + M2;
                // add the GMRF likelihood - is this correct? 
                l += pow( det(Q), 0.5 ) * mfexp(-0.5*omega*Q*omega) //nned to convert equation to ADMB
                
                
                 //NORMAL_PRIOR_FUNCTION void get_M(const dvariable& a)
                //int i,j;
                //dvar_matrix tmpM(1,n,1,n);
                //dvar_matrix tmpM(1,n,1,n);
                //for (i=1;i<=n;i++)
                //{
                //  tmpM(i,i)=1.0;
                //  for ( j=1;j<i;j++)
                //  {
                //    tmpM(i,j)=exp(-a*dd(i,j));			// Exponentially decaying correlation
                //    tmpM(j,i)=tmpM(i,j);
                //  }
                //}
                //M=tmpM;
              
              //FUNCTION void evaluate_M(void)
              
                //get_M(a);
  ////////////////////

REPORT_SECTION
  cout << "b =" << b << "tau=" << exp(log_tau) << endl;
  cout << "b =" << b << "kappa=" << exp(log_kappa) << endl;

TOP_OF_MAIN_SECTION
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(460404);

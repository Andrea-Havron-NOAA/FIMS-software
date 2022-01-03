#undef  TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR of    
template<class Type>
  vector<Type> test(vector<Type> mu, Type sd, objective_function<Type>* of){
    vector<Type> y( mu.size() );
    SIMULATE {
      y = rnorm(mu, sd);  // Simulate response
    }
    
    return y;
  }
#undef  TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this



#include "Common.hpp"
#include "simFns.hpp"

template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_VECTOR(y);
    DATA_VECTOR(x);
    
    PARAMETER(a);
    PARAMETER(b);
    PARAMETER(sd);
    
    vector<Type> mu = a + b * x;
    Type f = -sum(dnorm(y, mu, sd, true));
    y = test( mu, sd, this);
    
    REPORT(y);
    return f;
  }

#ifndef MODEL_HPP
#define MODEL_HPP
#include <vector>
#include "Common.hpp"
#include "likelihoods.hpp"
  
  template<class Type>
  class logisticGrowth{
  public:
    // /* Data section */
    typename model_traits<Type>::data_vector  y;
    // /* Parameter section */
    typename model_traits<Type>::parameter_vector u;
    typename model_traits<Type>::parameter_vector theta;
    Type ln_sig;
    Type ln_tau;
    static logisticGrowth<Type>* instance;

    
    logisticGrowth(){}
    
    static logisticGrowth<Type>* getinstance(){
      return logisticGrowth<Type>::instance;
    }
    
    Type evaluate(){
      Type r = exp(theta[0]);
      Type K = exp(theta[1]);
      Type sigma = exp(ln_sig);
      Type tau = exp(ln_tau);
      int t;
      int n=y.size();
      Type nll = 0;
      typename model_traits<Type>::parameter_vector eta(n);
      for(t=1; t<n; t++){
        eta(t) = u[t-1] + r * u[t-1] * (1-u[t-1]/K);
        nll -= dlognorm(u[t], log(eta[t]), sigma, true);
      }
      for(t=0; t<n; t++){
        nll -= dlognorm(y[t],log(u[t]),tau,true);
      }
 
      return nll;
    }

  };
  
  template<class Type>
  logisticGrowth<Type>* logisticGrowth<Type>::instance = new logisticGrowth<Type>();

#endif

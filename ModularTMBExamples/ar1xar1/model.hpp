#ifndef MODEL_HPP
#define MODEL_HPP
#include <vector>
#include "Common.hpp"
  
  
template<class Type>
class ar1xar1{
  
  using DataVector = typename model_traits<Type>::data_vector;
  using DataIndicator = typename model_traits<Type>::data_indicator;
  using ParameterArray = typename model_traits<Type>::array;
  //using ParameterVector = typename model_traits<Type>::parameter_vector;
  //using IndicatorVector = typename model_traits<Type>::indicator_vector;
  //using ParameterArray = typename model_traits<Type>::parameter_array;
  
public:
  // /* Data section */
  DataVector  y;
  DataIndicator keep;
  //data_indicator(keep,y);
  // /* Parameter section */
  ParameterArray eta;
  Type phi1;
  Type phi2;
  static ar1xar1<Type>* instance;
  
  #ifdef TMB_MODEL
  objective_function<Type>* of;
  #endif

  ar1xar1(){}
  
  static ar1xar1<Type>* getinstance(){
    return ar1xar1<Type>::instance;
  }
  
 
  Type evaluate(){
    Type nll = 0;
    vector<Type> ln_lambda(y.size());
    nll += SEPARABLE( AR1(phi2), AR1(phi1) )(eta);
    
    #ifdef TMB_MODEL
    SIMULATE_F(of){
      SEPARABLE(AR1(phi2), AR1(phi1)).simulate(eta);
      vector<Type> ln_lambda = eta;//convert array to vector for rpois simulation
      y = rpois(exp(ln_lambda));//rpois cannot accept array 
    }

    SIMULATE_F(of){
      REPORT_F(y, of);
    }
    #endif
    
    for(int i=0; i < y.size(); i++){
      nll -= keep[i] * dpois(y[i], exp(eta[i]), true);
      Type cdf = squeeze( ppois(y[i], exp(eta[i])) );
      nll -= keep.cdf_lower[i] * log( cdf );       // NaN protected
      nll -= keep.cdf_upper[i] * log( 1.0 - cdf ); // NaN protected
    }

    #ifdef TMB_MODEL
    REPORT_F(phi1, of);
    REPORT_F(phi2, of);
    REPORT_F(eta, of);
    ADREPORT_F(phi1, of);
    ADREPORT_F(phi2, of);
    #endif

    return nll;
  }
  
};

template<class Type>
ar1xar1<Type>* ar1xar1<Type>::instance = new ar1xar1<Type>();

#endif

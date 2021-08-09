#ifndef MODEL_HPP
#define MODEL_HPP
#include <vector>
#include "Commom.hpp"

template<class Type>
class ThetaLog{
public:
  // /* Data section */
  typename model_traits<Type>::data_vector  Y;
  // /* Parameter section */
  typename model_traits<Type>::parameter_vector X;
  Type logr0;
  Type logtheta;
  Type logK;
  Type logQ;
  Type logR;
  static ThetaLog<Type>* instance;
  
  
  ThetaLog(){}
  
  static ThetaLog<Type>* get_instance(){
    return ThetaLog<Type>::instance;
  }
  
  Type evaluate(){
    Type r0=exp(logr0);
    Type theta=exp(logtheta);
    Type K=exp(logK);
    Type Q=exp(logQ);
    Type R=exp(logR);
    int timeSteps=Y.size();
    Type ans=0;
    for(int i=1;i<timeSteps;i++){
      Type m=X[i-1]+r0*(1.0-pow(exp(X[i-1])/K,theta));
      ans-=dnorm(X[i],m,sqrt(Q),true);
    }
    for(int i=0;i<timeSteps;i++){
      ans-=dnorm(Y[i],X[i],sqrt(R),true);
    }
    return ans;
  }

};

template<class Type>
ThetaLog<Type>* ThetaLog<Type>::instance = new ThetaLog<Type>();

#endif
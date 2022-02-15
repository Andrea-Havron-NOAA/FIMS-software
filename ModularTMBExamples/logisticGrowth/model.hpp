/*
 * model.hpp defines the model to be optimized. The model is set up as a 
 * C++ templated class reliant on base C++ functions. By minimizing dependencies,
 * this model can be passed around to other software platforms as needed.
 * 
 * This example is not a module of FIMS and is for learning purposes only 
 * to demonstrate the concept of modularity and portability. 
 */


#ifndef MODEL_HPP
#define MODEL_HPP
#include <vector>
#include "Common.hpp" //sets up pre-processing macros
#include "likelihoods.hpp" //reads in distribution functions
  

/*
 * The logisticGrowth class is based on templated C++ programming. This 
 * essentially means that Types can be generic and are determined by C++
 * during the compilation stage based on data input.
 */
template<class Type>
class logisticGrowth{
  
  /*
  * using: C++ keyword used here to declare class members
  * Lines 29 and 30 declare DataVector and ParameterVector and assign
  * values using the model_traits from Common.hpp
  */
  using DataVector = typename model_traits<Type>::data_vector;
  using ParameterVector = typename model_traits<Type>::parameter_vector;
  
  /*
   * public: members of class that are available outside of the class
   */
    
public:
  // /* Data section */
  DataVector  y;
  // /* Parameter section */
  ParameterVector u;
  Type r;
  Type K;
  Type ln_sig;
  Type ln_tau;
  /*  
   * static: A static object is an object that persists from the time it's constructed 
              until the end of the program.
   * instance: The * symbol means we are initiating a pointer called 'instance'
   */ 
  static logisticGrowth<Type>* instance;

  /* Constructor section: 
   * The constructor is a special member function having the same name as that of its class
      and is used to initialize values to data/parameter members. C++ executes the constructor 
      automatically when an object of the class is created.
      Here, the constructor is empty as the logisticGrowth class is 
      initalized post construction when called in the TMB interface
  */
  logisticGrowth(){}
  
  /* 
   * Create static member function called getinstance() that initialize a singleton object
   * of logisticGrowth if needed, else return the singleton object
   * new: create object that is stored in heap memory (dynamically allocated)  
   * A singleton class is a way to insure the class is only created once 
   * - important for memory management
   */
  static logisticGrowth<Type>* getinstance(){
    if(logisticGrowth<Type>::instance == NULL){
      logisticGrowth<Type>::instance = new logisticGrowth<Type>();
    }
    return logisticGrowth<Type>::instance;
  }
  
  //calculateEta() 
  /**
   * @brief calculates mean random effects values
   * 
   * @param u random effect vector
   * @param r growth rate
   * @param K carrying capacity
   * @return Returns a vector of mean random effect values
   */
  ParameterVector calculateEta(ParameterVector u, Type r, Type K){
    
    ParameterVector eta(u.size());
    
    for(int t=1; t<u.size(); t++){
      eta[t] = u[t-1] + r * u[t-1] * (1-u[t-1]/K);
    }
  
    return eta; 
  }
  
  //calculateNLL()
  /**
   * @brief function that calculates the negative log likelihoods of observations and random effects
   * 
   * @param u random effects vector
   * @param eta mean of random effects vector
   * @param sigma process standard deviation
   * @param tau observation standard deviation
   * @param y observations
   * @return returns joint negative log likelihood of random effects and observations
   */
  Type calculateNll(ParameterVector u, ParameterVector eta, Type sigma, Type tau, DataVector y){
  
    Type nll = 0.0; 
    int n=y.size();
    
    for(int t=1; t<n; t++){
      nll -= dlognorm<Type>(u[t], log(eta[t]), sigma, true);
    }
    
    for(int t=0; t<n; t++){
      nll -= dlognorm<Type>(y[t], log(u[t]), tau, true);
    }
  
    return nll;  
  }
  
  //evaluate()
  /**
   * @brief  main function that runs the model. 
   * 
   * @return joint negative log likelihood to be minimized
   */
  Type evaluate(){
    Type sigma = exp(ln_sig);
    Type tau = exp(ln_tau);
    //auto is used to automatically figure out the Type of eta that is defined by calculateEta
    auto eta = calculateEta(u, r, K); 
    Type nll = calculateNll(u, eta, sigma, tau, y);

    return nll;
  }

};

//static members are defined outside class definition
template<class Type>
logisticGrowth<Type>* logisticGrowth<Type>::instance = NULL;

#endif

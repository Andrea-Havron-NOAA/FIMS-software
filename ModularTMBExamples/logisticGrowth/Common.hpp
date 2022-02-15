/*
 * Common.hpp sets up pre-processing macros that are run before compilation.
 * 
 * This example is not a module of FIMS and is for learning purposes only 
 * to demonstrate the concepts of modularity and portability. 
 */

/* 
 * #ifndef <token>: C++ header guard
 * checks if given token has been defined, if not defined, runs all code
 * between #ifndef and #endif
 */

//Lines 8-9: if COMMON_HPP not defined, define using code below
#ifndef COMMON_HPP 
#define COMMON_HPP 

/*
 * Based on whether we are running a TMB or ADMB model,
 * we will include different headers and type definitions specific to the 
 * software platform requirements. 
 * /
 
  

 * data_vector and parameter_vector (members of the struct model_traits) 
 * dependent on the requirements of TMB or ADMB
 */

/*if TMB_MODEL defined, run lines 14-30
 * #define TMB_MODEL is found at the top of the TMB version of the model,
 * logisticGrowth.cpp. If the TMB model is called from R using MakeADFun, the 
 * TMB_MODEL will be defined and this section of Common.hpp will be run.  
*/

#ifdef TMB_MODEL

#include <TMB.hpp>


/* typedef: C++ specifier used to indicate the declaration 
 *   is defining a type rather than a function or variable 
 * typename: C++ keyword used in a template to specify 
 *  that the identifier that follows is a type. */ 
template<typename Type>
/* 
 * Here we are creating a struct called model_traits with memebers 
 * data_vector and parameter_vector defined based on TMB requirements.
 * The data_vector and parameter_vector can be used in model.hpp regardless of
 * how the type is defined.
 */
struct model_traits{
 typedef typename CppAD::vector<Type> data_vector; 
 typedef typename CppAD::vector<Type> parameter_vector;
};

#endif

/* if ADMB_MODEL defined, run lines 35-52
* #define ADMB_MODEL is found at the top of the ADMB version of the model,
* logisticGrowth.tpl. If the ADMB model is called, the 
* ADMB_MODEL will be defined and this section of Common.hpp will be run.  
*/
#ifdef ADMB_MODEL

#include <admodel.h>
#include <vector>

template<typename Type>
/*  Here data_vector and parameter_vector are defined based on ADMB requirements.
 * This method allows data_vector and parameter_vector to be portable. We don't 
 * need to update their type defintions in model.hpp if we swicth to a 
 * different back-end.
 */
struct model_traits{
  typedef typename std::vector<Type> data_vector;
  typedef typename std::vector<Type> parameter_vector;
  typedef Type variable;
};



#endif

#ifdef STD_LIB

#include <cmath>
#include <vector>


template<typename T>
T exp(const T& x){
  return std::exp(x);
}

template <class T>
const T log(const T& x){return std::log(x);}

#endif

#endif

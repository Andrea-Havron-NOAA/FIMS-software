//Derived quantity functions

/**
 * @brief simple MSY calculation
 * 
 * @tparam Type 
 * @param r growth rate
 * @param K carrying capacity
 * @return derived MSY value
 */
template<class Type>
Type MSY(Type r, Type K){
  Type H = K*r/4;
  return H;
}

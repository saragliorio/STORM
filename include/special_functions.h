#pragma once

#include "config.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>

#ifndef SPECIAL_FUNCTIONS_NS
#define SPECIAL_FUNCTIONS_NS

/*
===============================================================================

  NAMESPACE: special_functions

  Collection of numerical routines for special functions used throughout the
  codebase.

  The functions implemented here include:
    - Gamma functions and ratios of Gamma functions
    - Gauss hypergeometric functions (2F1)
    - Confluent hypergeometric functions (1F1 and U)

  REFERENCE: Handbook of Mathematical Functions (Abramowitz and Stegun)
===============================================================================
*/

namespace special_functions
{  
  // Compute Gamma(z) for complex argument z
  void sf_compute_gamma_function_z (const complex_type& z, complex_type& gamma);

  // Compute the ratio Gamma(a+n)/Gamma(a) for complex a
  void sf_compute_gamma_ratio (const complex_type& a, const int n, complex_type& gamma_ratio);

  // Compute the ratio Gamma(a+n)/Gamma(a) for real a
  void sf_compute_gamma_ratio (const float_type& a, const int n, float_type& gamma_ratio);


  /*
  ---------------------------------------------------------------------------
  GAUSS HYPERGEOMETRIC FUNCTION 2F1
  ---------------------------------------------------------------------------
  */

  // Compute the Gauss hypergeometric function 2F1(a,b,c,z) for complex parameters
  void sf_compute_2F1_gauss(const complex_type& a , const complex_type& b , const complex_type& c , const float_type& z, complex_type& result);

  // Compute the Gauss hypergeometric function 2F1(a,b,c,z) for real parameters
  void sf_compute_2F1_gauss(const float_type& a , const float_type& b , const float_type& c , const float_type& z, complex_type& result);
  
  // Generalized evaluation of 2F1 for complex parameters
  void sf_compute_2F1_generalized (const complex_type& a , complex_type& b , const complex_type& c , const float_type& z, complex_type& result);

  // Generalized evaluation of 2F1 for real parameters
  void sf_compute_2F1_generalized (const float_type& a , const float_type& b , const float_type& c , const float_type& z, complex_type& result);
  
  // Special treatment when a or b is a negative integer (finite series)
  void sf_compute_2F1_a_or_b_negative_integer_numbers(const complex_type& a , const complex_type& b , const complex_type& c , const float_type& z, complex_type& result);

  // Same case for real parameters
  void sf_compute_2F1_a_or_b_negative_integer_numbers(const float_type& a , const float_type& b , const float_type& c , const float_type& z, complex_type& result);

  // Special case when (a − b) is an integer
  void sf_compute_2F1_gauss_when_a_minus_b_is_integer(const complex_type& a , const complex_type& b , const complex_type& c , const float_type& z, complex_type& result);
  
  // General interface to compute 2F1(a,b,c,z) for complex parameters
  void sf_compute_2F1(const complex_type& a , complex_type& b ,  complex_type& c , const float_type& z, complex_type& result);

  // General interface to compute 2F1(a,b,c,z) for real parameters
  void sf_compute_2F1(const float_type& a , float_type& b ,  float_type& c , const float_type& z, complex_type& result);
  

  /*
  ---------------------------------------------------------------------------
  CONFLUENT HYPERGEOMETRIC FUNCTIONS
  ---------------------------------------------------------------------------
  */

  // Compute the confluent hypergeometric function 1F1(a,b,z)
  void sf_compute_1F1_gauss (const complex_type& a , const complex_type& b , const complex_type& z, complex_type& result);

  // Compute the confluent hypergeometric function U(a,b,z)
  void sf_compute_confluent_hypergeometric_function(const complex_type& a , const complex_type& b ,const complex_type& z, complex_type& U);

}

#endif
#pragma once

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp> 

/*
===============================================================================
  Global configuration file defining:
    - numerical precision types
    - global constants used throughout the code
    - numerical tolerances and algorithm parameters
    - helper mathematical utilities
    - operator overloads for multiprecision arithmetic

  The code uses Boost multiprecision types to ensure high numerical accuracy
  in the computation of special functions, homogeneous solutions, geodesic
  motion, and scalar energy fluxes.

===============================================================================
*/


/*
-----------------------------------------------------------------------------
MULTIPRECISION NUMERICAL TYPES
-----------------------------------------------------------------------------
*/

/// Complex number with 100 decimal digits of precision (the numebr of digits can be set to any desired value)
typedef boost::multiprecision::cpp_complex<100, boost::multiprecision::backends::digit_base_10 > complex_type;

/// Real floating-point number with 100 decimal digits of precision (the numebr of digits can be set to any desired value)
typedef boost::multiprecision::number < boost::multiprecision::cpp_bin_float <100, boost::multiprecision::backends::digit_base_10 >> float_type;

/// Namespace aliases for convenience
namespace mp = boost::multiprecision;
namespace bm = boost::math;

#ifndef CONSTANTS_NS
#define CONSTANTS_NS


//  ####################################################################
//                        GLOBAL NUMERICAL CONSTANTS
//  ####################################################################

   namespace constants
{
    /*
    ------------------------------------------------------------------
    Global mathematical constants
    ------------------------------------------------------------------
    */

    /// Pi constant with multiprecision accuracy
    inline float_type pi = boost::math::constants::pi<float_type>(); 
    /// Imaginary unit with multiprecision accuracy
    inline complex_type I(0,1);
    
    /*
    ------------------------------------------------------------------
    Parameters for SpinWeightedSpheroidalHarmonics
    ------------------------------------------------------------------
    */
    
    /// Maximum number of iterations for Jacobi eigenvalue solver
    inline int max_iterations_jacobi_eig {1000}; 

    /// Tolerance for convergence of Jacobi eigenvalue computation
    inline float_type tolerance_jacobi_eig{float_type(1)/1e20};

    /// Half-width of the interval used to bracket the eigenvalue root around the initial guess
    inline float_type delta_guess{float_type(1)/2};

    /// Maximum iterations in Leaver's method
    inline int i_max_leaver{200};

    /*
    ------------------------------------------------------------------
    Parameters for continued fraction computations
    Used in:
      - SpinWeightedSpheroidalHarmonics
      - RadialHomogeneousSolution
    ------------------------------------------------------------------
    */

    /// Maximum number of terms in continued fraction expansions
    inline int maxterms_continued_fraction{1000000};
    
    /// Target accuracy for continued fraction evaluation
    inline float_type accuracy_continued_fraction_computation{float_type(1)/1e40};

    /// Maximum iterations allowed in Newton root-finding algorithm
    inline int max_interations_newton_algorithm {1000000};

    /// Number of digits of precision for Newton root-finding algorithm
    inline int digits_precision_newton_algorithm = 30;
    
    /*
    ------------------------------------------------------------------
    Parameters used in homogeneous solutions and special functions
    ------------------------------------------------------------------
    */

    /// Frequency threshold where Newton method is used to find the renormalized angular momentum root instead of the monodromy method
    inline float_type frequency_threshold_for_switching_to_newton_algorithm{float_type(4)/1};
    
    /// Maximum summation index in monodromy method
    inline int nmax_monodromy = 300;

    /// Maximum number of coefficients computed in the minimal-solution series
    /// expansion (positive and negative n) used in the radial homogeneous solution
    inline int imax_coefficients = 500; 

    /// Accuracy threshold for series summations
    inline float_type accuracy_sum_computation{float_type(1)/1e30};

    /// Required number of consecutive converged terms
    inline int cons_counts = 5;

    /// Step size used in numerical integration method
    inline float_type step_size_of_numerical_integration{float_type(1)/1e8};

    /// Integration error tolerance used in numerical integration method
    inline float_type error_tolerance_integration_of_numerical_integration{float_type(1)/1e25}; 
    
    /*
    ------------------------------------------------------------------
    Parameters for ScalarFluxMode
    ------------------------------------------------------------------
    */    
    
    /// Precision used in trapezoidal quadrature for flux integrals
    inline float_type precision_trapezoidal_quadrature{float_type(1)/1e20};

    /*
    ------------------------------------------------------------------
    Parameters for GeodesicOrbitalMotion
    ------------------------------------------------------------------
    */ 
    /// Numerical tolerance used when comparing angle variables (e.g. qr ≈ π)
    /// to avoid floating-point roundoff issues
    
    inline float_type precision_q{float_type(1)/1e15};
    
   }


//  ####################################################################
//                    HELPER POWER FUNCTION WRAPPERS
//  ####################################################################

/*
  These wrappers allow consistent usage of Boost multiprecision
  exponentiation with different combinations of integer and
  multiprecision arguments.
*/

inline float_type my_pow(float_type base, float_type exp) {
    return mp::pow(base, exp);
}
inline float_type my_pow(int base, int exp) {
    return mp::pow(float_type(base), float_type(exp));
}
inline float_type my_pow(int base, float_type exp) {
    return mp::pow(float_type(base), exp);
}
inline float_type my_pow(float_type base, int exp) {
    return mp::pow(base, float_type(exp));
}



//  ####################################################################
//                OPERATOR OVERLOADS FOR MULTIPRECISION TYPES
//  ####################################################################

/*
  These operator overloads allow arithmetic operations between
  multiprecision numbers and standard C++ numeric types (int, double),
  improving code readability and avoiding explicit type conversions.
*/


/* ---- complex_type with integers ---- */
complex_type operator+(const complex_type& lhs, int rhs);
complex_type operator+(int lhs, const complex_type& rhs);
complex_type operator-(const complex_type& lhs, int rhs);
complex_type operator-(int lhs, const complex_type& rhs);
complex_type operator*(const complex_type& lhs, int rhs);
complex_type operator*(int lhs, const complex_type& rhs);
complex_type operator/(const complex_type& lhs, int rhs);
complex_type operator/(int lhs, const complex_type& rhs);

/* ---- complex_type with double ---- */
complex_type operator+(const complex_type& lhs, double rhs);
complex_type operator+(double lhs, const complex_type& rhs);
complex_type operator-(const complex_type& lhs, double rhs);
complex_type operator-(double lhs, const complex_type& rhs);
complex_type operator*(const complex_type& lhs, double rhs);
complex_type operator*(double lhs, const complex_type& rhs);
complex_type operator/(const complex_type& lhs, double rhs);
complex_type operator/(double lhs, const complex_type& rhs);


/* ---- float_type with integers ---- */
float_type operator+(const float_type& lhs, int rhs);
float_type operator+(int lhs, const float_type& rhs);
float_type operator-(const float_type& lhs, int rhs);
float_type operator-(int lhs, const float_type& rhs);
float_type operator*(const float_type& lhs, int rhs);
float_type operator*(int lhs, const float_type& rhs);
float_type operator/(const float_type& lhs, int rhs);
float_type operator/(int lhs, const float_type& rhs);

/* ---- float_type with double ---- */
float_type operator*(const float_type& lhs, double rhs);
float_type operator*(double lhs, const float_type& rhs);
float_type operator+(const float_type& lhs, double rhs);
float_type operator+(double lhs, const float_type& rhs);
float_type operator-(const float_type& lhs, double rhs);
float_type operator-(double lhs, const float_type& rhs);
float_type operator/(const float_type& lhs, double rhs);
float_type operator/(double lhs, const float_type& rhs);


#endif


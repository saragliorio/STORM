#pragma once

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#ifndef STORM_DECIMAL_DIGITS
#define STORM_DECIMAL_DIGITS 100
#endif
static_assert(STORM_DECIMAL_DIGITS >= 50, "STORM requires at least 50 decimal digits");
typedef boost::multiprecision::cpp_complex<STORM_DECIMAL_DIGITS, boost::multiprecision::backends::digit_base_10 > complex_type;
typedef boost::multiprecision::number < boost::multiprecision::cpp_bin_float <STORM_DECIMAL_DIGITS, boost::multiprecision::backends::digit_base_10 >> float_type;
namespace mp = boost::multiprecision;
namespace bm = boost::math;

#ifndef CONSTANTS_NS
#define CONSTANTS_NS

   namespace constants
{
    inline float_type pi = boost::math::constants::pi<float_type>();
    inline complex_type I(0,1);
    inline int max_iterations_jacobi_eig {1000};
    inline float_type tolerance_jacobi_eig{float_type(1)/1e20};
    inline float_type delta_guess{float_type(1)/2};
    inline int i_max_leaver{200};
    inline int maxterms_continued_fraction{1000000};
    inline float_type accuracy_continued_fraction_computation{float_type(1)/1e40};
    inline int max_interations_newton_algorithm {1000000};
    inline int digits_precision_newton_algorithm = 30;
    inline float_type frequency_threshold_for_switching_to_newton_algorithm{float_type(4)/1};
    inline int nmax_monodromy = 300;
    inline int imax_coefficients = 500;
    inline float_type accuracy_sum_computation{float_type(1)/1e30};
    inline int cons_counts = 5;
    inline float_type step_size_of_numerical_integration{float_type(1)/1e8};
    inline float_type error_tolerance_integration_of_numerical_integration{float_type(1)/1e25};
    inline float_type precision_threshold{float_type(1)/1e15};
    inline float_type precision_trapezoidal_quadrature{float_type(1)/1e20};

    inline float_type precision_q{float_type(1)/1e15};

   }

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

complex_type operator+(const complex_type& lhs, int rhs);
complex_type operator+(int lhs, const complex_type& rhs);
complex_type operator-(const complex_type& lhs, int rhs);
complex_type operator-(int lhs, const complex_type& rhs);
complex_type operator*(const complex_type& lhs, int rhs);
complex_type operator*(int lhs, const complex_type& rhs);
complex_type operator/(const complex_type& lhs, int rhs);
complex_type operator/(int lhs, const complex_type& rhs);

complex_type operator+(const complex_type& lhs, double rhs);
complex_type operator+(double lhs, const complex_type& rhs);
complex_type operator-(const complex_type& lhs, double rhs);
complex_type operator-(double lhs, const complex_type& rhs);
complex_type operator*(const complex_type& lhs, double rhs);
complex_type operator*(double lhs, const complex_type& rhs);
complex_type operator/(const complex_type& lhs, double rhs);
complex_type operator/(double lhs, const complex_type& rhs);

float_type operator+(const float_type& lhs, int rhs);
float_type operator+(int lhs, const float_type& rhs);
float_type operator-(const float_type& lhs, int rhs);
float_type operator-(int lhs, const float_type& rhs);
float_type operator*(const float_type& lhs, int rhs);
float_type operator*(int lhs, const float_type& rhs);
float_type operator/(const float_type& lhs, int rhs);
float_type operator/(int lhs, const float_type& rhs);

float_type operator*(const float_type& lhs, double rhs);
float_type operator*(double lhs, const float_type& rhs);
float_type operator+(const float_type& lhs, double rhs);
float_type operator+(double lhs, const float_type& rhs);
float_type operator-(const float_type& lhs, double rhs);
float_type operator-(double lhs, const float_type& rhs);
float_type operator/(const float_type& lhs, double rhs);
float_type operator/(double lhs, const float_type& rhs);

#endif

#include "special_functions.h"
#include <boost/math/special_functions/bernoulli.hpp>
#include <limits>
#include <stdexcept>
#include <vector>

void special_functions::sf_compute_gamma_function_z (const complex_type& z, complex_type& gamma)
{
  if(z.imag() == 0)
  {
    gamma = bm::tgamma(z.real());
    return;
  }
  if(z.real() < 0)
  {
    const float_type integer_part = mp::floor(z.real() + float_type("0.5"));
    complex_type sine = mp::sin(constants::pi * complex_type(z.real() - integer_part, z.imag()));
    if(mp::fmod(integer_part, float_type(2)) != 0)
      sine = -sine;
    sf_compute_gamma_function_z(1 - z, gamma);
    gamma = constants::pi / (sine * gamma);
    return;
  }

  static const std::vector<float_type> coefficients = [] {
    std::vector<float_type> result(std::numeric_limits<float_type>::digits10 + 1);
    for(unsigned k = 1; k < result.size(); ++k)
      result[k] = bm::bernoulli_b2n<float_type>(k) / (float_type(2 * k) * float_type(2 * k - 1));
    return result;
  }();
  complex_type shifted = z, product = 1;
  while(shifted.real() < std::numeric_limits<float_type>::digits10)
  {
    product *= shifted;
    shifted += 1;
  }
  complex_type logarithm = (shifted - float_type("0.5")) * mp::log(shifted) - shifted +
                          mp::log(2 * constants::pi) / 2;
  const complex_type inverse = 1 / shifted;
  complex_type power = inverse;
  float_type bound_power = 1 / shifted.real();
  for(unsigned k = 1; k < coefficients.size(); ++k)
  {
    // DLMF 5.11(ii): the omitted log-Gamma remainder is bounded using Re(shifted).
    if(mp::abs(coefficients[k]) * bound_power < std::numeric_limits<float_type>::epsilon())
    {
      gamma = mp::exp(logarithm) / product;
      return;
    }
    logarithm += coefficients[k] * power;
    power *= inverse * inverse;
    bound_power /= shifted.real() * shifted.real();
  }
  throw std::runtime_error("Complex Gamma expansion did not converge");
}

void special_functions::sf_compute_gamma_ratio (const complex_type& a, const int n, complex_type& gamma_ratio)
{
  complex_type product = 1;
  if(n >= 0)
  {
    for(int i = 0; i < n; ++i)
      product *= a + i;
  }
  else
  {
    for(int i = -1; i >= n; --i)
      product *= a + i;
    product = 1 / product;
  }
  gamma_ratio = product;
}

void special_functions::sf_compute_gamma_ratio(const float_type& a, const int n, float_type& gamma_ratio)
{    
    gamma_ratio = 1.;
    if(n == 0) 
    {
      gamma_ratio = 1.;
    }
    else
    {
      gamma_ratio = 1.;
      for(int i = 0; i<n; ++i)
      {
        gamma_ratio = gamma_ratio*(a + i);
      }
    }
}


//!!!converges only for -1 < z < 1
void special_functions::sf_compute_2F1_gauss(const complex_type& a , const complex_type& b , const complex_type& c , const float_type& z, complex_type& result)
{
  result = 1.0; 
  complex_type term = 1.0;
  int n = 0;
  float_type accuracy = 10000.0;  
  const int max_iterations = 1000;  

  while(accuracy > constants::accuracy_sum_computation && n < max_iterations)
  {
      complex_type ratio = ((a + n) * (b + n) * z) / ((c + n) * (n + 1));
      term *= ratio;
      result += term;
      accuracy = mp::abs(term);
      n++;
  }
}

void special_functions::sf_compute_2F1_a_or_b_negative_integer_numbers(const float_type& a , const float_type& b , const float_type& c , const float_type& z, complex_type& result)
{

  result = 1.0;  
  complex_type term = 1.0; 
  int n = 0;

  int n_max;
  if(a<0)
  {
    n_max = -static_cast<int>(a);
  }
  else 
  {
    n_max = -static_cast<int>(b);
  }

  while(n <= n_max)
  {
      complex_type ratio = ((a + n) * (b + n) * z) / ((c + n) * (n + 1));
      term *= ratio;
      result += term;
      n++;
  }
}

void special_functions::sf_compute_2F1_a_or_b_negative_integer_numbers(const complex_type& a , const complex_type& b , const complex_type& c , const float_type& z, complex_type& result)
{
  result = 1.0;  
  complex_type term = 1.0; 
  int n = 0;

  int n_max;

  if(a.real()<0)
  {
    n_max = -static_cast<int>(a.real());
  }
  else 
  {
    n_max = -static_cast<int>(b.real());
  }

  while(n <= n_max)
  {
      complex_type ratio = ((a + n) * (b + n) * z) / ((c + n) * (n + 1));
      term *= ratio;
      result += term;
      n++;
  }
}

void special_functions::sf_compute_2F1_gauss(const float_type& a , const float_type& b , const float_type& c , const float_type& z, complex_type& result)
{
    complex_type a_complex = static_cast<complex_type>(a);
    complex_type b_complex = static_cast<complex_type>(b);
    complex_type c_complex = static_cast<complex_type>(c);

    sf_compute_2F1_gauss(a_complex, b_complex, c_complex, z, result);
}

void special_functions::sf_compute_2F1_generalized (const complex_type& a , complex_type& b , const complex_type& c , const float_type& z, complex_type& result)
{
    if( a.imag() - b.imag() == 0 && (a.real()-b.real()) == (int)(a.real() - b.real()))
    {
      b = b + 0.000001+ 0.*constants::I;
    }

    complex_type gamma_c;
    sf_compute_gamma_function_z(c, gamma_c);

    complex_type gamma_b_minus_a;
    sf_compute_gamma_function_z(b-a, gamma_b_minus_a);

    complex_type gamma_b;
    sf_compute_gamma_function_z(b, gamma_b);

    complex_type gamma_c_minus_a;
    sf_compute_gamma_function_z(c-a, gamma_c_minus_a);

    complex_type gamma_a_minus_b;
    sf_compute_gamma_function_z(a-b,  gamma_a_minus_b);

    complex_type gamma_c_minus_b;
    sf_compute_gamma_function_z(c-b, gamma_c_minus_b);

    complex_type gamma_a;
    sf_compute_gamma_function_z(a, gamma_a);

    complex_type hyper_geom_func1;
    sf_compute_2F1_gauss(a, a-c+1, a-b+1, 1./z, hyper_geom_func1);

    complex_type term1 =   mp::pow(- z , - a) * gamma_b_minus_a * gamma_c*hyper_geom_func1/gamma_b/gamma_c_minus_a;

    complex_type hyper_geom_func2;
    sf_compute_2F1_gauss(b, b-c+1, - a + b + 1, 1./z, hyper_geom_func2 );
    complex_type term2 =  mp::pow( - z , - b) * gamma_a_minus_b * gamma_c*hyper_geom_func2/gamma_a/gamma_c_minus_b;

    result = term1 + term2;
}

void special_functions::sf_compute_2F1_generalized (const float_type& a , const float_type& b , const float_type& c , const float_type& z, complex_type& result)
{
    complex_type a_complex = static_cast<complex_type>(a);
    complex_type b_complex = static_cast<complex_type>(b);
    complex_type c_complex = static_cast<complex_type>(c);

    sf_compute_2F1_generalized(a_complex, b_complex, c_complex, z, result);
  
}

void special_functions::sf_compute_2F1(const complex_type& a , complex_type& b , complex_type& c , const float_type& z, complex_type& result)
{
  if(a.real() == 0 && a.imag() == 0)
  {
    result = 1.;
  }
  else if(b.real() == 0 && b.imag() == 0)
  {
    result = 1.;
  }
  else if (a == b && a == c)
  {
    result = 1./mp::pow(1-z, a);
  } 
  else if ((a.imag() == 0 && a.real() == (int)a.real() && a.real() < 0 ) || ( b.imag() == 0 &&  b.real() == (int)b.real() && b.real() < 0) ) 
  {
    sf_compute_2F1_a_or_b_negative_integer_numbers(a , b , c , z, result);
  }
  else
  { 
    if(mp::abs(z) < 1)
    {
      sf_compute_2F1_gauss(a , b , c , z, result);
    }
    else
    {      
      sf_compute_2F1_generalized(a , b , c , z, result);
    }
  }
}

void special_functions::sf_compute_2F1(const float_type& a , float_type& b , float_type& c , const float_type& z, complex_type& result)
{
  if(a == 0)
  {
    result = 1.;
  }
  else if(b == 0)
  {
    result = 1.;
  }
  else if (a == b && a == c)
  {
    result = 1./mp::pow(1-z, a);
  }
  else if ((a == (int)a && a < 0 ) || ( b == (int)b && b < 0)) 
  {
    sf_compute_2F1_a_or_b_negative_integer_numbers(a , b , c , z, result);
  }
  else
  {
    if(mp::abs(z) < 1)
    {
      sf_compute_2F1_gauss(a , b , c , z, result);
    }
    else
    {
      sf_compute_2F1_generalized(a , b , c , z, result);
    }
  }
}


void special_functions::sf_compute_1F1_gauss(const complex_type& a , const complex_type& b , const complex_type& z, complex_type& result)
{
    result = 1.0;  
    complex_type term = 1.0; 
    int n = 0;
    float_type accuracy = 10000.0;  
    const int max_iterations = 1000;  

    while(accuracy > constants::accuracy_sum_computation && n < max_iterations)
    {
        complex_type ratio = ((a + n) * z) / ((b + n) * (n + 1));
        term *= ratio;
        result += term;
        accuracy = mp::abs(term);
        n++;
    }
}

void special_functions::sf_compute_confluent_hypergeometric_function(const complex_type& a , const complex_type& b ,const complex_type& z, complex_type& U)
{
    complex_type gamma_b_minus_1;
    sf_compute_gamma_function_z(b-1, gamma_b_minus_1);
    complex_type gamma_a;
    sf_compute_gamma_function_z(a, gamma_a);
    complex_type gamma_1_minus_b;
    sf_compute_gamma_function_z(1-b, gamma_1_minus_b);
    complex_type gamma_a_minus_b_plus_1;
    sf_compute_gamma_function_z(a-b+1, gamma_a_minus_b_plus_1);

    complex_type gauss1;
    complex_type gauss2;
    sf_compute_1F1_gauss(a-b+1, 2-b, z, gauss1);
    sf_compute_1F1_gauss(a, b, z, gauss2);

    U =  gamma_b_minus_1*mp::pow(z, 1-b)*gauss1/gamma_a + gamma_1_minus_b*gauss2/gamma_a_minus_b_plus_1;
}
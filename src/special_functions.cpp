#include "special_functions.h"

void special_functions::sf_compute_gamma_function_z (const complex_type& z, complex_type& gamma)
{
  if(mp::abs(z.imag()) != 0)
  {
    float_type g("22.618910");
    
    int n = 21;

    std::vector<complex_type> coeff;
    float_type coeff0("+2.0240434640140357514731512432760");
    coeff0 *= float_type{"1e-10"};
    float_type coeff1{"+1.5333183020199267370932516012553"};
    float_type coeff2{"-1.1640274608858812982567477805332"};
    coeff2 *= float_type{"1e1"}; 
    float_type coeff3{"+4.0053698000222503376927701573076"};
    coeff3 *= float_type{"1e1"};
    float_type coeff4{"-8.2667863469173479039227422723581"};
    coeff4 *= float_type{"1e1"};
    float_type coeff5{"+1.1414465885256804336106748692495"};
    coeff5 *= float_type{"1e2"};
    float_type coeff6{"-1.1135645608449754488425056563075"};
    coeff6 *= float_type{"1e2"};
    float_type coeff7{"+7.9037451549298877731413453151252"};
    coeff7 *= float_type{"1e1"}; 
    float_type coeff8{"-4.1415428804507353801947558814560"};
    coeff8 *=  float_type{"1e1"};
    float_type coeff9{"+1.6094742170165161102085734210327"};
    coeff9 *=  float_type{"1e1"};
    float_type coeff10{"-4.6223809979028638614212851576524"};  
    float_type coeff11{"+9.7030884294357827423006360746167"};
    coeff11 *= float_type{"1e-1"};
    float_type coeff12{"-1.4607332380456449418243363858893"};
    coeff12 *= float_type{"1e-1"};
    float_type coeff13{"+1.5330325530769204955496334450658"};
    coeff13 *=  float_type{"1e-2"};
    float_type coeff14{"-1.0773862404547660506042948153734"};
    coeff14 *=   float_type{"1e-3"};
    float_type coeff15{"+4.7911128916072940196391032755132"};
    coeff15 *=  float_type{"1e-5"};
    float_type coeff16{"-1.2437781042887028450811158692678"};
    coeff16 *= float_type{"1e-6"};
    float_type coeff17{"+1.6751019107496606112103160490729"};
    coeff17 *= float_type{"1e-8"};
    float_type coeff18{"-9.7674656970897286097939311684868"};
    coeff18 *=  float_type{"1e-11"};
    float_type coeff19{"+1.8326577220560509759575892664132"};
    coeff19 *= float_type{"1e-13"};
    float_type coeff20{"-6.4508377189118502115673823719605"};
    coeff20 *=  float_type{"1e-17"};
    float_type coeff21{"+1.3382662604773700632782310392171"};
    coeff21 *=  float_type{"1e-21"};

    coeff.push_back(coeff0);
    coeff.push_back(coeff1);
    coeff.push_back(coeff2);
    coeff.push_back(coeff3);
    coeff.push_back(coeff4);
    coeff.push_back(coeff5);
    coeff.push_back(coeff6);
    coeff.push_back(coeff7);
    coeff.push_back(coeff8);
    coeff.push_back(coeff9);
    coeff.push_back(coeff10);
    coeff.push_back(coeff11);
    coeff.push_back(coeff12);
    coeff.push_back(coeff13);
    coeff.push_back(coeff14);
    coeff.push_back(coeff15);
    coeff.push_back(coeff16);
    coeff.push_back(coeff17);
    coeff.push_back(coeff18);
    coeff.push_back(coeff19);
    coeff.push_back(coeff20);
    coeff.push_back(coeff21);

    complex_type result;

    if(z.real()>= 0)
    { 
        complex_type L_g = coeff[0];
        for(int k=1; k <= n; ++k)
        {
            L_g = L_g + coeff[k] / ( z - 1. + k );
        }

        result = 2. * mp::sqrt( mp::exp(float_type("1"))/constants::pi ) * mp::pow( (z + g - 0.5)/(mp::exp(float_type("1"))) , z - 0.5 ) * L_g;
    }
    
    else
    {
        complex_type L_g = coeff[0];
        complex_type t = -z;
        for(int k=1; k <= n; ++k)
        {
            L_g = L_g + coeff[k] / ( t + k );
        }
        complex_type gamma_t_plus_1 = 2. * mp::sqrt( mp::exp(float_type("1"))/constants::pi ) * mp::pow( (t + g + 0.5)/(mp::exp(float_type("1"))) , t + 0.5 ) * L_g;
        result = constants::pi / mp::sin(constants::pi * ( 1. + t ) ) / gamma_t_plus_1;
    }
    gamma = result;
  }
  else
  {
    gamma = bm::tgamma(z.real());
  }

  // float_type r = 1000;
  // complex_type sum = 0;
  // int nmax = 999; //r-1
  // for(int n = 1; n<=nmax; ++n)
  // {
  //   float_type fact = boost::math::factorial<float_type>(n-1);
  //   float_type coeff = (std::pow(-1,1 + n)*mp::exp(-n + r)*mp::pow(-n + r,-0.5 + n))/fact;
  //   sum = sum + coeff/(n+z-1);
  // }

  // gamma = mp::exp(-r - z + 1 )*mp::pow(r + z - 1, z - 0.5)*(mp::sqrt(2*constants::pi) + sum);
  
}

void special_functions::sf_compute_gamma_ratio (const complex_type& a, const int n, complex_type& gamma_ratio)
{   
  if(mp::abs(a.imag()) == 0 && n>=0)
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
  else if(mp::abs(a.imag()) == 0 && n<=0)
  {  
    gamma_ratio = 1.;
    if(n == 0) 
    {
      gamma_ratio = 1.;
    }
    else
    {
      gamma_ratio = 1.;
      for(int i = -1; i >= n ; --i)
      {
        gamma_ratio = gamma_ratio*(a + i);
      }
      gamma_ratio = 1./gamma_ratio;
    }
  }
  else
  {
    complex_type num;
    complex_type den;
    sf_compute_gamma_function_z(a+n, num);
    sf_compute_gamma_function_z(a, den);
    gamma_ratio = num/den;
  }
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
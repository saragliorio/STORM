#include "SpinWeightedSpheroidalHarmonics.h"

SpinWeightedSpheroidalHarmonics::SpinWeightedSpheroidalHarmonics(int ss, int ll, int mm, float_type obliteness) 
    : swshv_s(ss), 
      swshv_l(ll),
      swshv_m(mm),
      swshv_zeta(obliteness)
  {    
    if (std::abs(swshv_m) > swshv_l) 
    {
        std::cerr << "Error in SpinWeightedSpheroidalHarmonics class: |mode index m must satisfy |m| <=l. " << std::endl;
        std::exit(EXIT_FAILURE);
    }

    swshv_alpha =abs(swshv_m + swshv_s);
    swshv_beta = abs(swshv_m - swshv_s);
    swshv_k1 = swshv_beta/2.;
    swshv_k2 = swshv_alpha/2.;


    if(obliteness == 0)
    {
      swshv_lambda = swshv_l*(swshv_l+1) - swshv_s*(1 + swshv_s);
    }
    else
    {
      swshm_find_root_for_eigenvalue();
    }

    shwm_compute_a_n_vector_leaver();
    shwm_compute_normalization_vector_leaver();
    shwm_compute_sign_leaver();
  }

struct SpinWeightedSpheroidalHarmonics::swhs_continued_fraction_Rn 
  {
  private:
    int k;
    SpinWeightedSpheroidalHarmonics& parent; 
  public:
    const float_type& eigenvalue;
    typedef std::pair<float_type, float_type> result_type;
    swhs_continued_fraction_Rn(int n, SpinWeightedSpheroidalHarmonics& p, const float_type& eig_value) : k(n - 1), parent(p), eigenvalue(eig_value)
    {
    }
    
    result_type operator()()
    {
      ++k;
      float_type alpha_k;
      parent.swshm_compute_alpha_n(k, alpha_k);
      float_type beta_k_plus_1;
      parent.swshm_compute_beta_n(k + 1, eigenvalue, beta_k_plus_1);
      float_type gamma_k_plus_1;
      parent.swshm_compute_gamma_n(k + 1, gamma_k_plus_1);
      return result_type(-alpha_k * gamma_k_plus_1, beta_k_plus_1);
    }
  };
  
void SpinWeightedSpheroidalHarmonics::kHat(const float_type& l, float_type& result) 
{
    if (l == 0.0 && swshv_m == 0.0) 
    {
        result = swshv_zeta * swshv_zeta / 3.0;
    }

    if (l == 0.5 && std::abs(swshv_s) == 0.5) 
    {
      result = (-648.0 + 576.0*swshv_m*swshv_zeta + 216.0*swshv_zeta*swshv_zeta + 288.0*swshv_m*swshv_m*swshv_zeta*swshv_zeta) / 864.0;
    }

    float_type term1 = -l * (1.0 + l);

    float_type term2 = (2.0 * swshv_m * swshv_s * swshv_s * swshv_zeta) / (l + l*l);

    float_type term3 = (1.0/3.0) *(1.0 + (2.0 * (l + l*l - 3.0*swshv_m*swshv_m) * (l + l*l - 3.0*swshv_s*swshv_s)) /(l * (-3.0 + l + 8.0*l*l + 4.0*l*l*l))) *swshv_zeta * swshv_zeta;

    result = term1 + term2 + term3;
}

void SpinWeightedSpheroidalHarmonics::k2(const float_type& l, float_type& result) 
{
  float_type numerator =  (1.0 + l - swshv_m) * (2.0 + l - swshv_m) *
                          (1.0 + l + swshv_m) * (2.0 + l + swshv_m) *
                          (1.0 + l - swshv_s) * (2.0 + l - swshv_s) *
                          (1.0 + l + swshv_s) * (2.0 + l + swshv_s);

  float_type denominator = (1.0 + 2.0*l) * (5.0 + 2.0*l);

  float_type sqrt_part = mp::sqrt(numerator / denominator);

  result = sqrt_part * swshv_zeta * swshv_zeta /((1.0 + l) * (2.0 + l) * (3.0 + 2.0*l));
}

void SpinWeightedSpheroidalHarmonics::kTilde2(const float_type& l, float_type& result) 
{
    if (l == 0.0 && swshv_m == 0.0) 
    {
      result = -2.0 * swshv_s * mp::sqrt((float_type(1.0) - swshv_s*swshv_s) / 3.0) * swshv_zeta;
    }
    else
    {
      float_type sqrt_part = mp::sqrt(((1.0 + 2.0*l + l*l - swshv_m*swshv_m) *(1.0 + 2.0*l + l*l - swshv_s*swshv_s)) /(3.0 + 8.0*l + 4.0*l*l));

      float_type numerator = -2.0 * swshv_s * sqrt_part * swshv_zeta * (2.0*l + l*l + swshv_m*swshv_zeta);
      float_type denominator = l * (2.0 + 3.0*l + l*l);

      result = numerator / denominator;
  }
}

void SpinWeightedSpheroidalHarmonics::swshm_find_off_diagonal_max(const std::vector<std::vector<float_type>>& A, int& index_row, int& index_column, float_type& maxVal) 
{
  int n = A.size();
  int p = 0, q = 1;
  maxVal = 0.0;
  
  for (int i = 0; i < n; i++) 
  {
      for (int j = i + 1; j < n; j++) 
      {
          if (fabs(A[i][j]) > maxVal) 
          {
              maxVal = fabs(A[i][j]);
              p = i;
              q = j;
          }
      }
  }

  index_row = p;
  index_column = q; 

}

void SpinWeightedSpheroidalHarmonics::swshm_construct_identity_matrix(std::vector<std::vector<float_type>>& J) 
{
  int size = J.size();
  for (int i = 0; i < size; ++i)
  {
      for (int j = 0; j < size; ++j) 
      {
        if(i == j)
        {
          J[i][j] = 1.;
        }
        else
        {
          J[i][j] = 0.;
        }
      }
  }
}

void SpinWeightedSpheroidalHarmonics::swshm_multiply_matrices(const std::vector<std::vector<float_type>>& A, const std::vector<std::vector<float_type>>& B, std::vector<std::vector<float_type>>& result) 
{
  int size = A.size();
  for (int i = 0; i < size; ++i) 
  {
      for (int j = 0; j < size; ++j) 
      {
          result[i][j] = 0;
          for (int k = 0; k < size; ++k) 
          {
              result[i][j] += A[i][k] * B[k][j];
          }
      }
  }
}

void SpinWeightedSpheroidalHarmonics::swshm_transpose_matrix(const std::vector<std::vector<float_type>>& A, std::vector<std::vector<float_type>>& result) 
{
  int size = A.size();
  for (int i = 0; i < size; ++i) 
  {
      for (int j = 0; j < size; ++j) 
      {
          result[j][i] = A[i][j];
      }
  }
}

void SpinWeightedSpheroidalHarmonics::swshm_jacobi_eigenvalues(std::vector<std::vector<float_type>>& A, std::vector<float_type>& eigenvalues) 
{
  float_type tol = constants::tolerance_jacobi_eig;
  int maxIter = constants::max_iterations_jacobi_eig; 
  int n = A.size();
  int iter = 0;
  
  std::vector<std::vector<float_type>> Anew = A;

  while (iter < maxIter) 
  {
      int i, j;
      float_type max_val;

      swshm_find_off_diagonal_max(A, i , j, max_val);
      
      if (max_val < tol) break;
      
      float_type theta;
      if(A[i][i] == A[j][j]) 
      {
        theta = M_PI / 4;
      } 
      else 
      {
        theta = 0.5 * mp::atan(2. * A[i][j]/ (A[i][i] - A[j][j]));
      }    
            
      float_type c = mp::cos(theta);
      float_type s = mp::sin(theta);
      
      std::vector<std::vector<float_type>> J(n, std::vector<float_type>(n, 0.0));
      swshm_construct_identity_matrix(J);
      J[i][i] = c;
      J[i][j] = -s;
      J[j][i] = s;
      J[j][j] = c;

      std::vector<std::vector<float_type>> Jt(n, std::vector<float_type>(n, 0.0));
      swshm_transpose_matrix(J, Jt);

      std::vector<std::vector<float_type>> temp(n, std::vector<float_type>(n, 0.0));
      swshm_multiply_matrices(Jt, Anew, temp);
      swshm_multiply_matrices(temp, J, A);
      
      Anew = A;
      
      iter++;
  }

  for (int i = 0; i < n; i++) 
  {
      eigenvalues.push_back(Anew[i][i]);
  }

}

void SpinWeightedSpheroidalHarmonics::swshm_compute_eigenvalue_with_spectral_method(float_type& eigenvalue)
{
  int nUp =  std::ceil(std::abs(1.5 * static_cast<double>(swshv_zeta))) + 5; 
  int minl = std::max(std::abs(swshv_m), std::abs(swshv_s));
  int nDown = std::min (swshv_l - minl, nUp);
  int nterms = nUp + nDown + 1; 
  std::vector<std::vector<float_type>> M(nterms, std::vector<float_type>(nterms, 0.0));
  
  for (int i = 0; i < nterms; ++i) 
  {
      for (int j = 0; j < nterms; ++j) 
      {
          if (i == j) 
          {
            float_type Mii;
              kHat(swshv_l - nDown + i, Mii);

              M[i][i] = Mii;
          } 
          else if (j == i - 1) 
          {
            float_type Mij1;
              kTilde2(swshv_l - nDown + i - 1, Mij1);
              M[i][j] = Mij1;
          } 
          else if (j == i - 2) 
          {
            float_type Mij2;
              k2(swshv_l - nDown + i - 2, Mij2);

              M[i][j] = Mij2;
          }
          else if (j == i + 1) 
          {
            float_type Mij1;
              kTilde2(swshv_l - nDown + i, Mij1);
              M[i][j] = Mij1;
          } 
          else if (j == i + 2) 
          {
            float_type Mij2;
              k2(swshv_l - nDown + i, Mij2);

              M[i][j] = Mij2;
          }
      }
  }

  std::vector<float_type> eigenvalues;
  swshm_jacobi_eigenvalues(M, eigenvalues);

  std::sort(eigenvalues.begin(), eigenvalues.end(), std::greater<float_type>());

  eigenvalue = - eigenvalues[nDown];
}

void SpinWeightedSpheroidalHarmonics::swshm_compute_eigenvalue_low_frequency_expansion(float_type& eigenvalue)
{
  eigenvalue = swshv_l*(1 + swshv_l) - (2*swshv_zeta*swshv_m*my_pow(swshv_s,2))/(swshv_l + my_pow(swshv_l,2)) + 
        mp::pow(swshv_zeta,2)*(-1 - (2*(swshv_l - swshv_m)*(swshv_l + swshv_m)*my_pow(swshv_l - swshv_s,2)*my_pow(swshv_l + swshv_s,2))/
        (my_pow(swshv_l,3)*(-1 + 4*my_pow(swshv_l,2))) + (2*(1 + swshv_l - swshv_m)*(1 + swshv_l + swshv_m)*my_pow(1 + swshv_l - swshv_s,2)*my_pow(1 + swshv_l + swshv_s,2))/
        (my_pow(1 + swshv_l,3)*(3 + 4*swshv_l*(2 + swshv_l))));
}

void SpinWeightedSpheroidalHarmonics::shwm_compute_alpha_n_leaver(const int n, float_type& alpha_n)
{
  alpha_n = -2.*(n+1)*(n+2*swshv_k1+1);
}

void SpinWeightedSpheroidalHarmonics::shwm_compute_beta_n_leaver(const int n, float_type& beta_n)
{
  beta_n = (swshv_k1 + swshv_k2)*(1 + swshv_k1 + swshv_k2) + (-1 + n)*n - swshv_s*(1 + swshv_s) + 2*n*(1 + swshv_k1 + swshv_k2 - 2*swshv_zeta) - 2* swshv_m *swshv_zeta - 2*(1 + 2*swshv_k1 + swshv_s)*swshv_zeta - swshv_lambda;
}

void SpinWeightedSpheroidalHarmonics::shwm_compute_gamma_n_leaver(const int n, float_type& gamma_n)
{
  gamma_n = 2*(swshv_k1 + swshv_k2 + n + swshv_s)*swshv_zeta;
}

void SpinWeightedSpheroidalHarmonics::shwm_compute_a_n_vector_leaver()
{
  swshv_an_leaver.push_back(1.);

  float_type alpha_0;
  float_type beta_0;
  shwm_compute_alpha_n_leaver(0, alpha_0);
  shwm_compute_beta_n_leaver(0, beta_0);
  swshv_an_leaver.push_back(-beta_0/alpha_0);

  for(int i=2; i<constants::i_max_leaver; ++i)
  {
    float_type alpha_n_minus_1;
    float_type beta_n_minus_1;
    float_type gamma_n_minus_1;
    shwm_compute_alpha_n_leaver(i-1, alpha_n_minus_1);
    shwm_compute_beta_n_leaver(i-1, beta_n_minus_1);
    shwm_compute_gamma_n_leaver(i-1, gamma_n_minus_1);

    float_type a_i = (-swshv_an_leaver[i-1]*beta_n_minus_1 - swshv_an_leaver[i-2]*gamma_n_minus_1)/alpha_n_minus_1;
    swshv_an_leaver.push_back(a_i);
  }
}

void SpinWeightedSpheroidalHarmonics::shwm_compute_normalization_vector_leaver()
{
  std::vector<complex_type> normalization_vector_leaver;
  for(int i = 0; i < constants::i_max_leaver; ++i)
  {
    complex_type a = 1 + i + 2*swshv_k1;
    complex_type b = i + 2*(1 + swshv_k1 + swshv_k2);
    complex_type z = 4.*swshv_zeta;
    complex_type hyper1F1;
    special_functions::sf_compute_1F1_gauss(a , b , z, hyper1F1);
    complex_type x = i + 2*(1 + swshv_k1 + swshv_k2);
    int n = -swshv_alpha- 1;
    
    complex_type poc = 1.;
    for (int j = 1; j <= -n; ++j)
    {
      poc /= (x - j);
    }

    float_type sum = 0;
    for(int j = 0; j <= i; ++j)
    {
      sum = sum + swshv_an_leaver[j]*swshv_an_leaver[i-j];
    }

    complex_type norm_i_complex = my_pow(2,i)*poc*hyper1F1*sum;
    normalization_vector_leaver.push_back(norm_i_complex);
  }
    
    swshv_normalization_leaver = 0;
    for(int i = 0; i < constants::i_max_leaver; ++i)
    {
      swshv_normalization_leaver += normalization_vector_leaver[i];
    }    
  
    complex_type gamma;
    complex_type z = 1+2*swshv_k2;
    special_functions::sf_compute_gamma_function_z (z, gamma);

    swshv_normalization_leaver = mp::sqrt(2*constants::pi*mp::pow(2, 1+2*swshv_k1+2*swshv_k2)*mp::exp(-2*swshv_zeta)*gamma*swshv_normalization_leaver); 
    
}

void SpinWeightedSpheroidalHarmonics::shwm_compute_sign_leaver()
{
  bool is_l_odd = (swshv_l % 2 != 0); 
  bool is_l_even = !is_l_odd;
  bool is_m_plus_s_odd = ((swshv_m + swshv_s) % 2 != 0); 
  bool is_m_plus_s_even = !is_m_plus_s_odd;     
  bool is_m_minus_s_odd = ((swshv_m - swshv_s) % 2 != 0); 

  if ((is_l_odd && is_m_plus_s_even) || 
      (is_l_odd && is_m_plus_s_odd && swshv_m >= swshv_s) || 
      (is_l_even && is_m_minus_s_odd && swshv_m <= swshv_s)) 
  {
    swshv_sign_leaver =  -1;
  }
  else 
  {
    swshv_sign_leaver = 1;
  }
}

void SpinWeightedSpheroidalHarmonics::shwm_compute_S(const float_type& theta, const float_type& phi)
{ 
  float_type u = mp::cos(theta);
  float_type term_u_plus_1;
  float_type term_u_minus_1;
  
  if(swshv_k1 == 0)
  {
    term_u_plus_1 = 1;
  }
  else
  {
    term_u_plus_1 = mp::pow(1+u, swshv_k1);
  }

  if(swshv_k2 == 0)
  {
    term_u_minus_1 = 1;
  }
  else
  {
    term_u_minus_1 = mp::pow(1-u, swshv_k2);
  }

  float_type sum = 1;
  for(int i = 1; i<constants::i_max_leaver; ++i)
  {
    sum = sum + swshv_an_leaver[i]*mp::pow(1 + u, i);
  }
  
  swshv_S = swshv_sign_leaver/swshv_normalization_leaver*mp::exp(swshv_zeta*u)*term_u_plus_1*term_u_minus_1*sum*mp::exp(constants::I*phi);

}

void SpinWeightedSpheroidalHarmonics::shwm_compute_S_dS_ddS(const float_type& theta, const float_type& phi)
{ 
  float_type u = mp::cos(theta);
  float_type term_u_plus_1;
  float_type term_u_minus_1;
  
  if(swshv_k1 == 0)
  {
    term_u_plus_1 = 1;
  }
  else
  {
    term_u_plus_1 = mp::pow(1+u, swshv_k1);
  }

  if(swshv_k2 == 0)
  {
    term_u_minus_1 = 1;
  }
  else
  {
    term_u_minus_1 = mp::pow(1-u, swshv_k2);
  }

  float_type sum = 1;
  for(int i = 1; i<constants::i_max_leaver; ++i)
  {
    sum = sum + swshv_an_leaver[i]*mp::pow(1 + u, i);
  }
  swshv_S = swshv_sign_leaver/swshv_normalization_leaver*mp::exp(swshv_zeta*u)*term_u_plus_1*term_u_minus_1*sum*mp::exp(constants::I*phi);

  float_type d_term_u_plus_1;
  float_type d_term_u_minus_1;
  
  if(swshv_k1-1 == 0)
  {
    d_term_u_plus_1 = 1;
  }
  else
  {
    d_term_u_plus_1 = swshv_k1*mp::pow(1+u, swshv_k1-1);
  }

  if(swshv_k2-1 == 0)
  {
    d_term_u_minus_1 = -1;
  }
  else
  {
    d_term_u_minus_1 = -swshv_k2*mp::pow(1-u, swshv_k2-1);
  }

  float_type d_sum = 0;
  for(int i = 1; i < constants::i_max_leaver; ++i)
  {
    d_sum = d_sum + swshv_an_leaver[i]*i*mp::pow(1 + u, i-1);
  }

  float_type d_S_u = (
    swshv_zeta*mp::exp(swshv_zeta*u)*term_u_plus_1*term_u_minus_1*sum + 
    mp::exp(swshv_zeta*u)*d_term_u_plus_1*term_u_minus_1*sum +
    mp::exp(swshv_zeta*u)*term_u_plus_1*d_term_u_minus_1*sum +
    mp::exp(swshv_zeta*u)*term_u_plus_1*term_u_minus_1*d_sum
  );

  swshv_S_prime = swshv_sign_leaver/swshv_normalization_leaver*mp::exp(constants::I*phi)*(-mp::sin(theta)) * d_S_u;

  float_type dd_term_u_plus_1;
  float_type dd_term_u_minus_1;
  
  if(swshv_k1-2 == 0)
  {
    dd_term_u_plus_1 = 2;
  }
  else
  {
    dd_term_u_plus_1 = (-1+swshv_k1)*swshv_k1*mp::pow(1+u, swshv_k1-2);
  }

  if(swshv_k2-2 == 0)
  {
    dd_term_u_minus_1 = 2;
  }
  else
  {
    dd_term_u_minus_1 = (-1+swshv_k2)*swshv_k2*mp::pow(1-u, swshv_k2-2);
  }

  float_type dd_sum = 0;
  for(int i = 2; i<constants::i_max_leaver; ++i)
  {
    dd_sum = dd_sum + swshv_an_leaver[i]*i*(i-1)*mp::pow(1 + u, i-2);
  }

  float_type dd_S_u = 
  (2*(term_u_minus_1*d_term_u_plus_1 + term_u_plus_1*d_term_u_minus_1)*(sum*swshv_zeta*mp::exp(swshv_zeta*u) + mp::exp(swshv_zeta*u)*d_sum) 
  + 
  mp::exp(swshv_zeta*u) * sum * ( 2 * d_term_u_plus_1 * d_term_u_minus_1 + term_u_minus_1 * dd_term_u_plus_1 + term_u_plus_1 * dd_term_u_minus_1)
  + 
  term_u_plus_1*term_u_minus_1*(2*swshv_zeta*mp::exp(swshv_zeta*u)*d_sum +sum*swshv_zeta*swshv_zeta*mp::exp(swshv_zeta*u) + mp::exp(swshv_zeta*u)*dd_sum));
  
  swshv_S_double_prime = (dd_S_u*(mp::sin(theta)*mp::sin(theta)) + d_S_u * (-mp::cos(theta)))*swshv_sign_leaver/swshv_normalization_leaver*mp::exp(constants::I*phi);
}

void SpinWeightedSpheroidalHarmonics::swshm_compute_eigenvalue_high_frequency_expansion(float_type& eigenvalue)
{

  int slm = std::abs(swshv_m + std::abs(swshv_s)) + std::abs(swshv_s);

  int z0 = ((swshv_l + swshv_m) % 2 == 0) ? 0 : 1;

  int swshq;
  if (swshv_l >= slm)
  {
    swshq = swshv_l + 1 - z0;

  }
  else
  {
    swshq = 2*swshv_l + 1 - slm;
  }

  float_type A1 = -1.0/8.0 * (my_pow(swshq,3)
                - swshv_m*swshv_m*swshq
                + swshq
                - 2.0*swshv_s*swshv_s*(swshq + swshv_m));

  float_type A2 = 1.0/64.0 * (-my_pow(swshv_m,4)
                + 6.0*swshv_m*swshv_m*my_pow(swshq,2)
                + 2.0*swshv_m*swshv_m
                - 5.0*my_pow(swshq,4)
                - 10.0*my_pow(swshq,2)
                - 1.0
                + 4.0*swshv_s*swshv_s*(swshv_m*swshv_m + swshv_m*swshq + 3.0*my_pow(swshq,2) + 1.0));

  float_type A3 = 1.0/512.0 * (
                -swshq * (37.0 + 13.0*my_pow(swshv_m,4)
                + 114.0*my_pow(swshq,2)
                + 33.0*my_pow(swshq,4)
                - 2.0*swshv_m*swshv_m*(25.0 + 23.0*my_pow(swshq,2)))
                + 4.0*(13.0*swshv_m - my_pow(swshv_m,3)
                + 25.0*swshq
                + 9.0*swshv_m*swshv_m*swshq
                + 33.0*swshv_m*my_pow(swshq,2)
                + 23.0*my_pow(swshq,3))*swshv_s*swshv_s
                - 8.0*(swshv_m + swshq)*my_pow(swshv_s,4));

  float_type A4 = 1.0/1024.0 * (
                -14.0
                + 2.0*my_pow(swshv_m,6)
                - 239.0*my_pow(swshq,2)
                - 340.0*my_pow(swshq,4)
                - 63.0*my_pow(swshq,6)
                - 3.0*my_pow(swshv_m,4)*(6.0 + 13.0*my_pow(swshq,2))
                + 10.0*swshv_m*swshv_m*(3.0 + 23.0*my_pow(swshq,2) + 10.0*my_pow(swshq,4))
                + 4.0*(-my_pow(swshv_m,4)
                - 9.0*my_pow(swshv_m,3)*swshq
                + 5.0*swshv_m*swshv_m*(2.0 + 3.0*my_pow(swshq,2))
                + swshv_m*swshq*(93.0 + 73.0*my_pow(swshq,2))
                + 5.0*(3.0 + 23.0*my_pow(swshq,2) + 19.0*my_pow(swshq,4))*swshv_s*swshv_s
                - 8.0*(2.0 + 3.0*swshv_m*swshv_m + 9.0*swshv_m*swshq + 6.0*my_pow(swshq,2))*my_pow(swshv_s,4))
                );

  eigenvalue = -swshv_zeta*swshv_zeta
                  + 2.0*swshq*swshv_zeta
                  - (swshq*swshq - swshv_m*swshv_m - 2.0*swshv_s*swshv_s + 1.0)/2.0
                  + A1/swshv_zeta
                  + A2/(swshv_zeta*swshv_zeta)
                  + A3/my_pow(swshv_zeta,3)
                  + A4/my_pow(swshv_zeta,4);

}

void SpinWeightedSpheroidalHarmonics::swshm_compute_alpha_n(const int n, float_type& alpha_n)
{
  if(swshv_alpha == 0 && swshv_beta == 0)
  {
    alpha_n = 4.*swshv_zeta*my_pow(1+n, 3)/(2*n +2.)/(2*n +3.);
  }
  else
  {
    alpha_n=4*swshv_zeta*(n+swshv_alpha+1)*(n+swshv_beta+1)*(n+(swshv_alpha+swshv_beta)/2+1-swshv_s)/(2*n+swshv_alpha+swshv_beta+2)/(2*n+swshv_alpha+swshv_beta+3);
  }
}

void SpinWeightedSpheroidalHarmonics::swshm_compute_beta_n(const int n, const float_type& e, float_type& beta_n)
{
  if(swshv_alpha == 0 && swshv_beta == 0)
  {
    beta_n = e + mp::pow(swshv_zeta,2)-(n+1)*n;
  }
  else
  {
    beta_n = e + swshv_zeta*swshv_zeta-(n+(swshv_alpha+swshv_beta)/2)*(n+1+(swshv_alpha+swshv_beta)/2)+2*swshv_zeta*swshv_s*(swshv_alpha-swshv_beta)*(swshv_alpha+swshv_beta)/(2*n+swshv_alpha+swshv_beta)/(2*n+swshv_alpha+swshv_beta+2);
  }
}

void SpinWeightedSpheroidalHarmonics::swshm_compute_gamma_n(const int n, float_type& gamma_n)   
{
  if(swshv_alpha == 0 && swshv_beta == 0)
  {
    gamma_n = - 2*swshv_zeta*my_pow(n,2)/(2*n -1.);
  }
  else
  {
    gamma_n=-4*swshv_zeta*n*(n+swshv_alpha+swshv_beta)*(n+(swshv_alpha+swshv_beta)/2+swshv_s)/(2*n+swshv_alpha+swshv_beta-1.)/(2*n+swshv_alpha+swshv_beta);
  }
}

void SpinWeightedSpheroidalHarmonics::swshm_compute_Rn(const int n, const float_type& E, float_type& Rn)
{
    float_type gamma_n;
    swshm_compute_gamma_n (n, gamma_n);
    float_type beta_n;
    swshm_compute_beta_n (n, E, beta_n);

    swhs_continued_fraction_Rn gen_Rn(n, *this, E);
    std::uintmax_t maxterms_continued_fraction_computation = constants::maxterms_continued_fraction;

    Rn = - gamma_n / (beta_n + boost::math::tools::continued_fraction_a(gen_Rn, constants::accuracy_continued_fraction_computation , maxterms_continued_fraction_computation));
}

void SpinWeightedSpheroidalHarmonics::swshm_compute_Ln (const int n,  const float_type& eigenvalue, float_type& Ln)
{
  int k;
  float_type alpha0, beta0, L_k_minus_1, alpha_k, beta_k, gamma_k;
  swshm_compute_alpha_n(0, alpha0);
  swshm_compute_beta_n(0, eigenvalue, beta0);
 
  Ln = -alpha0/beta0;
  k = 0;
  while(k < n)
  {
    k++;
    L_k_minus_1 = Ln;
    swshm_compute_alpha_n(k, alpha_k);
    swshm_compute_beta_n(k, eigenvalue, beta_k);
    swshm_compute_gamma_n(k, gamma_k);

    Ln = -alpha_k/(beta_k+gamma_k*L_k_minus_1);
  }
}

void SpinWeightedSpheroidalHarmonics::swshm_compute_function_eigenvalue(const float_type& e, float_type& f_e)
{
    float_type beta_n; 
    int n = (swshv_l - (swshv_alpha+swshv_beta)*0.5);
    swshm_compute_beta_n(n, e, beta_n);
    float_type alpha_n;
    swshm_compute_alpha_n(n, alpha_n);
    float_type gamma_n;
    swshm_compute_gamma_n(n, gamma_n);

    float_type R_1;
    float_type L_minus_1;
    swshm_compute_Rn(n+1, e, R_1);
    swshm_compute_Ln(n-1, e, L_minus_1);

    f_e =  beta_n + alpha_n * R_1 + gamma_n * L_minus_1;
}

void SpinWeightedSpheroidalHarmonics::swshm_compute_derivative_function_eigenvalue(const float_type& e, float_type& f_prime) 
{  
  auto f_d = [this](float_type x)
  {
    float_type f_e;
    swshm_compute_function_eigenvalue(x, f_e);
    return f_e;
  };
  f_prime = boost::math::differentiation::finite_difference_derivative(f_d, e);
}

void SpinWeightedSpheroidalHarmonics::swshm_find_root_for_eigenvalue()
{
  
  // float_type guess_eigenvalue_low_freq;
  // swshm_compute_eigenvalue_low_frequency_expansion(guess_eigenvalue_low_freq);
  // std::cout << "Guess found with low frequency expansion: " << guess_eigenvalue_low_freq << std::endl;

  float_type eigenvalue_spectral_method;
  swshm_compute_eigenvalue_with_spectral_method(eigenvalue_spectral_method);
  // std::cout << "Guess found with spectral method: " << eigenvalue_spectral_method << std::endl;
  float_type guess_eigenvalue_spectral_method = eigenvalue_spectral_method;

  // float_type guess_eigenvalue_high_freq;
  // swshm_compute_eigenvalue_high_frequency_expansion(guess_eigenvalue_high_freq);
  // std::cout << "Guess found with high frequency expansion: " << guess_eigenvalue_high_freq << std::endl;

 
  float_type min = guess_eigenvalue_spectral_method - constants::delta_guess;
  float_type max = guess_eigenvalue_spectral_method + constants::delta_guess;
  std::uintmax_t maxit_spectral = constants::max_interations_newton_algorithm;  

  auto f_FT_spectral = [this](float_type e) 
  { 
    float_type f_e;
    swshm_compute_function_eigenvalue(e, f_e);
    float_type f_e_prime;
    swshm_compute_derivative_function_eigenvalue(e, f_e_prime);
    return std::make_pair(f_e, f_e_prime); 
  };
  
  int digits_precision = constants::digits_precision_newton_algorithm;  
  float_type eigenvalue_spectral = boost::math::tools::newton_raphson_iterate(f_FT_spectral, eigenvalue_spectral_method, min, max, digits_precision, maxit_spectral);
  swshv_eigenvalue = eigenvalue_spectral;
  float_type f_e;
  swshm_compute_function_eigenvalue(eigenvalue_spectral, f_e);
  float_type lambda_spectral = eigenvalue_spectral - swshv_s*(swshv_s + 1) + mp::pow(swshv_zeta, 2) - 2*swshv_m*swshv_zeta;

  swshv_lambda = lambda_spectral;
}

void SpinWeightedSpheroidalHarmonics::swshm_get_eigenvalue(float_type& eig)
{
  eig = swshv_eigenvalue;
}

void SpinWeightedSpheroidalHarmonics::swshm_get_lambda(float_type& lambda)
{
  lambda = swshv_lambda;
}

void SpinWeightedSpheroidalHarmonics::swshm_get_S(const float_type& x,const float_type& phi, complex_type& S)
{
  shwm_compute_S(mp::acos(x), phi);
  S = swshv_S;
}

void SpinWeightedSpheroidalHarmonics::swshm_get_S_dS_ddS(const float_type& x,const float_type& phi, complex_type& S,complex_type& dS, complex_type& ddS)
{
  shwm_compute_S_dS_ddS(mp::acos(x), phi);
  S = swshv_S;
  dS = swshv_S_prime;
  ddS = swshv_S_double_prime;
}

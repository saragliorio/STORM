#include "RadialHomogeneousSolution.h"

RadialHomogeneousSolution::RadialHomogeneousSolution(int s, float_type M, float_type a, float_type m, float_type lambda, float_type omega, int ll, int mm) : 
        rhsv_spin_weight_s (s), 
        rhsv_mass_black_hole_M (M), 
        rhsv_mass_test_particle_m (m), 
        rhsv_black_hole_spin_a (a),
        rhsv_angular_eingenvalue_lambda (lambda),
        rhsv_angular_frequency_omega (omega),
        rhsv_l {ll},
        rhsv_m {mm}
    {
        if (rhsv_angular_frequency_omega == 0) 
        {
            std::cerr << "Error in RadialHomogeneousSolution class: orbital angular freqeuncy paraleter omega cannot be zero (missing implementation)." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if (mp::abs(rhsv_black_hole_spin_a) > 0.999) 
        {
            std::cerr << "Error in RadialHomogeneousSolution class: black hole spin parameter a must be <= 0.999. " << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (mp::abs(rhsv_black_hole_spin_a) < -0.999) 
        {
            std::cerr << "Error in RadialHomogeneousSolution class: black hole spin parameter a must be >= -0.999. " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if (std::abs(rhsv_m) > rhsv_l) 
        {
            std::cerr << "Error in RadialHomogeneousSolution class: mode index m must satisfy |m| <=l. " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        rhsv_numerical_integration_in = 0;
        rhsv_numerical_integration_up = 0;

        if(rhsv_angular_frequency_omega < 0)
        {
            rhsv_angular_frequency_omega = - rhsv_angular_frequency_omega;
            rhsv_m = - rhsv_m;
            rhsv_negative_freq = 1;
        }
        else
        {
            rhsv_negative_freq = 0;
        }

        rhsv_black_hole_spin_per_unit_mass_q = rhsv_black_hole_spin_a/rhsv_mass_black_hole_M;
        rhsv_epsilon = 2.*rhsv_mass_black_hole_M*rhsv_angular_frequency_omega;
        rhsv_kappa = mp::pow( 1. - mp::pow(rhsv_black_hole_spin_per_unit_mass_q,2.),0.5 );
        rhsv_tau = ( rhsv_epsilon - rhsv_m*rhsv_black_hole_spin_per_unit_mass_q )/rhsv_kappa;    

        if(rhsv_angular_frequency_omega == 0)
        {
            rhsv_nu_root =  rhsv_l;
        }
        else
        {
            rhsm_monodromy_method_for_nu(); 
            
            if(mp::abs(rhsv_nu_root.imag()) < 1e-40 && rhsv_angular_frequency_omega <= constants::frequency_threshold_for_switching_to_newton_algorithm)
            {
                rhsm_newton_algorithm_complex();
            }

        }
        rhsm_compute_constants_quantities();
    }

    struct RadialHomogeneousSolution::continued_fraction_Rn 
    {
    private:
        int k;
        RadialHomogeneousSolution& parent; 
    public:
    const complex_type& nu;
        typedef std::pair<complex_type, complex_type> result_type;

        continued_fraction_Rn(int n, RadialHomogeneousSolution& p, const complex_type& ren_ang_mom)
            : k(n - 1), parent(p), nu(ren_ang_mom)
        {
        }

        result_type operator()()
        {
            ++k;
            complex_type alpha_k;
            parent.rhsm_compute_alpha_n(k, alpha_k, nu);
            complex_type beta_k_plus_1;
            parent.rhsm_compute_beta_n(k + 1, beta_k_plus_1, nu);
            complex_type gamma_k_plus_1;
            parent.rhsm_compute_gamma_n(k + 1, gamma_k_plus_1, nu);

            return result_type(-alpha_k * gamma_k_plus_1, beta_k_plus_1);
        }
    };

    struct RadialHomogeneousSolution::continued_fraction_Ln 
    {
    private:
    int k;
        RadialHomogeneousSolution& parent;
    public:
        const complex_type nu;
        typedef std::pair<complex_type, complex_type> result_type;

        continued_fraction_Ln(int n, RadialHomogeneousSolution& p, const complex_type& ren_ang_mon)
            : k(n + 1 ), parent(p), nu(ren_ang_mon)
        {
        }

        result_type operator()()
        {
            --k;
            complex_type alpha_k_minus_1;
            parent.rhsm_compute_alpha_n(k - 1, alpha_k_minus_1, nu);
            complex_type beta_k_minus_1;
            parent.rhsm_compute_beta_n(k - 1, beta_k_minus_1, nu);
            complex_type gamma_k;
            parent.rhsm_compute_gamma_n(k , gamma_k, nu);
            return result_type( - alpha_k_minus_1 * gamma_k, beta_k_minus_1);
        }
    };

    struct RadialHomogeneousSolution::numerical_integration 
    {
        private:
        RadialHomogeneousSolution& parent; 
        const int H;
        public:
        numerical_integration(RadialHomogeneousSolution& p, const int HH)
                : parent(p), H(HH)
            {
            }

            void operator()(const std::array<complex_type, 2> &x, std::array<complex_type, 2> &dxdr, const float_type r) const 
            {
            
            complex_type y2der;
            complex_type y2;
            complex_type y1;

            parent.rhsm_compute_coefficients_of_eq_numerical_integration(H, r, y2der, y2, y1);
            
            dxdr[0] = x[1];
            dxdr[1] = - y2/y2der*x[1] - y1/y2der*x[0];


            }
    };

    struct RadialHomogeneousSolution::numerical_integration_infinity
    {
        private:
        RadialHomogeneousSolution& parent; 
        const int H;
        const float_type num_inf;
        public:
        numerical_integration_infinity(RadialHomogeneousSolution& p, const int HH, const float_type& num_inf_temp)
                : parent(p), H(HH), num_inf(num_inf_temp)
            {
            }

            void operator()(const std::array<complex_type, 2> &x, std::array<complex_type, 2> &dxdr, const float_type rho) const 
            {
            
            complex_type y2der;float_type nu_real;
            complex_type y2;
            complex_type y1;

            float_type r = num_inf - rho;
            parent.rhsm_compute_coefficients_of_eq_numerical_integration(H, r, y2der, y2, y1);
            y2 = -y2;
            dxdr[0] =  x[1];
            dxdr[1] =  - y2/y2der*x[1] - y1/y2der*x[0];
            }
    };

void RadialHomogeneousSolution::rhsm_compute_alpha_n (const int n,  complex_type& alpha_n, const complex_type& nu)
{
    complex_type num = constants::I*rhsv_epsilon*rhsv_kappa* (n + nu + 1. + rhsv_spin_weight_s + constants::I*rhsv_epsilon )*( n + nu + 1. + rhsv_spin_weight_s - constants::I*rhsv_epsilon )*( n + nu + 1. + constants::I*rhsv_tau );
    complex_type den = ( n + nu + 1. )*( 2.*n + 2.*nu + 3. );
    alpha_n = num/den;
}

void RadialHomogeneousSolution::rhsm_compute_beta_n (const int n, complex_type& beta_n, const complex_type& nu)
{   
    complex_type num = rhsv_epsilon*( rhsv_epsilon - rhsv_m * rhsv_black_hole_spin_per_unit_mass_q )*( my_pow(rhsv_spin_weight_s,2.) + mp::pow(rhsv_epsilon,2.) );
    complex_type den = ( n + nu )*( n + nu + 1. );
    beta_n = - rhsv_angular_eingenvalue_lambda - rhsv_spin_weight_s*(rhsv_spin_weight_s + 1) + (n + nu)*(n + nu + 1) + mp::pow(rhsv_epsilon,2.) + rhsv_epsilon*(rhsv_epsilon - rhsv_m * rhsv_black_hole_spin_per_unit_mass_q) + num/den;
}

void RadialHomogeneousSolution::rhsm_compute_gamma_n (const int n, complex_type& gamma_n, const complex_type& nu)
{ 
    complex_type num = constants::I*rhsv_epsilon*rhsv_kappa*( n + nu - rhsv_spin_weight_s + constants::I*rhsv_epsilon )*(n + nu - rhsv_spin_weight_s - constants::I*rhsv_epsilon)*( n + nu - constants::I*rhsv_tau );
    complex_type den = ( n + nu )*( 2.*n + 2.*nu - 1. );
    gamma_n = - num/den;
}

void RadialHomogeneousSolution::rhsm_compute_Rn (const int n, complex_type& Rn, const complex_type& nu)
{
    complex_type gamma_n{0};
    rhsm_compute_gamma_n (n, gamma_n, nu);

    complex_type beta_n{0};
    rhsm_compute_beta_n (n, beta_n, nu);

    continued_fraction_Rn gen_Rn(n, *this, nu);
    std::uintmax_t maxterms_continued_fraction_computation = constants::maxterms_continued_fraction;

    Rn = - gamma_n / (beta_n + boost::math::tools::continued_fraction_a(gen_Rn, constants::accuracy_continued_fraction_computation , maxterms_continued_fraction_computation));

}

void RadialHomogeneousSolution::rhsm_compute_Ln (const int n, complex_type& Ln, const complex_type& nu)
{
    complex_type alpha_n{0};
    rhsm_compute_alpha_n (n, alpha_n, nu);

    complex_type beta_n{0};
    rhsm_compute_beta_n (n, beta_n, nu);

    continued_fraction_Ln gen_Ln(n, *this, nu);
    std::uintmax_t maxterms_continued_fraction_computation = constants::maxterms_continued_fraction;
    
    Ln = - alpha_n / (beta_n + boost::math::tools::continued_fraction_a(gen_Ln, constants::accuracy_continued_fraction_computation , maxterms_continued_fraction_computation));
    

}

void RadialHomogeneousSolution::rhsm_function_renormalized_angular_momentum(const complex_type& nu, complex_type& f_nu)
{
    complex_type beta_0{0}; 
    rhsm_compute_beta_n(0, beta_0, nu);
    
    complex_type alpha_0{0};
    rhsm_compute_alpha_n(0, alpha_0, nu);
    
    complex_type gamma_0{0};
    rhsm_compute_gamma_n(0, gamma_0, nu);
    
    complex_type R_1{0};
    rhsm_compute_Rn(1, R_1, nu);
    
    complex_type L_minus_1{0};
    rhsm_compute_Ln(-1, L_minus_1, nu);

    f_nu =  beta_0 + alpha_0 * R_1 + gamma_0 * L_minus_1;

}

void RadialHomogeneousSolution::rhsm_function_renormalized_angular_momentum_derivative(const complex_type& nu, complex_type& f_prime) 
{  
    float_type f_real_prime{0};
    auto f_r = [this, &nu](float_type x)
    {
        complex_type f_nu{0};
        complex_type temporary_nu = x + constants::I*nu.imag();
        rhsm_function_renormalized_angular_momentum(temporary_nu, f_nu);
        return f_nu.real();
    };
    f_real_prime = boost::math::differentiation::finite_difference_derivative(f_r, nu.real());

    float_type f_imag_prime{0};
    auto f_i = [this, &nu](float_type x)
    {
        complex_type f_nu{0};
        complex_type temporary_nu = x + constants::I*nu.imag();
        rhsm_function_renormalized_angular_momentum(temporary_nu, f_nu);
        return f_nu.imag();
    };
    f_imag_prime = boost::math::differentiation::finite_difference_derivative(f_i, nu.real());

    f_prime = f_real_prime + constants::I*f_imag_prime ;
}

void RadialHomogeneousSolution::rhsm_newton_algorithm_complex()
{
    float_type my_guess_real = rhsv_nu_root.real();
    float_type min = my_guess_real - 1.;
    float_type max= my_guess_real + 1.;
    std::uintmax_t maxit = constants::max_interations_newton_algorithm;  
    int digits_precision = constants::digits_precision_newton_algorithm;
    auto f = [this](complex_type nu) 
    { 
        complex_type f_nu{0};
        rhsm_function_renormalized_angular_momentum(nu, f_nu);
        complex_type f_nu_prime{0};
        rhsm_function_renormalized_angular_momentum_derivative(nu, f_nu_prime);
        return std::make_pair(f_nu.real(), f_nu_prime.real()); 
    };
    rhsv_nu_root = boost::math::tools::newton_raphson_iterate(f, my_guess_real, min, max, digits_precision, maxit);
}

void RadialHomogeneousSolution::rhsm_monodromy_method_for_nu()
{
    complex_type gammaCH = 1 - rhsv_spin_weight_s - constants::I*rhsv_epsilon - constants::I*rhsv_tau;
    complex_type deltaCH = 1 + rhsv_spin_weight_s + constants::I*rhsv_epsilon - constants::I*rhsv_tau;
    complex_type epsilonCH = 2*constants::I*rhsv_epsilon*rhsv_kappa;
    complex_type alphaCHepsilonCH = 1 - rhsv_spin_weight_s + constants::I*(rhsv_epsilon - rhsv_tau);
    complex_type qCH = rhsv_spin_weight_s*(1 + rhsv_spin_weight_s) - mp::pow(rhsv_epsilon,2.) - constants::I*(-1 + 2*rhsv_spin_weight_s)*rhsv_epsilon*rhsv_kappa + rhsv_angular_eingenvalue_lambda + rhsv_tau*(constants::I + rhsv_tau);
    complex_type mu1C = alphaCHepsilonCH - gammaCH - deltaCH;
    complex_type mu2C = - alphaCHepsilonCH;

    int nmax = constants::nmax_monodromy;

    std::vector<complex_type> a1_coefficients;
    a1_coefficients.push_back(0.); 
    a1_coefficients.push_back(1.);   

    for(int n = 1; n <= nmax + 1; ++n)
    {
        complex_type i_element =    ((1 - n + alphaCHepsilonCH - deltaCH)*(2 - n + alphaCHepsilonCH -gammaCH - deltaCH)*epsilonCH*a1_coefficients[-2 + n +1])/n - 
                                    ((my_pow(n,2.) - qCH + mp::pow(alphaCHepsilonCH,2.) + n*(-1 +gammaCH + deltaCH - epsilonCH) + epsilonCH - deltaCH*epsilonCH + alphaCHepsilonCH*(1 - 2*n -gammaCH - deltaCH + epsilonCH))*a1_coefficients[-1 + n + 1])/n;

        a1_coefficients.push_back(i_element);
    }

    std::vector<complex_type> a2_coefficients;
    a2_coefficients.push_back(0.); 
    a2_coefficients.push_back(1.);  
    for(int i = 1; i <= nmax + 1; ++i)
    {
        complex_type i_element = -(((-2 + i + alphaCHepsilonCH)*(-1 + i + alphaCHepsilonCH -gammaCH)*epsilonCH*a2_coefficients[-2 + i + 1])/i) + ((my_pow(i,2.) - qCH + mp::pow(alphaCHepsilonCH,2.) +gammaCH + deltaCH - i*(1 +gammaCH + deltaCH - epsilonCH) - epsilonCH + 
        alphaCHepsilonCH*(-1 + 2*i -gammaCH - deltaCH + epsilonCH))*a2_coefficients[-1 + i + 1])/i;
        a2_coefficients.push_back(i_element);
    }

    std::vector<complex_type> pochhammerp1m2;
    pochhammerp1m2.push_back(1); 
    for(int i = 1; i <= nmax; ++i)
    {
        complex_type i_element = (-mu2C+mu1C+i-1)*pochhammerp1m2[i-1];
        pochhammerp1m2.push_back(i_element);
    }

    std::vector<complex_type> pochhammerm1p2;

    pochhammerm1p2.push_back(1); 
    for(int i = 1; i <= nmax; ++i)
    {
        complex_type i_element = (mu2C-mu1C+i-1)*pochhammerm1p2[i-1];
        pochhammerm1p2.push_back(i_element);
    }

        
    int J_max = std::ceil(nmax*0.5);
    
    complex_type sum_a1 = 0;

    for(int j = 0; j<=J_max; ++j)
    {
      sum_a1 = sum_a1 + a1_coefficients[j+1]*pochhammerp1m2[nmax - j];
    }

    complex_type gamma;
    complex_type z_gamma = -mu2C + mu1C;
    
    complex_type sum = 0;
    float_type r = 1000;
    int nmax_gamma = 999; 
    for(int n = 1; n<=nmax_gamma; ++n)
    {
        float_type fact = boost::math::factorial<float_type>(n-1);
        float_type coeff = (my_pow(-1,1 + n)*mp::exp(-n + r)*mp::pow(-n + r,-0.5 + n))/fact;
        sum = sum + coeff/(n+z_gamma-1);
    }

    gamma = mp::exp(-r - z_gamma + 1 )*mp::pow(r + z_gamma - 1, z_gamma - 0.5)*(mp::sqrt(2*constants::pi) + sum);

    special_functions::sf_compute_gamma_function_z(-mu2C + mu1C, gamma);
    
    sum_a1 = sum_a1*gamma;
    
    complex_type sum_a2 = 0;

    for(int j = 0; j<=J_max; ++j)
    {
      sum_a2 = sum_a2 + my_pow(-1, j)*a2_coefficients[j+1]*pochhammerm1p2[nmax - j];
    }

    gamma = gamma.real() - constants::I*gamma.imag();
    
    sum_a2 = sum_a2*gamma;

    complex_type cos2pinu = (2.*my_pow(-1,-1 + nmax)*constants::pi*constants::pi*a1_coefficients[nmax+1]*a2_coefficients[nmax+1])/
    (sum_a1*sum_a2) + mp::cos(constants::pi*(mu1C - mu2C));

    complex_type Recos2pinu = cos2pinu.real();
    complex_type acos2pinu = - mp::acos(Recos2pinu);  
    complex_type z = Recos2pinu;  
        
    acos2pinu = - constants::I* mp::log(z + constants::I * mp::sqrt(1 - z*z))*constants::I;
    acos2pinu = mp::acos(z);

    if(Recos2pinu.real()<-1)
    {
        rhsv_nu_root = 0.5  - constants::I*acos2pinu.imag()/2./constants::pi;
    }
    else if(Recos2pinu.real() > 1)
    {
        rhsv_nu_root = - constants::I*acos2pinu.imag()/2./constants::pi;
    }
    else
    {
        rhsv_nu_root = rhsv_l - acos2pinu/2./constants::pi;
    }
}

void RadialHomogeneousSolution::rhsm_compute_minimal_coefficient_n_positive(const int n, const complex_type& nu, std::vector<complex_type>& rhsv_minimal_coefficient_n_positive)
{   
    int size = static_cast<int>(rhsv_minimal_coefficient_n_positive.size());
    if(size <= n)
    {
        for(int i = size; i <= n; ++i)
        {
            complex_type Rn_value;
            rhsm_compute_Rn(i, Rn_value, nu);
            complex_type last_computed_minimal_coeff_positive =rhsv_minimal_coefficient_n_positive.back();
            rhsv_minimal_coefficient_n_positive.push_back(Rn_value * last_computed_minimal_coeff_positive);
        }
    }
}

void RadialHomogeneousSolution::rhsm_compute_minimal_coefficient_n_negative(const int n, const complex_type& nu, std::vector<complex_type>& rhsv_minimal_coefficient_n_negative)
{ 
int size = static_cast<int>(rhsv_minimal_coefficient_n_negative.size());

    if(size <= n)
    {
        for(int i = size; i <= n; ++i)
        {
            complex_type Ln_value;
            rhsm_compute_Ln(-i, Ln_value, nu);
            complex_type last_computed_minimal_coeff_negative = rhsv_minimal_coefficient_n_negative.back();
            rhsv_minimal_coefficient_n_negative.push_back(Ln_value * last_computed_minimal_coeff_negative);
            }
    }
}

void RadialHomogeneousSolution::rhsm_sum_minimal_coefficients(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative , complex_type& result)
{   
    float_type accuracy {1000};
    complex_type sum {rhsv_minimal_coefficient_n_positive[0]};
    int i {1};
    
    int consecutive_counts{0};
    
    while(consecutive_counts < constants::cons_counts)
    {   
        complex_type addend_1 {0};
        complex_type addend_2 {0};
        addend_1 = rhsv_minimal_coefficient_n_positive[i];
        addend_2 = rhsv_minimal_coefficient_n_negative[i];
        sum = sum + addend_1 + addend_2;

        if(mp::abs(addend_1/sum) > mp::abs(addend_2/sum))
        {
            accuracy = mp::abs(addend_1/sum);
        }
        else
        {
            accuracy = mp::abs(addend_2/sum);
        }

        if(accuracy < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }
        i++;  
    }
    result = sum;
}

void RadialHomogeneousSolution::rhsm_compute_A_minus(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, complex_type& A_minus)
{
    int i {1};
    float_type accuracy {10};
    complex_type sum = 0;
    complex_type gamma_ratio1;
    complex_type gamma_ratio2;
    special_functions::sf_compute_gamma_ratio(nu + 1. + rhsv_spin_weight_s - constants::I*rhsv_epsilon, 0, gamma_ratio1);
    special_functions::sf_compute_gamma_ratio(nu + 1. - rhsv_spin_weight_s + constants::I*rhsv_epsilon, 0, gamma_ratio2);
    
    complex_type arg_gamma_ratio_1 = nu + 1. + rhsv_spin_weight_s - constants::I*rhsv_epsilon;
    complex_type arg_gamma_ratio_2 = nu + 1. - rhsv_spin_weight_s + constants::I*rhsv_epsilon;

    complex_type gamma_ratio3;
    complex_type gamma_ratio4;

    complex_type gamma_ratio1_plus = gamma_ratio1;
    complex_type gamma_ratio2_plus = gamma_ratio2;

    complex_type gamma_ratio1_minus = gamma_ratio1;
    complex_type gamma_ratio2_minus = gamma_ratio2;

    sum = sum + my_pow(-1, 0) * ( gamma_ratio1/ gamma_ratio2 ) * rhsv_minimal_coefficient_n_positive[0];

    int consecutive_counts{0};
    
    while(consecutive_counts < constants::cons_counts)
    {    
        complex_type addend_1 {0};
        complex_type addend_2 {0};

        gamma_ratio1_plus = gamma_ratio1_plus*(arg_gamma_ratio_1 + i - 1);
        gamma_ratio2_plus = gamma_ratio2_plus*(arg_gamma_ratio_2 + i - 1);
        
        gamma_ratio1_minus = gamma_ratio1_minus/(arg_gamma_ratio_1 - i );
        gamma_ratio2_minus = gamma_ratio2_minus/(arg_gamma_ratio_2 - i );
        
        addend_1 = my_pow(-1, i) * ( gamma_ratio1_plus / gamma_ratio2_plus ) * rhsv_minimal_coefficient_n_positive[i];
        addend_2 = my_pow(-1, i) * ( gamma_ratio1_minus / gamma_ratio2_minus ) * rhsv_minimal_coefficient_n_negative[i];
        sum = sum + addend_1 + addend_2;
       
        if(mp::abs(addend_1/sum) > mp::abs(addend_2/sum))
        {
            accuracy = mp::abs(addend_1/sum);
        }
        else
        {
            accuracy = mp::abs(addend_2/sum);
        }
        accuracy = mp::abs(addend_1 + addend_2);

        if(accuracy < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }
        i++;

    }
    A_minus = sum * mp::pow(2, -1. - rhsv_spin_weight_s + constants::I*rhsv_epsilon) * mp::exp(- (constants::pi * rhsv_epsilon * 0.5) ) * mp::exp(- constants::pi*constants::I*(nu + 1. + rhsv_spin_weight_s)*0.5);   
}

void RadialHomogeneousSolution::rhsm_compute_K(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, complex_type& K)
{   
    complex_type z1 = 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon - constants::I*rhsv_tau;
    complex_type gamma1;
    special_functions::sf_compute_gamma_function_z(z1, gamma1);
    complex_type z2 = 2*nu + 2.;
    complex_type gamma2;
    special_functions::sf_compute_gamma_function_z(z2, gamma2);
    complex_type z3 = nu + 1. - rhsv_spin_weight_s + constants::I*rhsv_epsilon;
    complex_type gamma3;
    special_functions::sf_compute_gamma_function_z(z3, gamma3);
    complex_type z4 = nu + 1. +  constants::I*rhsv_tau;
    complex_type gamma4;
    special_functions::sf_compute_gamma_function_z(z4, gamma4);
    complex_type z5 = nu + 1. + rhsv_spin_weight_s + constants::I*rhsv_epsilon;
    complex_type gamma5;
    special_functions::sf_compute_gamma_function_z(z5, gamma5);


    complex_type num = mp::exp(constants::I*rhsv_epsilon*rhsv_kappa)*mp::pow(2.*rhsv_epsilon*rhsv_kappa, rhsv_spin_weight_s - nu ) * my_pow(2. , - rhsv_spin_weight_s) * gamma1 * gamma2;
    complex_type den = gamma3 * gamma4 * gamma5;

    float_type accuracy1 {10};
    complex_type sum1 {0};
    
    int consecutive_counts{0};
    
    complex_type sum1z1 = 2*nu + 1. ;
    complex_type sum1gamma1;
    special_functions::sf_compute_gamma_function_z(sum1z1, sum1gamma1);
    complex_type sum1z2 = 1;
    complex_type sum1gamma2 = 1.;
    complex_type sum1z3 = nu + 1 + rhsv_spin_weight_s + constants::I*rhsv_epsilon;
    complex_type sum1gamma3;
    special_functions::sf_compute_gamma_function_z(sum1z3, sum1gamma3);
    complex_type sum1z4 = nu + 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon;
    complex_type sum1gamma4;
    special_functions::sf_compute_gamma_function_z(sum1z4, sum1gamma4);
    complex_type sum1z5 = nu + 1. + constants::I*rhsv_tau;
    complex_type sum1gamma5;
    special_functions::sf_compute_gamma_function_z(sum1z5, sum1gamma5);
    complex_type sum1z6 = nu + 1. - constants::I*rhsv_tau;
    complex_type sum1gamma6;
    special_functions::sf_compute_gamma_function_z(sum1z6, sum1gamma6);

    complex_type addend1 =  ((sum1gamma1)/(sum1gamma2)) * ((sum1gamma3)/(sum1gamma4)) * ((sum1gamma5)/(sum1gamma6)) * rhsv_minimal_coefficient_n_positive[0];
    sum1 = sum1 + addend1;

    int n = 1;
    while(consecutive_counts < constants::cons_counts)
    {
        sum1gamma1 *= sum1z1;
        sum1gamma2 *= sum1z2;
        sum1gamma3 *= sum1z3;
        sum1gamma4 *= sum1z4;
        sum1gamma5 *= sum1z5;
        sum1gamma6 *= sum1z6;

        sum1z1 +=1.;
        sum1z2 +=1.;
        sum1z3 +=1.;
        sum1z4 +=1.;
        sum1z5 +=1.;
        sum1z6 +=1.;

        addend1 = my_pow(-1., n) * ((sum1gamma1)/(sum1gamma2)) * ((sum1gamma3)/(sum1gamma4)) * ((sum1gamma5)/(sum1gamma6)) * rhsv_minimal_coefficient_n_positive[n];
        accuracy1 = mp::abs(addend1/sum1);
        sum1 = sum1 + addend1;
        
        if(accuracy1 < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }
        n++;

    }

    float_type accuracy2 {10};
    complex_type sum2 {0};
    
    complex_type sum2z = 1.;
    complex_type sum2gamma = 1.;

    complex_type sum2z1 = 2*nu + 2.;
    complex_type sum2ratio1 = 1.;

    complex_type sum2z2 = nu + 1. + rhsv_spin_weight_s - constants::I*rhsv_epsilon;
    complex_type sum2ratio2 = 1.;

    complex_type sum2z3 = nu + 1. - rhsv_spin_weight_s + constants::I*rhsv_epsilon;
    complex_type sum2ratio3 = 1.;

    complex_type addend2 = ((1.)/(sum2gamma)) * (1./sum2ratio1) * (sum2ratio2/sum2ratio3) * rhsv_minimal_coefficient_n_negative[0];
    sum2 = sum2 + addend2;


    n = -1;

    consecutive_counts = 0;

    while(consecutive_counts < constants::cons_counts)
    {    
        sum2gamma *= sum2z;
        sum2z += 1.;

        sum2ratio1 = sum2ratio1/(sum2z1 + n );
        sum2ratio2 = sum2ratio2/(sum2z2 + n );
        sum2ratio3 = sum2ratio3/(sum2z3 + n );
       

        complex_type addend2 = my_pow(-1, n) * ((1.)/(sum2gamma)) * (1./sum2ratio1) * (sum2ratio2/sum2ratio3) * rhsv_minimal_coefficient_n_negative[-n];
        accuracy2 = mp::abs(addend2/sum2);
        sum2 = sum2 + addend2;
        
        if(accuracy2 < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }
        n--;

    }
    K = (num/den) * sum1 * 1./sum2;
}

void RadialHomogeneousSolution::rhsm_compute_A_plus(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, const complex_type& sum_minimal_coefficients, complex_type& A_plus)
{   
    complex_type arg_num = nu + 1. - rhsv_spin_weight_s + constants::I*rhsv_epsilon;
    complex_type gamma_num;
    special_functions::sf_compute_gamma_function_z(arg_num, gamma_num);
    complex_type arg_den = nu + 1. + rhsv_spin_weight_s - constants::I*rhsv_epsilon;
    complex_type gamma_den;
    special_functions::sf_compute_gamma_function_z(arg_den, gamma_den);

    A_plus = sum_minimal_coefficients * mp::pow(2, -1. + rhsv_spin_weight_s - constants::I*rhsv_epsilon) * mp::exp(- (constants::pi * rhsv_epsilon * 0.5) ) * mp::exp( constants::pi*constants::I*(nu + 1. - rhsv_spin_weight_s)*0.5)*gamma_num/gamma_den;
}

void RadialHomogeneousSolution::rhsm_compute_constants_quantities() 
{
    complex_type sum_min_coeff_nu;
    complex_type A_minus_nu;
    rhsv_minus_nu_root_minus_1 = - rhsv_nu_root - 1.;
    complex_type sum_min_coeff_nu1;
    complex_type A_minus_nu1;

    rhsv_minimal_coefficient_n_positive_nu.push_back(1.);
    rhsv_minimal_coefficient_n_negative_nu.push_back(1.);
    
    rhsv_minimal_coefficient_n_positive_nu1.push_back(1.);
    rhsv_minimal_coefficient_n_negative_nu1.push_back(1.);

    for(int i = 0; i < constants::imax_coefficients; ++i)
    {
        rhsm_compute_minimal_coefficient_n_positive(i, rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu);
    }

    for(int i = 0; i < constants::imax_coefficients; ++i)
    {
        rhsm_compute_minimal_coefficient_n_negative(i, rhsv_nu_root, rhsv_minimal_coefficient_n_negative_nu);
    }

    for(int i = 1; i < constants::imax_coefficients; ++ i)
    {
        rhsv_minimal_coefficient_n_negative_nu1.push_back(rhsv_minimal_coefficient_n_positive_nu[i]);
        rhsv_minimal_coefficient_n_positive_nu1.push_back(rhsv_minimal_coefficient_n_negative_nu[i]);
    }
    
    rhsm_sum_minimal_coefficients(rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu, sum_min_coeff_nu);
    
    rhsm_compute_A_minus(rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu, A_minus_nu);

    sum_min_coeff_nu1 = sum_min_coeff_nu;

    rhsm_compute_A_minus(rhsv_minus_nu_root_minus_1, rhsv_minimal_coefficient_n_positive_nu1, rhsv_minimal_coefficient_n_negative_nu1, A_minus_nu1);

    complex_type K_nu;
    rhsm_compute_K(rhsv_nu_root ,rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu, K_nu);

    complex_type K_nu1;
    rhsm_compute_K(rhsv_minus_nu_root_minus_1 ,rhsv_minimal_coefficient_n_positive_nu1, rhsv_minimal_coefficient_n_negative_nu1, K_nu1);

    complex_type A_plus_nu;
    rhsm_compute_A_plus(rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu, sum_min_coeff_nu, A_plus_nu);

    rhsv_B_trans = mp::pow( (rhsv_epsilon * rhsv_kappa )/( rhsv_angular_frequency_omega ) , 2*rhsv_spin_weight_s ) * 
                  mp::exp( constants::I * rhsv_kappa * (rhsv_epsilon + rhsv_tau)* 0.5 * (1 + (2*mp::log(rhsv_kappa))/(1+rhsv_kappa)) ) * 
                  sum_min_coeff_nu;


    rhsv_C_trans = mp::pow(rhsv_angular_frequency_omega , -1. - 2*rhsv_spin_weight_s ) * 
                  mp::exp(constants::I * ( rhsv_epsilon * mp::log(rhsv_epsilon) - (1 - rhsv_kappa)*0.5*rhsv_epsilon) ) * 
                  A_minus_nu; 

    rhsv_B_inc = (1./rhsv_angular_frequency_omega)*(K_nu - constants::I * mp::exp(- constants::I * constants::pi * rhsv_nu_root)*K_nu1*mp::sin(constants::pi*(rhsv_nu_root - rhsv_spin_weight_s + constants::I * rhsv_epsilon))/mp::sin(constants::pi*(rhsv_nu_root + rhsv_spin_weight_s - constants::I * rhsv_epsilon)))*A_plus_nu*mp::exp(-constants::I*(rhsv_epsilon*mp::log(rhsv_epsilon) - (1-rhsv_kappa)*0.5*rhsv_epsilon));

    rhsv_B_ref = mp::pow(rhsv_angular_frequency_omega, -1 - 2 * rhsv_spin_weight_s)*(K_nu + constants::I*mp::exp(constants::I*constants::pi*rhsv_nu_root) * K_nu1 )*A_minus_nu * mp::exp(constants::I*(rhsv_epsilon*mp::log(rhsv_epsilon) - (1-rhsv_kappa)*0.5*rhsv_epsilon));

    rhsv_B_inc = rhsv_B_inc/rhsv_B_trans;

    rhsv_B_ref = rhsv_B_ref/rhsv_B_trans;

    //BTrans  = Bhole 
    //Bref = Bout 
    //Bin = Binc
    //Dout = Cout
    //Din = Cref
    //Binf = Ctrans
}

//R_in and R_up functions MST method
void RadialHomogeneousSolution::rhsm_compute_R0(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, complex_type& R0, const float_type& x_var )
{
    int n = 0;

    complex_type a = - n - nu - constants::I * rhsv_tau;
    complex_type b = - n - nu - rhsv_spin_weight_s - constants::I * rhsv_epsilon;
    complex_type c = - 2*n - 2*nu; 

    complex_type value2f1;
    float_type arg2f1 = 1/(1 - x_var);

    special_functions::sf_compute_2F1(a, b, c , arg2f1, value2f1);

    complex_type gamma1;
    special_functions::sf_compute_gamma_function_z(1 - rhsv_spin_weight_s - constants::I * rhsv_epsilon - constants::I * rhsv_tau, gamma1);
    complex_type gamma2;
    special_functions::sf_compute_gamma_function_z(2*n+ 2*nu+1., gamma2);
    complex_type gamma3;
    special_functions::sf_compute_gamma_function_z(n + nu + 1. - constants::I*rhsv_tau, gamma3);
    complex_type gamma4;
    special_functions::sf_compute_gamma_function_z(n + nu + 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon, gamma4);



    complex_type z2 = 2*nu + 1.;
    complex_type z3 =  nu + 1. - constants::I*rhsv_tau;
    complex_type z4 = nu + 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon;
    
    complex_type z2_plus = z2;
    complex_type z3_plus = z3;
    complex_type z4_plus = z4;

    complex_type z2_minus = z2;
    complex_type z3_minus = z3;
    complex_type z4_minus = z4;

    complex_type gamma2_plus, gamma3_plus, gamma4_plus;
    gamma2_plus = gamma2;
    gamma3_plus = gamma3;
    gamma4_plus = gamma4;

    complex_type gamma2_minus, gamma3_minus, gamma4_minus;
    gamma2_minus = gamma2;
    gamma3_minus = gamma3;
    gamma4_minus = gamma4;

    complex_type coeff_sum = gamma1*gamma2/gamma3/gamma4; 
    complex_type x_term_sum_R0 = mp::pow(1-x_var, n)* value2f1;
    complex_type sum_R0 = x_term_sum_R0 * coeff_sum;
    

    complex_type addend_1 = 0;
    complex_type addend_2 = 0;
    
    n = 1;

    float_type accuracy {10};

    int consecutive_counts{0};
    
    while(consecutive_counts < constants::cons_counts) 
    {  
        a = - n - nu - constants::I * rhsv_tau;
        b = - n - nu - rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        c = - 2*n - 2*nu;   
        
        gamma3_plus *= z3_plus;
        gamma4_plus *= z4_plus;
        gamma2_plus = (z2_plus+1)*z2_plus*gamma2_plus;

        coeff_sum = rhsv_minimal_coefficient_n_positive[n]*gamma1*gamma2_plus/gamma3_plus/gamma4_plus;
        
        z2_plus +=2.;
        z3_plus += 1.;
        z4_plus += 1.;
        special_functions::sf_compute_2F1(a, b, c , arg2f1, value2f1);
        
        x_term_sum_R0 = mp::pow(1 - x_var, n)* value2f1;
        addend_1 = x_term_sum_R0 * coeff_sum;

        n = -n;

        a = - n - nu - constants::I * rhsv_tau;
        b = - n - nu - rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        c = - 2*n - 2*nu; 
        
        z2_minus -=2.;
        z3_minus -= 1.;
        z4_minus -= 1.;
        gamma3_minus = gamma3_minus/z3_minus;
        gamma4_minus = gamma4_minus/z4_minus;
        gamma2_minus = gamma2_minus/z2_minus/(z2_minus+1);

        coeff_sum = rhsv_minimal_coefficient_n_negative[-n]*gamma1*gamma2_minus/gamma3_minus/gamma4_minus;
        
        special_functions::sf_compute_2F1(a, b, c , arg2f1, value2f1);

        x_term_sum_R0 = mp::pow(1-x_var, n)* value2f1;
        addend_2 = x_term_sum_R0 * coeff_sum;

        n = -n;

        sum_R0 = sum_R0 + addend_1 + addend_2;

        accuracy = mp::abs(addend_1/sum_R0);
        if(mp::abs(addend_2/sum_R0) > accuracy)
        {
            accuracy = mp::abs(addend_2/sum_R0);
        }
           
        if(accuracy < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }
        n++;
    }
    complex_type factor_R0 = mp::exp(constants::I * rhsv_epsilon * rhsv_kappa * x_var )* 
         mp::pow(-x_var, - rhsv_spin_weight_s - constants::I * 0.5 * (rhsv_epsilon + rhsv_tau) )* 
         mp::pow( 1- x_var, constants::I * 0.5 * (rhsv_epsilon + rhsv_tau ) + nu);

    R0 = factor_R0 * sum_R0;
}

void RadialHomogeneousSolution::rhsm_compute_R0_and_1der(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, complex_type& R0, complex_type& R0_der, const float_type& x_var )
{
    int n = 0;

    complex_type a = - n - nu - constants::I * rhsv_tau;
    complex_type b = - n - nu - rhsv_spin_weight_s - constants::I * rhsv_epsilon;
    complex_type c = - 2*n - 2*nu; 

    complex_type value2f1;
    float_type arg2f1 = 1/(1 - x_var);
    complex_type value2f1_der;

    special_functions::sf_compute_2F1(a, b, c , arg2f1, value2f1);
    special_functions::sf_compute_2F1(a + 1,  b, c , arg2f1, value2f1_der);

    complex_type gamma1;
    special_functions::sf_compute_gamma_function_z(1 - rhsv_spin_weight_s - constants::I * rhsv_epsilon - constants::I * rhsv_tau, gamma1);
    complex_type gamma2;
    special_functions::sf_compute_gamma_function_z(2*n+ 2*nu+1., gamma2);
    complex_type gamma3;
    special_functions::sf_compute_gamma_function_z(n + nu + 1. - constants::I*rhsv_tau, gamma3);
    complex_type gamma4;
    special_functions::sf_compute_gamma_function_z(n + nu + 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon, gamma4);

    complex_type coeff_sum = gamma1*gamma2/gamma3/gamma4; 
    
    complex_type x_term_sum_R0 = mp::pow(1-x_var, n)* value2f1;
    complex_type sum_R0 = x_term_sum_R0 * coeff_sum;

    complex_type x_term_sum_R0_der =  mp::pow(1-x_var, -1+n)*((a + n)*value2f1 - a*value2f1_der);
    complex_type sum_R0_der = - coeff_sum * x_term_sum_R0_der; 

    complex_type addend_1 = 0;
    complex_type addend_2 = 0;
    complex_type addend_1_der = 0;
    complex_type addend_2_der = 0;

    n = 1;

    float_type accuracy {10};

    int consecutive_counts{0};
    
    while(consecutive_counts < constants::cons_counts) 
    {  
        a = - n - nu - constants::I * rhsv_tau;
        b = - n - nu - rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        c = - 2*n - 2*nu;
    
        special_functions::sf_compute_gamma_function_z(1 - rhsv_spin_weight_s - constants::I * rhsv_epsilon - constants::I * rhsv_tau, gamma1);
        special_functions::sf_compute_gamma_function_z(2*n+ 2*nu+1., gamma2);
        special_functions::sf_compute_gamma_function_z(n + nu + 1. - constants::I*rhsv_tau, gamma3);
        special_functions::sf_compute_gamma_function_z(n + nu + 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon, gamma4);
    
        coeff_sum = rhsv_minimal_coefficient_n_positive[n]*gamma1*gamma2/gamma3/gamma4;

        special_functions::sf_compute_2F1(a, b, c , arg2f1, value2f1);
        special_functions::sf_compute_2F1(a+1, b, c , arg2f1, value2f1_der);
        
        x_term_sum_R0 = mp::pow(1 - x_var, n)* value2f1;
        addend_1 = x_term_sum_R0 * coeff_sum;

        x_term_sum_R0_der =  mp::pow(1-x_var, -1+n)*((a + n)*value2f1 - a*value2f1_der);
        addend_1_der = - coeff_sum * x_term_sum_R0_der;

        n = -n;

        a = - n - nu - constants::I * rhsv_tau;
        b = - n - nu - rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        c = - 2*n - 2*nu; 
    
        special_functions::sf_compute_gamma_function_z(1 - rhsv_spin_weight_s - constants::I * rhsv_epsilon - constants::I * rhsv_tau, gamma1);
        special_functions::sf_compute_gamma_function_z(2*n+ 2*nu+1., gamma2);
        special_functions::sf_compute_gamma_function_z(n + nu + 1. - constants::I*rhsv_tau, gamma3);
        special_functions::sf_compute_gamma_function_z(n + nu + 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon, gamma4);
        
        coeff_sum = rhsv_minimal_coefficient_n_negative[-n]*gamma1*gamma2/gamma3/gamma4;

        special_functions::sf_compute_2F1(a, b, c , arg2f1, value2f1);
        special_functions::sf_compute_2F1(a+1, b, c , arg2f1,value2f1_der);

        x_term_sum_R0 = mp::pow(1-x_var, n)* value2f1;
        addend_2 = x_term_sum_R0 * coeff_sum;

        x_term_sum_R0_der =  mp::pow(1-x_var, -1+n)*((a + n)*value2f1 - a*value2f1_der);
        addend_2_der = - coeff_sum * x_term_sum_R0_der;

        n = -n;

        sum_R0 = sum_R0 + addend_1 + addend_2;
        sum_R0_der = sum_R0_der + addend_1_der + addend_2_der;

        accuracy = mp::abs(addend_1);
        if(mp::abs(addend_2) > accuracy)
        {
            accuracy = mp::abs(addend_2);
        }
        if(mp::abs(addend_1_der) > accuracy)
        {
            accuracy = mp::abs(addend_1_der);
        }
        if(mp::abs(addend_2_der) > accuracy)
        {
            accuracy = mp::abs(addend_2_der);
        }
           
        if(accuracy < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }
        n++;
    }

    complex_type factor_R0 = mp::exp(constants::I * rhsv_epsilon * rhsv_kappa * x_var )* 
         mp::pow(-x_var, - rhsv_spin_weight_s - constants::I * 0.5 * (rhsv_epsilon + rhsv_tau) )* 
         mp::pow( 1- x_var, constants::I * 0.5 * (rhsv_epsilon + rhsv_tau ) + nu);

    complex_type factor_R0_der =    0.5*mp::exp(constants::I*x_var*rhsv_epsilon*rhsv_kappa)*
                                    mp::pow(1-x_var, -1. + 0.5*rhsv_epsilon*constants::I + nu + 0.5*rhsv_tau*constants::I)*
                                    mp::pow(-x_var, -1. - rhsv_spin_weight_s - 0.5*constants::I*( rhsv_epsilon + rhsv_tau ))*
                                    (-2.*(-1. + x_var)*(rhsv_spin_weight_s - constants::I*x_var*rhsv_epsilon*rhsv_kappa) + 2*x_var*nu + constants::I*(rhsv_epsilon+rhsv_tau));

    R0 = factor_R0 * sum_R0;
    R0_der = (factor_R0_der * sum_R0 + factor_R0 * sum_R0_der)*(- rhsv_angular_frequency_omega/rhsv_epsilon/rhsv_kappa);   
}

void RadialHomogeneousSolution::rhsm_compute_R0_and_1der_and_2der(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, complex_type& R0, complex_type& R0_der, complex_type& R0_2der, const float_type& x_var )
{
    int n = 0;

    complex_type a = - n - nu - constants::I * rhsv_tau;
    complex_type b = - n - nu - rhsv_spin_weight_s - constants::I * rhsv_epsilon;
    complex_type c = - 2*n - 2*nu; 

    complex_type value2f1;
    float_type arg2f1 = 1/(1 - x_var);
    complex_type value2f1_der;

    special_functions::sf_compute_2F1(a, b, c , arg2f1, value2f1);
    special_functions::sf_compute_2F1(a + 1,  b, c , arg2f1, value2f1_der);

    complex_type gamma1;
    special_functions::sf_compute_gamma_function_z(1 - rhsv_spin_weight_s - constants::I * rhsv_epsilon - constants::I * rhsv_tau, gamma1);
    complex_type gamma2;
    special_functions::sf_compute_gamma_function_z(2*n+ 2*nu+1., gamma2);
    complex_type gamma3;
    special_functions::sf_compute_gamma_function_z(n + nu + 1. - constants::I*rhsv_tau, gamma3);
    complex_type gamma4;
    special_functions::sf_compute_gamma_function_z(n + nu + 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon, gamma4);

    complex_type coeff_sum = gamma1*gamma2/gamma3/gamma4; 
    
    complex_type x_term_sum_R0 = mp::pow(1-x_var, n)* value2f1;
    complex_type sum_R0 = x_term_sum_R0 * coeff_sum;

    complex_type x_term_sum_R0_der =  mp::pow(1-x_var, -1+n)*((a + n)*value2f1 - a*value2f1_der);
    complex_type sum_R0_der = - coeff_sum * x_term_sum_R0_der; 

    complex_type x_term_sum_R0_2der = (mp::pow(1 - x_var,-2 + n)*((a*(1 + a - c) + (-1 + n)*n*x_var + a*(-2 + c + 2*n)*x_var)*value2f1 - a*(1 + a + b - c + (-2 + c + 2*n)*x_var)*value2f1_der))/x_var;
    complex_type sum_R0_2der =  coeff_sum * x_term_sum_R0_2der;

    complex_type addend_1 = 0;
    complex_type addend_2 = 0;
    complex_type addend_1_der = 0;
    complex_type addend_2_der = 0;
    complex_type addend_1_2der = 0;
    complex_type addend_2_2der = 0;

    n = 1;

    float_type accuracy {10};

    int consecutive_counts{0};
    
    while(consecutive_counts < constants::cons_counts) 
    {  
        a = - n - nu - constants::I * rhsv_tau;
        b = - n - nu - rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        c = - 2*n - 2*nu;
    
        special_functions::sf_compute_gamma_function_z(1 - rhsv_spin_weight_s - constants::I * rhsv_epsilon - constants::I * rhsv_tau, gamma1);
        special_functions::sf_compute_gamma_function_z(2*n+ 2*nu+1., gamma2);
        special_functions::sf_compute_gamma_function_z(n + nu + 1. - constants::I*rhsv_tau, gamma3);
        special_functions::sf_compute_gamma_function_z(n + nu + 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon, gamma4);
    
        coeff_sum = rhsv_minimal_coefficient_n_positive[n]*gamma1*gamma2/gamma3/gamma4;

        special_functions::sf_compute_2F1(a, b, c , arg2f1, value2f1);
        special_functions::sf_compute_2F1(a+1, b, c , arg2f1, value2f1_der);
        
        x_term_sum_R0 = mp::pow(1 - x_var, n)* value2f1;
        addend_1 = x_term_sum_R0 * coeff_sum;

        x_term_sum_R0_der =  mp::pow(1-x_var, -1+n)*((a + n)*value2f1 - a*value2f1_der);
        addend_1_der = - coeff_sum * x_term_sum_R0_der;

        x_term_sum_R0_2der = (mp::pow(1 - x_var,-2 + n)*((a*(1 + a - c) + (-1 + n)*n*x_var + a*(-2 + c + 2*n)*x_var)*value2f1 - a*(1 + a + b - c + (-2 + c + 2*n)*x_var)*value2f1_der))/x_var;
        addend_1_2der =  coeff_sum * x_term_sum_R0_2der;
    
        n = -n;

        a = - n - nu - constants::I * rhsv_tau;
        b = - n - nu - rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        c = - 2*n - 2*nu; 
    
        special_functions::sf_compute_gamma_function_z(1 - rhsv_spin_weight_s - constants::I * rhsv_epsilon - constants::I * rhsv_tau, gamma1);
        special_functions::sf_compute_gamma_function_z(2*n+ 2*nu+1., gamma2);
        special_functions::sf_compute_gamma_function_z(n + nu + 1. - constants::I*rhsv_tau, gamma3);
        special_functions::sf_compute_gamma_function_z(n + nu + 1. - rhsv_spin_weight_s - constants::I*rhsv_epsilon, gamma4);
        
        coeff_sum = rhsv_minimal_coefficient_n_negative[-n]*gamma1*gamma2/gamma3/gamma4;

        special_functions::sf_compute_2F1(a, b, c , arg2f1, value2f1);
        special_functions::sf_compute_2F1(a+1, b, c , arg2f1,value2f1_der);

        x_term_sum_R0 = mp::pow(1-x_var, n)* value2f1;
        addend_2 = x_term_sum_R0 * coeff_sum;

        x_term_sum_R0_der =  mp::pow(1-x_var, -1+n)*((a + n)*value2f1 - a*value2f1_der);
        addend_2_der = - coeff_sum * x_term_sum_R0_der;

        x_term_sum_R0_2der = (mp::pow(1 - x_var,-2 + n)*((a*(1 + a - c) + (-1 + n)*n*x_var + a*(-2 + c + 2*n)*x_var)*value2f1 - a*(1 + a + b - c + (-2 + c + 2*n)*x_var)*value2f1_der))/x_var;
        addend_2_2der =  coeff_sum * x_term_sum_R0_2der;

        n = -n;

        sum_R0 = sum_R0 + addend_1 + addend_2;
        sum_R0_der = sum_R0_der + addend_1_der + addend_2_der;
        sum_R0_2der = sum_R0_2der + addend_1_2der + addend_2_2der;

        accuracy = mp::abs(addend_1);
        if(mp::abs(addend_2) > accuracy)
        {
            accuracy = mp::abs(addend_2);
        }
        if(mp::abs(addend_1_der) > accuracy)
        {
            accuracy = mp::abs(addend_1_der);
        }
        if(mp::abs(addend_2_der) > accuracy)
        {
            accuracy = mp::abs(addend_2_der);
        }
        if(mp::abs(addend_1_2der) > accuracy)
        {
            accuracy = mp::abs(addend_1_2der);
        }
        if(mp::abs(addend_2_2der) > accuracy)
        {
            accuracy = mp::abs(addend_2_2der);
        }
           
        if(accuracy < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }
        n++;
    }

    complex_type factor_R0 = mp::exp(constants::I * rhsv_epsilon * rhsv_kappa * x_var )* 
         mp::pow(-x_var, - rhsv_spin_weight_s - constants::I * 0.5 * (rhsv_epsilon + rhsv_tau) )* 
         mp::pow( 1- x_var, constants::I * 0.5 * (rhsv_epsilon + rhsv_tau ) + nu);

    complex_type factor_R0_der =    0.5*mp::exp(constants::I*x_var*rhsv_epsilon*rhsv_kappa)*
                                    mp::pow(1-x_var, -1. + 0.5*rhsv_epsilon*constants::I + nu + 0.5*rhsv_tau*constants::I)*
                                    mp::pow(-x_var, -1. - rhsv_spin_weight_s - 0.5*constants::I*( rhsv_epsilon + rhsv_tau ))*
                                    (-2.*(-1. + x_var)*(rhsv_spin_weight_s - constants::I*x_var*rhsv_epsilon*rhsv_kappa) + 2*x_var*nu + constants::I*(rhsv_epsilon+rhsv_tau));

    complex_type factor_R0_2der =   -0.25*(mp::exp(constants::I*x_var*rhsv_epsilon*rhsv_kappa)*mp::pow(1 - x_var,(constants::I*(4*constants::I + rhsv_epsilon - 2*constants::I*nu + rhsv_tau))/2.)*
                                    mp::pow(-x_var,-2 - rhsv_spin_weight_s - 0.5*constants::I*(rhsv_epsilon + rhsv_tau))*(-4*my_pow(rhsv_spin_weight_s,2.)*mp::pow(-1 + x_var,2.) + mp::pow(rhsv_epsilon + 
                                    2*(-1 + x_var)*x_var*rhsv_epsilon*rhsv_kappa,2.) - 4*mp::pow(x_var,2.)*(-1 + nu)*nu - 2*constants::I*rhsv_epsilon*(1 + 2*x_var*(-1 + nu + 2*(-1 + x_var)*x_var*rhsv_kappa*nu)) + 
                                    2*rhsv_epsilon*(1 + 2*(-1 + x_var)*x_var*rhsv_kappa)*rhsv_tau - 2*constants::I*(1 + 2*x_var*(-1 + nu))*rhsv_tau + mp::pow(rhsv_tau,2.) + 4*constants::I*rhsv_spin_weight_s*(-1 + x_var)*
                                    (-1*constants::I + rhsv_epsilon + x_var*(constants::I + 2*(-1 + x_var)*rhsv_epsilon*rhsv_kappa - 2*constants::I*nu) + rhsv_tau)));

    R0 = factor_R0 * sum_R0;
    R0_der = (factor_R0_der * sum_R0 + factor_R0 * sum_R0_der)*(- rhsv_angular_frequency_omega/rhsv_epsilon/rhsv_kappa);
    R0_2der = (factor_R0_2der*sum_R0 + 2*factor_R0_der*sum_R0_der + factor_R0*sum_R0_2der)*mp::pow( - rhsv_angular_frequency_omega/rhsv_epsilon/rhsv_kappa, 2);
   
}

void RadialHomogeneousSolution::rhsm_compute_Rup(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, const float_type& z_var, complex_type& R_up)
{
    complex_type sum_Rup {0};
    
    int n = 0;
    complex_type a = n + nu + 1 + rhsv_spin_weight_s - constants::I * rhsv_epsilon;
    complex_type b = 2*n + 2*nu + 2;
    complex_type z = -2.*constants::I*z_var;
    
    complex_type U1;
    special_functions::sf_compute_confluent_hypergeometric_function(a, b, z, U1);

    complex_type arg_gamma_ratio_1 = nu + 1 + rhsv_spin_weight_s - constants::I * rhsv_epsilon;
    complex_type arg_gamma_ratio_2 = nu + 1 - rhsv_spin_weight_s + constants::I * rhsv_epsilon;

    complex_type gamma_ratio1 = 1.;
    complex_type gamma_ratio2 = 1.;

    complex_type gamma_ratio1_plus = gamma_ratio1;
    complex_type gamma_ratio2_plus = gamma_ratio2;

    complex_type gamma_ratio1_minus = gamma_ratio1;
    complex_type gamma_ratio2_minus = gamma_ratio2;

    complex_type coeff_sum = mp::pow(constants::I, n)*gamma_ratio1/gamma_ratio2;
    complex_type z_term_Rup = mp::pow(2*z_var, n)*U1;
    sum_Rup = coeff_sum*z_term_Rup;

    complex_type addend_1 = 0;
    complex_type addend_2 = 0;

    n = 1;

    float_type accuracy {10};

    int consecutive_counts{0};
    while(consecutive_counts < constants::cons_counts)
    {   
        a = n + nu + 1. + rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        b = 2*n + 2*nu + 2;
       
        special_functions::sf_compute_confluent_hypergeometric_function(a, b, z, U1);
       
        gamma_ratio1_plus = gamma_ratio1_plus*(arg_gamma_ratio_1 + n - 1);
        gamma_ratio2_plus = gamma_ratio2_plus*(arg_gamma_ratio_2 + n - 1);
        
        coeff_sum = rhsv_minimal_coefficient_n_positive[n]*mp::pow(constants::I, n)*gamma_ratio1_plus/gamma_ratio2_plus;
        z_term_Rup = mp::pow(2*z_var, n)*U1;
        addend_1 = coeff_sum*z_term_Rup;

        n = -n;

        a = n + nu + 1 +rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        b = 2*n + 2.*nu + 2;

        special_functions::sf_compute_confluent_hypergeometric_function(a, b, z, U1);
        
        gamma_ratio1_minus = gamma_ratio1_minus/(arg_gamma_ratio_1 + n );
        gamma_ratio2_minus = gamma_ratio2_minus/(arg_gamma_ratio_2 + n );
        
        coeff_sum = rhsv_minimal_coefficient_n_negative[-n]*mp::pow(constants::I, n)*gamma_ratio1_minus/gamma_ratio2_minus;
        z_term_Rup = mp::pow(2*z_var, n)*U1;
        addend_2 = coeff_sum*z_term_Rup;

        n = -n;

        sum_Rup = sum_Rup + addend_1 + addend_2;

        accuracy = mp::abs(addend_1/sum_Rup);
        if(mp::abs(addend_2/sum_Rup) > accuracy)
        {
            accuracy = mp::abs(addend_2/sum_Rup);
        }
           
        if(accuracy < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }

        n++;
    }


    complex_type factor_Rup =   mp::pow(2, nu)*mp::exp( - constants::pi*rhsv_epsilon - constants::I*constants::pi*(nu+1.+rhsv_spin_weight_s ) + constants::I*z_var)*
                                mp::pow(z_var, nu + constants::I*(rhsv_epsilon+rhsv_tau)*0.5)*
                                mp::pow(z_var-rhsv_epsilon*rhsv_kappa, -rhsv_spin_weight_s-constants::I*(rhsv_epsilon+rhsv_tau)*0.5);

    R_up = factor_Rup*sum_Rup;
}

void RadialHomogeneousSolution::rhsm_compute_Rup_and_1der(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, const float_type& z_var, complex_type& R_up, complex_type& R_up_der)
{
    complex_type sum_Rup {0};
    complex_type sum_Rup_der {0};
    int n = 0;
    complex_type a = n + nu + 1 + rhsv_spin_weight_s - constants::I * rhsv_epsilon;
    complex_type b = 2*n + 2*nu + 2;
    complex_type z = -2.*constants::I*z_var;
   
    complex_type U1;
    special_functions::sf_compute_confluent_hypergeometric_function(a, b, z, U1);

    complex_type gamma_ratio1;
    special_functions::sf_compute_gamma_ratio(nu + 1 + rhsv_spin_weight_s - constants::I * rhsv_epsilon, n, gamma_ratio1);
    complex_type gamma_ratio2;
    special_functions::sf_compute_gamma_ratio(nu + 1 - rhsv_spin_weight_s + constants::I * rhsv_epsilon,n, gamma_ratio2);
    
    complex_type coeff_sum = mp::pow(constants::I, n)*gamma_ratio1/gamma_ratio2;
    complex_type z_term_Rup = mp::pow(2*z_var, n)*U1;
    sum_Rup = coeff_sum*z_term_Rup;

    complex_type U2;
    special_functions::sf_compute_confluent_hypergeometric_function(a + 1, b + 1, z, U2);
    complex_type z_term_Rup_der =   my_pow(2, n)*mp::pow(z_var,-1+n)*
                                        (n*U1 + 2.*constants::I*a*z_var* U2);
    sum_Rup_der = coeff_sum*z_term_Rup_der;

    complex_type addend_1 = 0;
    complex_type addend_2 = 0;
    complex_type addend_1_der = 0;
    complex_type addend_2_der = 0;

    n = 1;

    float_type accuracy {10};

    int consecutive_counts{0};
    while(consecutive_counts < constants::cons_counts)
    {   
        a = n + nu + 1. + rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        b = 2*n + 2*nu + 2;
        special_functions::sf_compute_confluent_hypergeometric_function(a, b, z, U1);
        special_functions::sf_compute_gamma_ratio(nu + 1 + rhsv_spin_weight_s - constants::I * rhsv_epsilon, n, gamma_ratio1);
        special_functions::sf_compute_gamma_ratio(nu + 1 - rhsv_spin_weight_s + constants::I * rhsv_epsilon,n, gamma_ratio2);
        
        coeff_sum = rhsv_minimal_coefficient_n_positive[n]*mp::pow(constants::I, n)*gamma_ratio1/gamma_ratio2;
        z_term_Rup = mp::pow(2*z_var, n)*U1;
        addend_1 = coeff_sum*z_term_Rup;

        special_functions::sf_compute_confluent_hypergeometric_function(a + 1, b+1, z, U2);
        z_term_Rup_der =    my_pow(2, n)*mp::pow(z_var,-1+n)*
                            (n*U1 + 2.*constants::I*a*z_var* U2);
        addend_1_der =coeff_sum*z_term_Rup_der;


        n = -n;

        a = n + nu + 1 +rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        b = 2*n + 2.*nu + 2;

        special_functions::sf_compute_gamma_ratio(nu + 1. + rhsv_spin_weight_s - constants::I * rhsv_epsilon, n, gamma_ratio1);
        special_functions::sf_compute_gamma_ratio(nu + 1. - rhsv_spin_weight_s + constants::I * rhsv_epsilon,n, gamma_ratio2);
        special_functions::sf_compute_confluent_hypergeometric_function(a, b, z, U1);

        coeff_sum = rhsv_minimal_coefficient_n_negative[-n]*mp::pow(constants::I, n)*gamma_ratio1/gamma_ratio2;
        z_term_Rup = mp::pow(2*z_var, n)*U1;
        addend_2 = coeff_sum*z_term_Rup;

        special_functions::sf_compute_confluent_hypergeometric_function(a + 1., b + 1., z, U2);
        z_term_Rup_der =    my_pow(2, n)*mp::pow(z_var,-1.+ n)*
                            (n*U1 + 2.*constants::I*a*z_var* U2);
        addend_2_der =coeff_sum*z_term_Rup_der;

        n = -n;

        sum_Rup = sum_Rup + addend_1 + addend_2;
        sum_Rup_der = sum_Rup_der + addend_1_der + addend_2_der;

        accuracy = mp::abs(addend_1);
        if(mp::abs(addend_2) > accuracy)
        {
            accuracy = mp::abs(addend_2);
        }
        if(mp::abs(addend_1_der) > accuracy)
        {
            accuracy = mp::abs(addend_1_der);
        }
        if(mp::abs(addend_2_der) > accuracy)
        {
            accuracy = mp::abs(addend_2_der);
        }
           
        if(accuracy < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }
        n++;
    }

    complex_type factor_Rup =   mp::pow(2, nu)*mp::exp( - constants::pi*rhsv_epsilon - constants::I*constants::pi*(nu+1.+rhsv_spin_weight_s ) + constants::I*z_var)*
                                mp::pow(z_var, nu + constants::I*(rhsv_epsilon+rhsv_tau)*0.5)*
                                mp::pow(z_var-rhsv_epsilon*rhsv_kappa, -rhsv_spin_weight_s-constants::I*(rhsv_epsilon+rhsv_tau)*0.5);

    complex_type factor_Rup_der =   mp::pow(2, -1. + nu )* mp::exp(constants::I*z_var - constants::pi*rhsv_epsilon - constants::I * constants::pi * (1 + rhsv_spin_weight_s + nu))*
                                    mp::pow(z_var, 0.5*constants::I*(2.*constants::I + rhsv_epsilon - 2.*constants::I*nu + rhsv_tau))*
                                    mp::pow(z_var - rhsv_epsilon*rhsv_kappa, -rhsv_spin_weight_s - 0.5*constants::I*(-2.*constants::I + rhsv_epsilon + rhsv_tau))*
                                    (-2*rhsv_spin_weight_s*z_var + 2.*constants::I*mp::pow(z_var,2.) + 2*z_var*(-constants::I*rhsv_epsilon*rhsv_kappa + nu) - constants::I*rhsv_epsilon*rhsv_kappa*(rhsv_epsilon-2.*constants::I*nu+rhsv_tau));

    R_up = factor_Rup*sum_Rup;
    R_up_der = (factor_Rup * sum_Rup_der + factor_Rup_der * sum_Rup)*(rhsv_angular_frequency_omega);
}

void RadialHomogeneousSolution::rhsm_compute_Rup_and_1der_and_2der(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, const float_type& z_var, complex_type& R_up, complex_type& R_up_der, complex_type& R_up_2der)
{
    complex_type sum_Rup {0};
    complex_type sum_Rup_der {0};
    complex_type sum_Rup_2der {0};
    int n = 0;
    complex_type a = n + nu + 1 + rhsv_spin_weight_s - constants::I * rhsv_epsilon;
    complex_type b = 2*n + 2*nu + 2;
    complex_type z = -2.*constants::I*z_var;
   
    complex_type U1;
    special_functions::sf_compute_confluent_hypergeometric_function(a, b, z, U1);

    complex_type gamma_ratio1;
    special_functions::sf_compute_gamma_ratio(nu + 1 + rhsv_spin_weight_s - constants::I * rhsv_epsilon, n, gamma_ratio1);
    complex_type gamma_ratio2;
    special_functions::sf_compute_gamma_ratio(nu + 1 - rhsv_spin_weight_s + constants::I * rhsv_epsilon,n, gamma_ratio2);
    
    complex_type coeff_sum = mp::pow(constants::I, n)*gamma_ratio1/gamma_ratio2;
    complex_type z_term_Rup = mp::pow(2*z_var, n)*U1;
    sum_Rup = coeff_sum*z_term_Rup;

    complex_type U2;
    special_functions::sf_compute_confluent_hypergeometric_function(a + 1, b + 1, z, U2);
    complex_type z_term_Rup_der =   my_pow(2, n)*mp::pow(z_var,-1+n)*
                                        (n*U1 + 2.*constants::I*a*z_var* U2);
    sum_Rup_der = coeff_sum*z_term_Rup_der;

    complex_type U3;
    special_functions::sf_compute_confluent_hypergeometric_function(a + 2, b + 2, z, U3);
    
    complex_type z_term_Rup_2der =  my_pow(2, n)*mp::pow(z_var, -2 + n)*((-1 + n)*n*U1 + 4*z_var*(1 + n + rhsv_spin_weight_s - constants::I * rhsv_epsilon + nu)*
                                    (constants::I * n * U2 - z_var * (2 + n + rhsv_spin_weight_s - constants::I * rhsv_epsilon + nu )* U3));

    
    sum_Rup_2der = coeff_sum*z_term_Rup_2der;

    complex_type addend_1 = 0;
    complex_type addend_2 = 0;
    complex_type addend_1_der = 0;
    complex_type addend_2_der = 0;
    complex_type addend_1_2der = 0;
    complex_type addend_2_2der = 0;

    n = 1;

    float_type accuracy {10};

    int consecutive_counts{0};
    while(consecutive_counts < constants::cons_counts)
    {   
        a = n + nu + 1. + rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        b = 2*n + 2*nu + 2;
        special_functions::sf_compute_confluent_hypergeometric_function(a, b, z, U1);
        special_functions::sf_compute_gamma_ratio(nu + 1 + rhsv_spin_weight_s - constants::I * rhsv_epsilon, n, gamma_ratio1);
        special_functions::sf_compute_gamma_ratio(nu + 1 - rhsv_spin_weight_s + constants::I * rhsv_epsilon,n, gamma_ratio2);
        
        coeff_sum = rhsv_minimal_coefficient_n_positive[n]*mp::pow(constants::I, n)*gamma_ratio1/gamma_ratio2;
        z_term_Rup = mp::pow(2*z_var, n)*U1;
        addend_1 = coeff_sum*z_term_Rup;

        special_functions::sf_compute_confluent_hypergeometric_function(a + 1, b+1, z, U2);
        z_term_Rup_der =    my_pow(2, n)*mp::pow(z_var,-1+n)*
                            (n*U1 + 2.*constants::I*a*z_var* U2);
        addend_1_der =coeff_sum*z_term_Rup_der;

        special_functions::sf_compute_confluent_hypergeometric_function(a + 2, b + 2, z, U3);
        z_term_Rup_2der =   my_pow(2, n)*mp::pow(z_var, -2 + n)*((-1 + n)*n*U1 + 4*z_var*(1 + n + rhsv_spin_weight_s - constants::I * rhsv_epsilon + nu)*
                            (constants::I * n * U2 - z_var * (2 + n + rhsv_spin_weight_s - constants::I * rhsv_epsilon + nu )* U3));
        addend_1_2der =coeff_sum*z_term_Rup_2der;

        n = -n;

        a = n + nu + 1 +rhsv_spin_weight_s - constants::I * rhsv_epsilon;
        b = 2*n + 2.*nu + 2;

        special_functions::sf_compute_gamma_ratio(nu + 1. + rhsv_spin_weight_s - constants::I * rhsv_epsilon, n, gamma_ratio1);
        special_functions::sf_compute_gamma_ratio(nu + 1. - rhsv_spin_weight_s + constants::I * rhsv_epsilon,n, gamma_ratio2);
        special_functions::sf_compute_confluent_hypergeometric_function(a, b, z, U1);

        coeff_sum = rhsv_minimal_coefficient_n_negative[-n]*mp::pow(constants::I, n)*gamma_ratio1/gamma_ratio2;
        z_term_Rup = mp::pow(2*z_var, n)*U1;
        addend_2 = coeff_sum*z_term_Rup;

        special_functions::sf_compute_confluent_hypergeometric_function(a + 1., b + 1., z, U2);
        z_term_Rup_der =    my_pow(2, n)*mp::pow(z_var,-1.+ n)*
                            (n*U1 + 2.*constants::I*a*z_var* U2);
        addend_2_der =coeff_sum*z_term_Rup_der;

        special_functions::sf_compute_confluent_hypergeometric_function(a + 2, b + 2, z, U3);
        z_term_Rup_2der =   my_pow(2, n)*mp::pow(z_var, -2 + n)*((-1 + n)*n*U1 + 4*z_var*(1 + n + rhsv_spin_weight_s - constants::I * rhsv_epsilon + nu)*
                            (constants::I * n * U2 - z_var * (2 + n + rhsv_spin_weight_s - constants::I * rhsv_epsilon + nu )* U3));
        addend_2_2der =coeff_sum*z_term_Rup_2der;

        n = -n;

        sum_Rup = sum_Rup + addend_1 + addend_2;
        sum_Rup_der = sum_Rup_der + addend_1_der + addend_2_der;
        sum_Rup_2der = sum_Rup_2der + addend_1_2der + addend_2_2der;

        accuracy = mp::abs(addend_1);
        if(mp::abs(addend_2) > accuracy)
        {
            accuracy = mp::abs(addend_2);
        }
        if(mp::abs(addend_1_der) > accuracy)
        {
            accuracy = mp::abs(addend_1_der);
        }
        if(mp::abs(addend_2_der) > accuracy)
        {
            accuracy = mp::abs(addend_2_der);
        }
        if(mp::abs(addend_1_2der) > accuracy)
        {
            accuracy = mp::abs(addend_1_2der);
        }
        if(mp::abs(addend_2_2der) > accuracy)
        {
            accuracy = mp::abs(addend_2_2der);
        }
           
        if(accuracy < constants::accuracy_sum_computation) 
        {
            consecutive_counts++;  
        } 
        else 
        {
            consecutive_counts = 0;  
        }
        n++;
    }


    complex_type factor_Rup =   mp::pow(2, nu)*mp::exp( - constants::pi*rhsv_epsilon - constants::I*constants::pi*(nu+1.+rhsv_spin_weight_s ) + constants::I*z_var)*
                                mp::pow(z_var, nu + constants::I*(rhsv_epsilon+rhsv_tau)*0.5)*
                                mp::pow(z_var-rhsv_epsilon*rhsv_kappa, -rhsv_spin_weight_s-constants::I*(rhsv_epsilon+rhsv_tau)*0.5);

    complex_type factor_Rup_der =   mp::pow(2, -1. + nu )* mp::exp(constants::I*z_var - constants::pi*rhsv_epsilon - constants::I * constants::pi * (1 + rhsv_spin_weight_s + nu))*
                                    mp::pow(z_var, 0.5*constants::I*(2.*constants::I + rhsv_epsilon - 2.*constants::I*nu + rhsv_tau))*
                                    mp::pow(z_var - rhsv_epsilon*rhsv_kappa, -rhsv_spin_weight_s - 0.5*constants::I*(-2.*constants::I + rhsv_epsilon + rhsv_tau))*
                                    (-2*rhsv_spin_weight_s*z_var + 2.*constants::I*mp::pow(z_var,2.) + 2*z_var*(-constants::I*rhsv_epsilon*rhsv_kappa + nu) - constants::I*rhsv_epsilon*rhsv_kappa*(rhsv_epsilon-2.*constants::I*nu+rhsv_tau));

    complex_type factor_Rup_2der =  mp::pow(2,-2 + nu)*mp::exp(constants::I*z_var - constants::pi*rhsv_epsilon - constants::I*constants::pi*(1 + rhsv_spin_weight_s + nu))*
                                    mp::pow(z_var,(constants::I*(4*constants::I + rhsv_epsilon - 2*constants::I*nu + rhsv_tau))/2.)*mp::pow(z_var - rhsv_epsilon*rhsv_kappa,- rhsv_spin_weight_s - 0.5*constants::I*(-4*constants::I + rhsv_epsilon + rhsv_tau))*
                                    (4*my_pow(rhsv_spin_weight_s,2.)*mp::pow(z_var,2.) - 4*mp::pow(z_var,4) + 8*mp::pow(z_var,3)*(rhsv_epsilon*rhsv_kappa + constants::I*nu) - 4*z_var*rhsv_epsilon*rhsv_kappa*(rhsv_epsilon*rhsv_kappa + constants::I*(-1 + nu))*
                                    (rhsv_epsilon - 2*constants::I*nu + rhsv_tau) - mp::pow(rhsv_epsilon,2.)*mp::pow(rhsv_kappa,2.)*(rhsv_epsilon - 2*constants::I*nu + rhsv_tau)*(2*constants::I + rhsv_epsilon - 2*constants::I*nu + rhsv_tau) + 
                                    4*rhsv_spin_weight_s*z_var*(-2*constants::I*mp::pow(z_var,2.) + z_var*(1 + 2*constants::I*rhsv_epsilon*rhsv_kappa - 2*nu) + constants::I*rhsv_epsilon*rhsv_kappa*(rhsv_epsilon - 2*constants::I*nu + rhsv_tau)) + 
                                    4*mp::pow(z_var,2.)*(-(mp::pow(rhsv_epsilon,2.)*(-1 + rhsv_kappa)*rhsv_kappa) + (-1 + nu)*nu + rhsv_epsilon*rhsv_kappa*(-4*constants::I*nu + rhsv_tau)));

    R_up = factor_Rup*sum_Rup;
    R_up_der = (factor_Rup * sum_Rup_der + factor_Rup_der * sum_Rup)*(rhsv_angular_frequency_omega);
    R_up_2der = mp::pow(rhsv_angular_frequency_omega, 2)*(factor_Rup*sum_Rup_2der + 2*factor_Rup_der*sum_Rup_der + factor_Rup_2der*sum_Rup);
}
 
void RadialHomogeneousSolution::rhsm_compute_R_in_mst(const float_type& r, complex_type& R_in) 
{
    float_type x_var;
    float_type r_plus = rhsv_mass_black_hole_M + mp::sqrt(mp::pow(rhsv_mass_black_hole_M,2.) - mp::pow(rhsv_black_hole_spin_a, 2));
    
    if(rhsv_angular_frequency_omega == 0)
    {
        x_var = (r_plus - r)/ 2./rhsv_mass_black_hole_M/rhsv_kappa;
    }
    else
    {
        x_var = rhsv_angular_frequency_omega*(r_plus - r)/rhsv_epsilon/rhsv_kappa;
    }    

    complex_type R0_nu;
    complex_type R0_nu1;
    
    rhsm_compute_R0(rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu, R0_nu, x_var);

    rhsm_compute_R0(rhsv_minus_nu_root_minus_1, rhsv_minimal_coefficient_n_positive_nu1, rhsv_minimal_coefficient_n_negative_nu1, R0_nu1, x_var);
    
    R_in = R0_nu + R0_nu1; 
    
    R_in = R_in/rhsv_B_trans;
}

void RadialHomogeneousSolution::rhsm_compute_R_in_and_1der_mst(const float_type& r, complex_type& R_in,complex_type& R_in_der) 
{
    float_type x_var;
    float_type r_plus = rhsv_mass_black_hole_M + mp::sqrt(mp::pow(rhsv_mass_black_hole_M,2.) - mp::pow(rhsv_black_hole_spin_a, 2));
    
    if(rhsv_angular_frequency_omega == 0)
    {
        x_var = (r_plus - r)/ 2./rhsv_mass_black_hole_M/rhsv_kappa;
    }
    else
    {
        x_var = rhsv_angular_frequency_omega*(r_plus - r)/rhsv_epsilon/rhsv_kappa;
    }    

    complex_type R0_nu;
    complex_type R0_der_nu;

    complex_type R0_nu1;
    complex_type R0_der_nu1;
            
    rhsm_compute_R0_and_1der(rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu, R0_nu, R0_der_nu, x_var);

    rhsm_compute_R0_and_1der(rhsv_minus_nu_root_minus_1, rhsv_minimal_coefficient_n_positive_nu1, rhsv_minimal_coefficient_n_negative_nu1, R0_nu1, R0_der_nu1, x_var);
    
    R_in = R0_nu + R0_nu1; 
    R_in_der = R0_der_nu + R0_der_nu1; 

    R_in = R_in/rhsv_B_trans;
    R_in_der = R_in_der/rhsv_B_trans;
}

void RadialHomogeneousSolution::rhsm_compute_R_in_and_1der_and_2der_mst(const float_type& r, complex_type& R_in,complex_type& R_in_der, complex_type& R_in_2der) 
{
    float_type x_var;
    float_type r_plus = rhsv_mass_black_hole_M + mp::sqrt(mp::pow(rhsv_mass_black_hole_M,2.) - mp::pow(rhsv_black_hole_spin_a, 2));
    
    if(rhsv_angular_frequency_omega == 0)
    {
        x_var = (r_plus - r)/ 2./rhsv_mass_black_hole_M/rhsv_kappa;
    }
    else
    {
        x_var = rhsv_angular_frequency_omega*(r_plus - r)/rhsv_epsilon/rhsv_kappa;
    }    

    complex_type R0_nu;
    complex_type R0_der_nu;
    complex_type R0_2der_nu;

    complex_type R0_nu1;
    complex_type R0_der_nu1;
    complex_type R0_2der_nu1;
            
    rhsm_compute_R0_and_1der_and_2der(rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu, R0_nu, R0_der_nu,R0_2der_nu, x_var);

    rhsm_compute_R0_and_1der_and_2der(rhsv_minus_nu_root_minus_1, rhsv_minimal_coefficient_n_positive_nu1, rhsv_minimal_coefficient_n_negative_nu1, R0_nu1, R0_der_nu1, R0_2der_nu1, x_var);
    
    R_in = R0_nu + R0_nu1; 
    R_in_der = R0_der_nu + R0_der_nu1;
    R_in_2der = R0_2der_nu + R0_2der_nu1; 

    R_in = R_in/rhsv_B_trans;
    R_in_der = R_in_der/rhsv_B_trans;
    R_in_2der = R_in_2der/rhsv_B_trans;
}

void RadialHomogeneousSolution::rhsm_compute_R_up_mst(const float_type& r, complex_type& R_up) 
{
    float_type r_minus = rhsv_mass_black_hole_M - mp::sqrt(mp::pow(rhsv_mass_black_hole_M,2.) - mp::pow(rhsv_black_hole_spin_a, 2));
    float_type z_var = rhsv_angular_frequency_omega*(r - r_minus );
            
    rhsm_compute_Rup(rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu , z_var, R_up);
    
    R_up = R_up/rhsv_C_trans;
}

void RadialHomogeneousSolution::rhsm_compute_R_up_and_1der_mst(const float_type& r, complex_type& R_up, complex_type& R_up_der) 
{

    float_type r_minus = rhsv_mass_black_hole_M - mp::sqrt(mp::pow(rhsv_mass_black_hole_M,2.) - mp::pow(rhsv_black_hole_spin_a, 2));
    float_type z_var = rhsv_angular_frequency_omega*(r - r_minus );
            
    rhsm_compute_Rup_and_1der(rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu , z_var, R_up, R_up_der);
    
    R_up = R_up/rhsv_C_trans;
    R_up_der = R_up_der/rhsv_C_trans;
}

void RadialHomogeneousSolution::rhsm_compute_R_up_and_1der_and_2der_mst(const float_type& r, complex_type& R_up, complex_type& R_up_der, complex_type& R_up_2der) 
{

    float_type r_minus = rhsv_mass_black_hole_M - mp::sqrt(mp::pow(rhsv_mass_black_hole_M,2.) - mp::pow(rhsv_black_hole_spin_a, 2));
    float_type z_var = rhsv_angular_frequency_omega*(r - r_minus );
            
    rhsm_compute_Rup_and_1der_and_2der(rhsv_nu_root, rhsv_minimal_coefficient_n_positive_nu, rhsv_minimal_coefficient_n_negative_nu , z_var, R_up, R_up_der, R_up_2der);
    
    R_up = R_up/rhsv_C_trans;
    R_up_der = R_up_der/rhsv_C_trans;
    R_up_2der = R_up_2der/rhsv_C_trans;
}

//R_in and R_up functions numerical integration method
void RadialHomogeneousSolution::rhsm_compute_coefficients_of_eq_numerical_integration(const int H, const float_type& r, complex_type& y2der, complex_type& y2, complex_type& y1)
{

    y2der = mp::pow(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r)*r,2)/mp::pow(r,4);

    y2 = (2*(mp::pow(rhsv_black_hole_spin_a, 2) + (-2 + r)*r)*(-mp::pow(rhsv_black_hole_spin_a, 2) + r + constants::I*rhsv_black_hole_spin_a*rhsv_m*r + r*rhsv_spin_weight_s - mp::pow(r,2)*rhsv_spin_weight_s + constants::I*H*r*(mp::pow(rhsv_black_hole_spin_a,2) + mp::pow(r,2))*rhsv_angular_frequency_omega))/mp::pow(r,5);

    y1 = ((mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r)*r)*(2*mp::pow(rhsv_black_hole_spin_a,2) - r*(2 + 2*rhsv_spin_weight_s + r*rhsv_angular_eingenvalue_lambda)) - 2*rhsv_black_hole_spin_a*(1 + H)*rhsv_m*mp::pow(r,2)*(mp::pow(rhsv_black_hole_spin_a,2) + mp::pow(r,2))*rhsv_angular_frequency_omega + 
    2.*constants::I*mp::pow(r,2)*(mp::pow(r,2)*(-3 + H + r - H*r) + mp::pow(rhsv_black_hole_spin_a,2)*(1 + H + r - H*r))*rhsv_spin_weight_s*rhsv_angular_frequency_omega - 
    (-1 + my_pow(H,2))*mp::pow(r,2)*mp::pow(mp::pow(rhsv_black_hole_spin_a,2) + mp::pow(r,2),2)*mp::pow(rhsv_angular_frequency_omega,2) - 2.*constants::I*rhsv_black_hole_spin_a*r*(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r)*r)*(rhsv_m + rhsv_black_hole_spin_a*H*rhsv_angular_frequency_omega))/mp::pow(r,6);
}

void RadialHomogeneousSolution::rhsm_compute_R_in_numerical_integration(const float_type& r_start_integration, const float_type& r_end_integration)
{
    int H = -1;
    float_type step_size = constants::step_size_of_numerical_integration;

    complex_type R_in_mst;
    complex_type R_in_der_mst;
    complex_type R_in_2der_mst;

    rhsm_compute_R_in_and_1der_mst(r_start_integration, R_in_mst , R_in_der_mst);
        
    float_type r_bc = r_start_integration;
    
    complex_type y1 = (mp::exp(constants::I*rhsv_angular_frequency_omega*(r_bc + ((1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*
                        mp::log((-1 + r_bc - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))/2.) - 
                       (1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*
                        mp::log((-1 + r_bc + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))/2.))/
                     mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2))))*r_bc*R_in_mst*
                mp::pow(-2*r_bc + mp::pow(r_bc,2) + mp::pow(rhsv_black_hole_spin_a,2),rhsv_spin_weight_s))
               /mp::pow((-1 + r_bc - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))/
                (-1 + r_bc + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2))),
               (constants::I*rhsv_m*rhsv_black_hole_spin_a)/(2.*mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2))));

                        
    complex_type y2 =   (mp::exp(constants::I*rhsv_angular_frequency_omega*(r_bc + (1 + 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc)/2.) + (1 - 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc)/2.)))*
                        mp::pow(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r_bc)*r_bc,-1 + rhsv_spin_weight_s)*(((-1.)*constants::I*rhsv_black_hole_spin_a*rhsv_m*r_bc + mp::pow(rhsv_black_hole_spin_a,2)*(1 + constants::I*r_bc*rhsv_angular_frequency_omega) + r_bc*(-2 + r_bc - 2*rhsv_spin_weight_s + 2*r_bc*rhsv_spin_weight_s + constants::I*mp::pow(r_bc,2)*rhsv_angular_frequency_omega))*
                        R_in_mst + r_bc*(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r_bc)*r_bc)*R_in_der_mst))/
                        mp::pow((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc)/(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc),(constants::I*rhsv_black_hole_spin_a*rhsv_m)/(2.*mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2))));
                       
                        std::array<complex_type, 2> x = {{y1, y2}};

               
    boost::numeric::odeint::bulirsch_stoer_dense_out<std::array<complex_type, 2>, float_type, std::array<complex_type, 2>, float_type> stepper(constants::error_tolerance_integration_of_numerical_integration, constants::error_tolerance_integration_of_numerical_integration, 1.0, 1.0);

    stepper.initialize(x, r_start_integration, step_size);
    
    while (stepper.current_time() < r_end_integration) {
        stepper.do_step(numerical_integration(*this, H));
        DenseStepData data;
        data.t_start = stepper.previous_time();
        data.dt = stepper.current_time() - stepper.previous_time();
        data.x = stepper.current_state(); 
        data.stepper_snapshot = stepper;  
        rhsv_dense_data_horizon.push_back(data);
    }

}

void RadialHomogeneousSolution::rhsm_compute_R_up_numerical_integration(const float_type& r_start_integration, const float_type& r_end_integration)
{
    int H = 1;
    float_type start = 0.;
    float_type end =  r_end_integration - r_start_integration;
    float_type step_size = constants::step_size_of_numerical_integration;

    complex_type R_up_mst;
    complex_type R_up_der_mst;
    complex_type R_up_2der_mst;
    rhsm_compute_R_up_and_1der_mst(r_end_integration, R_up_mst , R_up_der_mst);
    
    float_type r_bc = r_end_integration;
    
    complex_type y1 =   (r_bc*mp::pow(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r_bc)*r_bc,rhsv_spin_weight_s)*R_up_mst)/
    (mp::exp(constants::I*rhsv_angular_frequency_omega*(r_bc + (1 + 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc)/2.) + 
          (1 - 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc)/2.)))*
      mp::pow((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc)/(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc),
       (constants::I*rhsv_black_hole_spin_a*rhsv_m)/(2.*mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))));

    complex_type y2 =   (mp::pow(2,2*constants::I*rhsv_angular_frequency_omega)*mp::pow(-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc,-2 - constants::I*(1 + 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*rhsv_angular_frequency_omega)*
    mp::pow((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc)/(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc),
     1 - ((0.5)*constants::I*rhsv_black_hole_spin_a*rhsv_m)/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*
    mp::pow(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_bc,((-1.)*constants::I + constants::I/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*rhsv_angular_frequency_omega)*
    mp::pow(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r_bc)*r_bc,rhsv_spin_weight_s)*
    (((-1.)*constants::I*rhsv_black_hole_spin_a*rhsv_m*r_bc + mp::pow(rhsv_black_hole_spin_a,2)*(1 - constants::I*r_bc*rhsv_angular_frequency_omega) + 
         r_bc*(-2 + r_bc - 2*rhsv_spin_weight_s + 2*r_bc*rhsv_spin_weight_s - constants::I*mp::pow(r_bc,2)*rhsv_angular_frequency_omega))*R_up_mst + 
      r_bc*(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r_bc)*r_bc)*R_up_der_mst))/mp::exp(constants::I*r_bc*rhsv_angular_frequency_omega);

      y2 = -y2;
   
    std::array<complex_type, 2> x = {{y1, y2}};
   
    boost::numeric::odeint::bulirsch_stoer_dense_out<std::array<complex_type, 2>, float_type, std::array<complex_type, 2>, float_type> stepper(constants::error_tolerance_integration_of_numerical_integration, constants::error_tolerance_integration_of_numerical_integration, 1.0, 1.0);

    stepper.initialize(x, start, step_size);

    while (stepper.current_time() < end) {
        stepper.do_step(numerical_integration_infinity(*this, H, r_end_integration));
        DenseStepData data;
        data.t_start = stepper.previous_time();
        data.dt = stepper.current_time() - stepper.previous_time();
        data.x = stepper.current_state(); 
        data.stepper_snapshot = stepper;  
        rhsv_dense_data_infinity.push_back(data);
    }

}

//get functions
void RadialHomogeneousSolution::rhsm_get_nu(complex_type& res)
{
    res = rhsv_nu_root;
}

void RadialHomogeneousSolution::rhsm_get_B_trans(complex_type& res)
{
    res = rhsv_B_trans;
}

void RadialHomogeneousSolution::rhsm_get_C_trans(complex_type& res)
{
    res = rhsv_C_trans;
}

void RadialHomogeneousSolution::rhsm_get_B_inc(complex_type& res)
{
    res = rhsv_B_inc;
    if(rhsv_negative_freq == 1)
    {
        res = rhsv_B_inc.real() - rhsv_B_inc.imag()*constants::I;
    }
}

void RadialHomogeneousSolution::rhsm_get_B_ref(complex_type& res)
{
    res = rhsv_B_ref;
}

void RadialHomogeneousSolution::rhsm_get_R_in_mst(const float_type& r, complex_type& res)
{
    rhsm_compute_R_in_mst(r, res); 
}

void RadialHomogeneousSolution::rhsm_get_R_up_mst(const float_type& r, complex_type& res)
{
    rhsm_compute_R_up_mst(r, res); 
}

void RadialHomogeneousSolution::rhsm_get_R_in_and_1der_mst(const float_type&r, complex_type& R_in, complex_type& R_in_der)
{
    rhsm_compute_R_in_and_1der_mst(r, R_in, R_in_der); 
}

void RadialHomogeneousSolution::rhsm_get_R_up_and_1der_mst(const float_type&r, complex_type& R_up, complex_type& R_up_der)
{
    rhsm_compute_R_up_and_1der_mst(r, R_up, R_up_der); 
}

void RadialHomogeneousSolution::rhsm_get_R_in_and_1der_and_2der_mst(const float_type&r, complex_type& R_in, complex_type& R_in_der, complex_type& R_in_2der)
{
    rhsm_compute_R_in_and_1der_and_2der_mst(r, R_in, R_in_der, R_in_2der); 
}

void RadialHomogeneousSolution::rhsm_get_R_up_and_1der_and_2der_mst(const float_type&r, complex_type& R_up, complex_type& R_up_der, complex_type& R_up_2der)
{
     rhsm_compute_R_up_and_1der_and_2der_mst(r, R_up, R_up_der, R_up_2der);
}

void RadialHomogeneousSolution::rhsm_get_R_in_numerical_integration(const float_type& r_start_integration, const float_type& r_end_integration, const float_type&r, complex_type& R_in)
{
    if(rhsv_numerical_integration_in == 0)
    {
        rhsm_compute_R_in_numerical_integration(r_start_integration, r_end_integration);
        rhsv_numerical_integration_in = 1;
    }

    complex_type coeff = (mp::pow((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r),
    (constants::I*rhsv_black_hole_spin_a*rhsv_m)/(2.*mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2))))/
    (mp::exp(constants::I*rhsv_angular_frequency_omega*(r + (1 + 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/2.) + 
         (1 - 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/2.)))*r*
     mp::pow(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r)*r,rhsv_spin_weight_s)));

    std::array<complex_type, 2> x_interp;

    for (const auto& d : rhsv_dense_data_horizon) {
        if (r >= d.t_start && r <= d.t_start + d.dt) {
            d.stepper_snapshot.calc_state(r, x_interp);
            break;
        }
    }
    
    R_in = x_interp[0] * coeff; 

    if(rhsv_negative_freq == 1)
    {
        R_in = R_in.real() - R_in.imag()*constants::I;
    }
}

void RadialHomogeneousSolution::rhsm_get_R_up_numerical_integration(const float_type& r_start_integration, const float_type& r_end_integration, const float_type&r, complex_type& R_up)
{
    if(rhsv_numerical_integration_up == 0)
    {
        rhsm_compute_R_up_numerical_integration(r_start_integration, r_end_integration);
        rhsv_numerical_integration_up = 1;
    }

    float_type r_coeff = r;
    complex_type coeff = (mp::exp(constants::I*rhsv_angular_frequency_omega*(r_coeff + ((1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_coeff)/2.) - 
    (1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_coeff)/2.))/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2))))*
    mp::pow((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_coeff)/(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_coeff),
    (constants::I*rhsv_black_hole_spin_a*rhsv_m)/(2.*mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))))/(r_coeff*mp::pow(mp::pow(rhsv_black_hole_spin_a,2) - 2*r_coeff + mp::pow(r_coeff,2),rhsv_spin_weight_s));

    std::array<complex_type, 2> x_interp;

    float_type rho =  r_end_integration - r;
    for (const auto& d : rhsv_dense_data_infinity) {
        if (rho >= d.t_start && rho <= d.t_start + d.dt) {
            d.stepper_snapshot.calc_state(rho, x_interp);
            break;
        }
    }
    
    R_up = x_interp[0] * coeff; 
        
    if(rhsv_negative_freq == 1)
    {
        R_up = R_up.real() - R_up.imag()*constants::I;
    }
}

void RadialHomogeneousSolution::rhsm_get_R_in_numerical_integration_and_1der_and_2der(const float_type& r_start_integration, const float_type& r_end_integration, const float_type&r, complex_type& R_in, complex_type& R_in_der, complex_type& R_in_2der)
{
    if(rhsv_numerical_integration_in == 0)
    {
        rhsm_compute_R_in_numerical_integration(r_start_integration, r_end_integration);
        rhsv_numerical_integration_in = 1;
    }

    complex_type coeff = (mp::pow((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r),
    (constants::I*rhsv_black_hole_spin_a*rhsv_m)/(2.*mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2))))/
    (mp::exp(constants::I*rhsv_angular_frequency_omega*(r + (1 + 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/2.) + 
         (1 - 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/2.)))*r*
     mp::pow(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r)*r,rhsv_spin_weight_s)));

    std::array<complex_type, 2> x_interp;

    for (const auto& d : rhsv_dense_data_horizon) {
        if (r >= d.t_start && r <= d.t_start + d.dt) {
            d.stepper_snapshot.calc_state(r, x_interp);
            break;
        }
    }
    
    R_in = x_interp[0] * coeff; 

    complex_type coeff_der = (mp::pow(2,2*constants::I*rhsv_angular_frequency_omega)*mp::pow(-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r,-2 - constants::I*(1 + 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*rhsv_angular_frequency_omega)*
    mp::pow((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r),1 + ((0.5)*constants::I*rhsv_black_hole_spin_a*rhsv_m)/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*
    mp::pow(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r,((-1)*constants::I + constants::I/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*rhsv_angular_frequency_omega)*
    (constants::I*rhsv_black_hole_spin_a*rhsv_m*r + mp::pow(rhsv_black_hole_spin_a,2)*(-1 - constants::I*r*rhsv_angular_frequency_omega) + r*(2 - r + 2*rhsv_spin_weight_s - 2*r*rhsv_spin_weight_s - constants::I*mp::pow(r,2)*rhsv_angular_frequency_omega)))/
  (mp::exp(constants::I*r*rhsv_angular_frequency_omega)*mp::pow(r,2)*mp::pow(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r)*r,rhsv_spin_weight_s));
  
    R_in_der = x_interp[0]*coeff_der + x_interp[1]*coeff;

    R_in_2der = (-(R_in_der*(-2 + 2*r)*(1 + rhsv_spin_weight_s)) + R_in*((-4.)*constants::I*rhsv_angular_frequency_omega*r*rhsv_spin_weight_s - 
    (mp::pow(-(rhsv_black_hole_spin_a*rhsv_m) + rhsv_angular_frequency_omega*(mp::pow(rhsv_black_hole_spin_a,2) + mp::pow(r,2)),2) - 
       (2.)*constants::I*(-1 + r)*(-(rhsv_black_hole_spin_a*rhsv_m) + rhsv_angular_frequency_omega*(mp::pow(rhsv_black_hole_spin_a,2) + mp::pow(r,2)))*rhsv_spin_weight_s)/
     (mp::pow(rhsv_black_hole_spin_a,2) - 2*r + mp::pow(r,2)) + rhsv_angular_eingenvalue_lambda))/(mp::pow(rhsv_black_hole_spin_a,2) - 2*r + mp::pow(r,2));
    
    if(rhsv_negative_freq == 1)
    {
        R_in = R_in.real() - R_in.imag()*constants::I;
        R_in_der = R_in_der.real() - R_in_der.imag()*constants::I;
        R_in_2der = R_in_2der.real() - R_in_2der.imag()*constants::I;
    }
}

void RadialHomogeneousSolution::rhsm_get_R_up_numerical_integration_and_1der_and_2der(const float_type& r_start_integration, const float_type& r_end_integration, const float_type&r, complex_type& R_up, complex_type& R_up_der, complex_type& R_up_2der)
{
    if(rhsv_numerical_integration_up == 0)
    {
        rhsm_compute_R_up_numerical_integration(r_start_integration, r_end_integration);
        rhsv_numerical_integration_up = 1;
    }

    float_type r_coeff = r;
    complex_type coeff = (mp::exp(constants::I*rhsv_angular_frequency_omega*(r_coeff + ((1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_coeff)/2.) - 
    (1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_coeff)/2.))/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2))))*
    mp::pow((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_coeff)/(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r_coeff),
    (constants::I*rhsv_black_hole_spin_a*rhsv_m)/(2.*mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))))/(r_coeff*mp::pow(mp::pow(rhsv_black_hole_spin_a,2) - 2*r_coeff + mp::pow(r_coeff,2),rhsv_spin_weight_s));

    std::array<complex_type, 2> x_interp;

    float_type rho =  r_end_integration - r;
    for (const auto& d : rhsv_dense_data_infinity) {
        if (rho >= d.t_start && rho <= d.t_start + d.dt) {
            d.stepper_snapshot.calc_state(rho, x_interp);
            break;
        }
    }
    
    R_up = x_interp[0] * coeff; 
    
    complex_type coeff_der = (mp::exp(constants::I*rhsv_angular_frequency_omega*(r + (1 + 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/2.) + (1 - 1/mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)))*mp::log((-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/2.)))*
    mp::pow((-1 - mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r)/(-1 + mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2)) + r),(constants::I*rhsv_black_hole_spin_a*rhsv_m)/(2.*mp::sqrt(1 - mp::pow(rhsv_black_hole_spin_a,2))))*mp::pow(mp::pow(rhsv_black_hole_spin_a,2) + (-2 + r)*r,-1 - rhsv_spin_weight_s)*
    (constants::I*rhsv_black_hole_spin_a*rhsv_m*r + mp::pow(rhsv_black_hole_spin_a,2)*(-1 + constants::I*r*rhsv_angular_frequency_omega) + r*(2 - r + 2*rhsv_spin_weight_s - 2*r*rhsv_spin_weight_s + constants::I*mp::pow(r,2)*rhsv_angular_frequency_omega)))/mp::pow(r,2);

    R_up_der = x_interp[0]*coeff_der - x_interp[1]*coeff;

    R_up_2der = (-(R_up_der*(-2 + 2*r)*(1 + rhsv_spin_weight_s)) + R_up*((-4.)*constants::I*rhsv_angular_frequency_omega*r*rhsv_spin_weight_s - 
        (mp::pow(-(rhsv_black_hole_spin_a*rhsv_m) + rhsv_angular_frequency_omega*(mp::pow(rhsv_black_hole_spin_a,2) + mp::pow(r,2)),2) - 
           (2.)*constants::I*(-1 + r)*(-(rhsv_black_hole_spin_a*rhsv_m) + rhsv_angular_frequency_omega*(mp::pow(rhsv_black_hole_spin_a,2) + mp::pow(r,2)))*rhsv_spin_weight_s)/
         (mp::pow(rhsv_black_hole_spin_a,2) - 2*r + mp::pow(r,2)) + rhsv_angular_eingenvalue_lambda))/(mp::pow(rhsv_black_hole_spin_a,2) - 2*r + mp::pow(r,2));

    if(rhsv_negative_freq == 1)
    {
        R_up = R_up.real() - R_up.imag()*constants::I;
        R_up_der = R_up_der.real() - R_up_der.imag()*constants::I;
        R_up_2der = R_up_2der.real() - R_up_2der.imag()*constants::I;
    }
}

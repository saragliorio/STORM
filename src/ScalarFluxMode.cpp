#include "ScalarFluxMode.h"

ScalarFluxMode::ScalarFluxMode(int s, float_type M, float_type m, float_type a, float_type p, float_type e, float_type theta_inc, int ll, int mm, int kk, int nn) : 
        sfmv_spin_weight_s (s),
        sfmv_mass_black_hole_M (M), 
        sfmv_mass_test_particle_m (m), 
        sfmv_black_hole_spin_a (a),
        sfmv_semilatus_rectum_p (p),
        sfmv_eccentricity_e(e),
        sfmv_inclination_angle(theta_inc*constants::pi/180),
        sfmv_l(ll),
        sfmv_m(mm),
        sfmv_k(kk),
        sfmv_n(nn)
    {
        //IF s!=0 , a > 0.999 and a < -0.999, m > l, m, l, n, k not integers, inclination angle != , e >1 and e < 0, p > separatrix 
        
        if (mp::abs(sfmv_black_hole_spin_a) > 0.999) 
        {
            std::cerr << "Error in ScalarFluxMode class: black hole spin parameter a must be <= 0.999" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (mp::abs(sfmv_black_hole_spin_a) < -0.999) 
        {
            std::cerr << "Error in ScalarFluxMode class: black hole spin parameter a must be >= -0.999" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if (std::abs(sfmv_m) > sfmv_l) 
        {
            std::cerr << "Error in ScalarFluxMode class: mode index m must satisfy |m| <=l" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (sfmv_semilatus_rectum_p <= 0) 
        {
            std::cerr << "Error in ScalarFluxMode class: semilatus rectum parameter p must be > 0 (the method to compute the separatrix is missing, pay attention to use physical values of p (even when p > 0))" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (sfmv_eccentricity_e < 0) 
        {
            std::cerr << "Error in ScalarFluxMode class: eccentricity parameter e must be >= 0. " << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (sfmv_eccentricity_e > 1) 
        {
            std::cerr << "Error in ScalarFluxMode class: eccentricity parameter e must be < 1. " << std::endl;
            std::exit(EXIT_FAILURE);
        }


        sfmv_orbmot = std::make_unique<GeodesicOrbitalMotion>(M, m, a, p, e, theta_inc, ll, mm, kk, nn);
        
        sfmv_orbmot->gomm_get_constants_of_motion(sfmv_energy, sfmv_angular_momentum, sfmv_carter_constant);
        
        sfmv_orbmot->gomm_get_frequencies_mino_time(sfmv_Upsilon_r, sfmv_Upsilon_theta, sfmv_Upsilon_phi, sfmv_Gamma);
        
        sfmv_orbmot->gomm_get_frequencies_standard_time(sfmv_Omega_r, sfmv_Omega_theta, sfmv_Omega_phi);

        sfmv_orbmot->gomm_get_angular_frequency_omega(sfmv_angular_frequency_omega);


        sfmv_swh = std::make_unique<SpinWeightedSpheroidalHarmonics>(sfmv_spin_weight_s, sfmv_l, sfmv_m, sfmv_black_hole_spin_a*sfmv_angular_frequency_omega);
        sfmv_swh->swshm_get_lambda(sfmv_lambda);
              
        sfmv_r_min = sfmv_semilatus_rectum_p/(1 + sfmv_eccentricity_e);
        sfmv_r_max = sfmv_semilatus_rectum_p/(1 - sfmv_eccentricity_e);
    
        sfmv_homsol = std::make_unique<RadialHomogeneousSolution>(sfmv_spin_weight_s, sfmv_mass_black_hole_M, sfmv_black_hole_spin_a, sfmv_mass_test_particle_m, sfmv_lambda, sfmv_angular_frequency_omega, sfmv_l, sfmv_m);

        complex_type nu;    

        sfmv_homsol->rhsm_get_nu(nu);
        complex_type B_inc;
        sfmv_homsol->rhsm_get_B_inc(B_inc);

        if(mp::abs(sfmv_angular_frequency_omega)  <= constants::precision_threshold)
        {
            sfmv_W = - (2.*sfmv_l + 1); // to be verified
        }
        else
        {
            sfmv_W = 2.*constants::I*sfmv_angular_frequency_omega*B_inc;
        }
    }

void ScalarFluxMode::sfmm_compute_I_plus(const float_type& qr, complex_type& I_plus)
{  
    float_type tr;
    sfmv_orbmot->gomm_get_tr_function(qr, tr);
    float_type phi_r;
    sfmv_orbmot->gomm_get_phir_function(qr, phi_r);
    float_type r;
    sfmv_orbmot->gomm_get_r_function(qr, r);
    complex_type R_in;

    sfmv_homsol->rhsm_get_R_in_numerical_integration(sfmv_r_min, sfmv_r_max, r, R_in);
    
    I_plus = mp::cos(sfmv_n*qr + sfmv_angular_frequency_omega*tr-sfmv_m*phi_r)*R_in;
}

void ScalarFluxMode::sfmm_compute_I_minus(const float_type& qr, complex_type& I_minus)
{   
    float_type tr;
    sfmv_orbmot->gomm_get_tr_function(qr, tr);
    float_type phi_r;
    sfmv_orbmot->gomm_get_phir_function(qr, phi_r);
    float_type r;
    sfmv_orbmot->gomm_get_r_function(qr, r);

    complex_type R_up;

    sfmv_homsol->rhsm_get_R_up_numerical_integration(sfmv_r_min, sfmv_r_max, r, R_up);
    
    I_minus = mp::cos(sfmv_n*qr + sfmv_angular_frequency_omega*tr-sfmv_m*phi_r)*R_up;
}

void ScalarFluxMode::sfmm_compute_I_z(const float_type& qz, float_type& Iz)
{    
    float_type tz;
    sfmv_orbmot->gomm_get_tz_function(qz, tz);
    float_type phi_z;
    sfmv_orbmot->gomm_get_phiz_function(qz, phi_z);
    float_type z;
    sfmv_orbmot->gomm_get_z_function(qz, z);
   
    float_type phi = 0;
    
    complex_type S{0};
    sfmv_swh->swshm_get_S(z, phi, S);
    float_type S_real = S.real();
    
    Iz = mp::cos(sfmv_k*qz + sfmv_angular_frequency_omega*tz-sfmv_m*phi_z)*S_real;
}

void ScalarFluxMode::sfmm_compute_Iz1(float_type& Iz1)
{
    float_type qz_start = 0.;
    float_type qz_end = constants::pi;
    
    auto f_Iz1 = [this](float_type qz) 
    {
        float_type Iz;
        sfmm_compute_I_z(qz, Iz);
        return Iz; 
    };

    Iz1 = boost::math::quadrature::trapezoidal(f_Iz1, qz_start, qz_end, constants::precision_trapezoidal_quadrature);
}

void ScalarFluxMode::sfmm_compute_Iz2(float_type& Iz2)
{
    float_type qz_start = 0.;
    float_type qz_end = constants::pi;

    auto f_Iz2 = [this](float_type qz) 
    {
        float_type Iz;
        sfmm_compute_I_z(qz, Iz);

        float_type z;
        sfmv_orbmot->gomm_get_z_function(qz, z);

        return Iz*z*z; 
    };

    Iz2 = sfmv_black_hole_spin_a*sfmv_black_hole_spin_a*boost::math::quadrature::trapezoidal(f_Iz2, qz_start, qz_end, constants::precision_trapezoidal_quadrature);  
}

void ScalarFluxMode::sfmm_compute_Ir1_plus(complex_type& Ir1_plus)
{
    float_type qr_start = 0.;
    float_type qr_end = constants::pi;

    auto f_Ir1_plus = [this](float_type qr) 
    {
        complex_type I_plus;
        sfmm_compute_I_plus(qr, I_plus);
        return I_plus; 
    };

    Ir1_plus = boost::math::quadrature::trapezoidal(f_Ir1_plus, qr_start, qr_end, constants::precision_trapezoidal_quadrature);
}

void ScalarFluxMode::sfmm_compute_Ir2_plus(complex_type& Ir2_plus)
{
    float_type qr_start = 0.;
    float_type qr_end = constants::pi;

    auto f_Ir2_plus = [this](float_type qr) 
    {
        complex_type I_plus;
        sfmm_compute_I_plus(qr, I_plus);

        float_type r;
        sfmv_orbmot->gomm_get_r_function(qr, r);

        return I_plus*r*r; 
    };

    Ir2_plus = boost::math::quadrature::trapezoidal(f_Ir2_plus, qr_start, qr_end, constants::precision_trapezoidal_quadrature);
}

void ScalarFluxMode::sfmm_compute_Ir1_minus(complex_type& Ir1_minus)
{
    float_type qr_start = 0.;
    float_type qr_end = constants::pi;

    auto f_Ir1_minus = [this](float_type qr) 
    {
        complex_type I_minus;
        sfmm_compute_I_minus(qr, I_minus);
        return I_minus; 
    };

    Ir1_minus = boost::math::quadrature::trapezoidal(f_Ir1_minus, qr_start, qr_end, constants::precision_trapezoidal_quadrature);
    
}

void ScalarFluxMode::sfmm_compute_Ir2_minus(complex_type& Ir2_minus)
{
    float_type qr_start = 0.;
    float_type qr_end = constants::pi;

    auto f_Ir2_minus = [this](float_type qr) 
    {
        complex_type I_minus;
        sfmm_compute_I_minus(qr, I_minus);

        float_type r;
        sfmv_orbmot->gomm_get_r_function(qr, r);

        return I_minus*r*r; 
    };

    Ir2_minus = boost::math::quadrature::trapezoidal(f_Ir2_minus, qr_start, qr_end, constants::precision_trapezoidal_quadrature); 
}


//flux on generic orbit
void ScalarFluxMode::sfmm_compute_energy_flux_horizon_generic_orbit()
{
    float_type Iz1;
    sfmm_compute_Iz1(Iz1);

    float_type Iz2;
    sfmm_compute_Iz2(Iz2);

    complex_type Ir1_minus;
    sfmm_compute_Ir1_minus(Ir1_minus);

    complex_type Ir2_minus;
    sfmm_compute_Ir2_minus(Ir2_minus);

    complex_type coeff = (4./constants::pi/sfmv_Gamma/sfmv_W);
    complex_type delta_phi_minus = coeff * (Iz1*Ir2_minus + Iz2*Ir1_minus);
    
    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;
     
    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );
}

void ScalarFluxMode::sfmm_compute_energy_flux_infinity_generic_orbit()
{
    float_type Iz1;
    sfmm_compute_Iz1(Iz1);

    float_type Iz2;
    sfmm_compute_Iz2(Iz2);

    complex_type Ir1_plus;
    sfmm_compute_Ir1_plus(Ir1_plus);

    complex_type Ir2_plus;
    sfmm_compute_Ir2_plus(Ir2_plus);

    complex_type coeff = (4./constants::pi/sfmv_Gamma/sfmv_W);
    complex_type delta_phi_plus = coeff * (Iz1*Ir2_plus + Iz2*Ir1_plus);
         
    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
}

void ScalarFluxMode::sfmm_compute_energy_flux_horizon_infinity_generic_orbit()
{   
    float_type Iz1;
    sfmm_compute_Iz1(Iz1);

    float_type Iz2;
    sfmm_compute_Iz2(Iz2);

    complex_type Ir1_plus;
    sfmm_compute_Ir1_plus(Ir1_plus);

    complex_type Ir2_plus;
    sfmm_compute_Ir2_plus(Ir2_plus);

    complex_type Ir1_minus;
    sfmm_compute_Ir1_minus(Ir1_minus);

    complex_type Ir2_minus;
    sfmm_compute_Ir2_minus(Ir2_minus);

    complex_type coeff = (4./constants::pi/sfmv_Gamma/sfmv_W);
    complex_type delta_phi_plus = coeff * (Iz1*Ir2_plus + Iz2*Ir1_plus);
    complex_type delta_phi_minus = coeff * (Iz1*Ir2_minus + Iz2*Ir1_minus);
    
    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;
     
    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );

}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_generic_orbit(float_type& flux_horizon)
{
    if(mp::abs(sfmv_angular_frequency_omega)  <= constants::precision_threshold)
    {
        flux_horizon = 0;
    }
    else
    {
        if(sfmv_eccentricity_e <= constants::precision_threshold && sfmv_inclination_angle  <= constants::precision_threshold)
        {
            sfmm_compute_energy_flux_horizon_circular_orbit();
        }
        else if(sfmv_eccentricity_e > constants::precision_threshold && sfmv_inclination_angle  <= constants::precision_threshold)
        {
            sfmm_compute_energy_flux_horizon_eccentric_orbit();
        }
        else
        {
            sfmm_compute_energy_flux_horizon_generic_orbit();
        }
        flux_horizon = sfmv_scalar_energy_flux_horizon;
    }
}

void ScalarFluxMode::sfmm_get_energy_flux_infinity_generic_orbit(float_type& flux_infinity)
{
    if(mp::abs(sfmv_angular_frequency_omega)  <= constants::precision_threshold)
    {
        flux_infinity = 0;
    }
    else
    {
    if(sfmv_eccentricity_e <= constants::precision_threshold && sfmv_inclination_angle  <= constants::precision_threshold)
    {
        sfmm_compute_energy_flux_infinity_circular_orbit();
    }
    else if(sfmv_eccentricity_e > constants::precision_threshold && sfmv_inclination_angle  <= constants::precision_threshold)
    {
        sfmm_compute_energy_flux_infinity_eccentric_orbit();
    }
    else
    {
        sfmm_compute_energy_flux_infinity_generic_orbit();
    }
    flux_infinity = sfmv_scalar_energy_flux_infinity;
}
}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_infinity_generic_orbit(float_type& flux_horizon, float_type& flux_infinity)
{
    if(mp::abs(sfmv_angular_frequency_omega)  <= constants::precision_threshold)
    {
        flux_horizon = 0.;
        flux_infinity = 0.;
    }
    else
    {
        if(sfmv_eccentricity_e <= constants::precision_threshold && sfmv_inclination_angle  <= constants::precision_threshold)
        {
            sfmm_compute_energy_flux_horizon_infinity_circular_orbit();
        }
        else if(sfmv_eccentricity_e > constants::precision_threshold && sfmv_inclination_angle  <= constants::precision_threshold)
        {
            sfmm_compute_energy_flux_horizon_infinity_eccentric_orbit();
        }
        else
        {
            sfmm_compute_energy_flux_horizon_infinity_generic_orbit();
        }
        flux_horizon = sfmv_scalar_energy_flux_horizon;
        flux_infinity = sfmv_scalar_energy_flux_infinity;
    }
}


//flux on eccentric equatorial orbit
void ScalarFluxMode::sfmm_compute_energy_flux_horizon_eccentric_orbit()
{
    complex_type Ir2_minus;
    sfmm_compute_Ir2_minus(Ir2_minus);

    complex_type coeff = (4./constants::pi/sfmv_Gamma/sfmv_W);

    complex_type S{0};
    sfmv_swh->swshm_get_S(0., 0., S);
    float_type S_real = S.real();

    complex_type delta_phi_minus = coeff * (Ir2_minus*S_real*constants::pi);

    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;
     
    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );
}

void ScalarFluxMode::sfmm_compute_energy_flux_infinity_eccentric_orbit()
{
    complex_type Ir2_plus;
    sfmm_compute_Ir2_plus(Ir2_plus);

    complex_type coeff = (4./constants::pi/sfmv_Gamma/sfmv_W);

    complex_type S{0};
    sfmv_swh->swshm_get_S(0., 0., S);
    float_type S_real = S.real();

    complex_type delta_phi_plus = coeff * (Ir2_plus*S_real*constants::pi);
     
    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
}

void ScalarFluxMode::sfmm_compute_energy_flux_horizon_infinity_eccentric_orbit()
{
    complex_type Ir2_plus;
    sfmm_compute_Ir2_plus(Ir2_plus);

    complex_type Ir2_minus;
    sfmm_compute_Ir2_minus(Ir2_minus);

    complex_type coeff = (4./constants::pi/sfmv_Gamma/sfmv_W);

    complex_type S{0};
    sfmv_swh->swshm_get_S(0., 0., S);
    float_type S_real = S.real();

    complex_type delta_phi_plus = coeff * (Ir2_plus*S_real*constants::pi);
    complex_type delta_phi_minus = coeff * (Ir2_minus*S_real*constants::pi);

    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;
     
    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );
}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_eccentric_orbit(float_type& flux_horizon)
{
    sfmm_compute_energy_flux_horizon_eccentric_orbit();
    flux_horizon = sfmv_scalar_energy_flux_horizon;
}

void ScalarFluxMode::sfmm_get_energy_flux_infinity_eccentric_orbit(float_type& flux_infinity)
{
    sfmm_compute_energy_flux_infinity_eccentric_orbit();
    flux_infinity = sfmv_scalar_energy_flux_infinity;
}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_infinity_eccentric_orbit(float_type& flux_horizon, float_type& flux_infinity)
{
    sfmm_compute_energy_flux_horizon_infinity_eccentric_orbit();
    flux_horizon = sfmv_scalar_energy_flux_horizon;
    flux_infinity = sfmv_scalar_energy_flux_infinity;
}


//flux on circular equatorial orbit
void ScalarFluxMode::sfmm_compute_energy_flux_horizon_circular_orbit()
{
    float_type r_radius = sfmv_semilatus_rectum_p;
    
    complex_type R_up;
    sfmv_homsol->rhsm_get_R_up_mst(r_radius, R_up);
    
    float_type phi = 0;
    float_type x = 0.;  
    complex_type S{0};
    sfmv_swh->swshm_get_S(x, phi, S);

    complex_type alpha = - 4. * constants::pi * r_radius * r_radius * S/sfmv_Gamma/sfmv_W;
    

    complex_type delta_phi_minus = alpha*R_up;

    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;
     
    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );
}

void ScalarFluxMode::sfmm_compute_energy_flux_infinity_circular_orbit()
{
    float_type r_radius = sfmv_semilatus_rectum_p;
    complex_type R_in;
    sfmv_homsol->rhsm_get_R_in_mst(r_radius, R_in);
    
    float_type phi = 0;
    float_type x = 0.; 
    complex_type S{0};
    sfmv_swh->swshm_get_S(x, phi, S);

    complex_type alpha = - 4. * constants::pi * r_radius * r_radius * S/sfmv_Gamma/sfmv_W;
    

    complex_type delta_phi_plus = alpha*R_in;

     
    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
}

void ScalarFluxMode::sfmm_compute_energy_flux_horizon_infinity_circular_orbit()
{
    float_type r_radius = sfmv_semilatus_rectum_p;
    complex_type R_in;
    sfmv_homsol->rhsm_get_R_in_mst(r_radius, R_in);
    
    complex_type R_up;
    sfmv_homsol->rhsm_get_R_up_mst(r_radius, R_up);
    
    float_type phi = 0;
    float_type x = 0.; 
    complex_type S{0};
    sfmv_swh->swshm_get_S(x, phi, S);


    complex_type alpha = - 4. * constants::pi * r_radius * r_radius * S/sfmv_Gamma/sfmv_W;
    

    complex_type delta_phi_plus = alpha*R_in;
    complex_type delta_phi_minus = alpha*R_up;

    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;
     
    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );
}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_circular_orbit(float_type& flux_horizon)
{
    sfmm_compute_energy_flux_horizon_circular_orbit();
    flux_horizon = sfmv_scalar_energy_flux_horizon;
}

void ScalarFluxMode::sfmm_get_energy_flux_infinity_circular_orbit(float_type& flux_infinity)
{
    sfmm_compute_energy_flux_infinity_circular_orbit();
    flux_infinity = sfmv_scalar_energy_flux_infinity;
}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_infinity_circular_orbit(float_type& flux_horizon, float_type& flux_infinity)
{
    sfmm_compute_energy_flux_horizon_infinity_circular_orbit();
    flux_horizon = sfmv_scalar_energy_flux_horizon;
    flux_infinity = sfmv_scalar_energy_flux_infinity;
}
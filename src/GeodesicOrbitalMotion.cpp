#include "GeodesicOrbitalMotion.h"

GeodesicOrbitalMotion::GeodesicOrbitalMotion(float_type M, float_type m, float_type a, float_type p, float_type e, float_type theta_inc, int ll, int mm, int kk, int nn) : 
        gomv_mass_black_hole_M (M), 
        gomv_mass_test_particle_m (m), 
        gomv_black_hole_spin_a (a),
        gomv_semilatus_rectum_p (p),
        gomv_eccentricity_e(e),
        gomv_inclination_angle(theta_inc*constants::pi/180),
        gomv_l(ll),
        gomv_m(mm),
        gomv_k(kk),
        gomv_n(nn)
    {
        //IF a > 0.999 and a < -0.999, m > l, m, l, n, k not integers, inclination angle != , e >1 and e < 0, p > separatrix 
        if (mp::abs(gomv_black_hole_spin_a) > 0.999) 
        {
            std::cerr << "Error in GeodesicOrbitalMotion class: black hole spin parameter a must be <= 0.999" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (mp::abs(gomv_black_hole_spin_a) < -0.999) 
        {
            std::cerr << "Error in GeodesicOrbitalMotion class: black hole spin parameter a must be >= -0.999" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if (std::abs(gomv_m) > gomv_l) 
        {
            std::cerr << "Error in GeodesicOrbitalMotion class: mode index m must satisfy |m| <=l" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (gomv_semilatus_rectum_p <= 0) 
        {
            std::cerr << "Error in GeodesicOrbitalMotion class: semilatus rectum parameter p must be > 0 (the method to compute the separatrix is missing, pay attention to use physical values of p (even when p > 0))" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (gomv_eccentricity_e < 0) 
        {
            std::cerr << "Error in GeodesicOrbitalMotion class: eccentricity parameter e must be >= 0. " << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (gomv_eccentricity_e > 1) 
        {
            std::cerr << "Error in GeodesicOrbitalMotion class: eccentricity parameter e must be < 1. " << std::endl;
            std::exit(EXIT_FAILURE);
        }


        float_type signcos = mp::cos(gomv_inclination_angle)/mp::abs(mp::cos(gomv_inclination_angle));
        if(signcos>0)
        {
            gomv_prograde_retrograde_orbit = +1;
        }
        else
        {
            gomv_prograde_retrograde_orbit = -1;
        }

        gomv_theta_min = mp::abs(constants::pi*0.5 - gomv_inclination_angle);
        gomv_z1 = mp::cos(gomv_theta_min);
        gomv_r1 = gomv_semilatus_rectum_p/(1.- gomv_eccentricity_e);                                                                               
        gomv_r2 = gomv_semilatus_rectum_p/(1.+ gomv_eccentricity_e);                                                                               
    
        gomm_compute_constants_of_motion();
    
        float_type AplusB = 2./(1. - mp::pow(gomv_energy, 2) ) - gomv_r1 - gomv_r2;                                                                
        float_type AtimesB = mp::pow(gomv_black_hole_spin_a,2) * gomv_carter_constant / ( (1. - mp::pow(gomv_energy, 2) ) * gomv_r1 * gomv_r2);      
        gomv_r3 = 0.5*(AplusB + mp::pow ( mp::pow(AplusB, 2) - 4.*AtimesB, 0.5));                                                                
        gomv_r4 = AtimesB/gomv_r3;     
        
        gomv_kr =(gomv_r1-gomv_r2)/(gomv_r1-gomv_r3)*(gomv_r3-gomv_r4)/(gomv_r2-gomv_r4);                                 
    
        gomv_z2 = mp::sqrt(mp::pow(gomv_black_hole_spin_a,2)*(1.-gomv_energy*gomv_energy) + gomv_angular_momentum*gomv_angular_momentum/( 1. - gomv_z1*gomv_z1 ) );
        
        gomv_kz = gomv_black_hole_spin_a*gomv_black_hole_spin_a * (1.-gomv_energy*gomv_energy)*gomv_z1*gomv_z1/gomv_z2/gomv_z2;                                                          
        
        gomv_rp = gomv_mass_black_hole_M + mp::sqrt(gomv_mass_black_hole_M*gomv_mass_black_hole_M - gomv_black_hole_spin_a*gomv_black_hole_spin_a);
        gomv_rm = gomv_mass_black_hole_M - mp::sqrt(gomv_mass_black_hole_M*gomv_mass_black_hole_M - gomv_black_hole_spin_a*gomv_black_hole_spin_a);
        
        gomv_hr = (gomv_r1 - gomv_r2)/(gomv_r1 - gomv_r3);
        gomv_hp = gomv_hr*(gomv_r3 - gomv_rp)/(gomv_r2 - gomv_rp);
        gomv_hm = gomv_hr*(gomv_r3 - gomv_rm)/(gomv_r2 - gomv_rm);

        gomm_compute_frequencies_mino_time();
        gomm_compute_frequencies_standard_time();
        gomm_compute_angular_frequency_omega();
    }

void GeodesicOrbitalMotion::gomm_compute_delta_function(const float_type& r, float_type& delta)
{
    delta = mp::pow(r, 2)- 2.*r + mp::pow(gomv_black_hole_spin_a, 2);
}

void GeodesicOrbitalMotion::gomm_compute_constants_of_motion()
{   
    if(gomv_eccentricity_e == 0 && gomv_inclination_angle == 0)
    {
        gomv_energy = (gomv_black_hole_spin_a + (-2 + gomv_semilatus_rectum_p)*mp::sqrt(gomv_semilatus_rectum_p))/mp::sqrt(2*gomv_black_hole_spin_a*mp::pow(gomv_semilatus_rectum_p,1.5) + (-3 + gomv_semilatus_rectum_p)*mp::pow(gomv_semilatus_rectum_p,2));
        gomv_angular_momentum = (mp::pow(gomv_semilatus_rectum_p,2) - 2*mp::sqrt(gomv_semilatus_rectum_p)*gomv_black_hole_spin_a + mp::pow(gomv_black_hole_spin_a,2))/(mp::pow(gomv_semilatus_rectum_p,0.75)*mp::sqrt((-3 + gomv_semilatus_rectum_p)*mp::sqrt(gomv_semilatus_rectum_p) + 2*gomv_black_hole_spin_a));
        gomv_carter_constant = 0;
    }
    else if(gomv_eccentricity_e != 0 && gomv_inclination_angle == 0)
    {
        float_type F = (-4*mp::pow(gomv_black_hole_spin_a,2)*mp::pow(1 - mp::pow(gomv_eccentricity_e,2),2)*gomv_mass_black_hole_M + mp::pow(3 + mp::pow(gomv_eccentricity_e,2),2)*mp::pow(gomv_mass_black_hole_M,2)*gomv_semilatus_rectum_p - 
        2*(3 + mp::pow(gomv_eccentricity_e,2))*gomv_mass_black_hole_M*mp::pow(gomv_semilatus_rectum_p,2) + mp::pow(gomv_semilatus_rectum_p,3))/mp::pow(gomv_semilatus_rectum_p,3);

        float_type N = (2*(-(mp::pow(gomv_black_hole_spin_a,2)*(1 + 3*mp::pow(gomv_eccentricity_e,2))*gomv_mass_black_hole_M) + (-mp::pow(gomv_black_hole_spin_a,2) + (3 + mp::pow(gomv_eccentricity_e,2))*mp::pow(gomv_mass_black_hole_M,2))*gomv_semilatus_rectum_p - 
        gomv_mass_black_hole_M*mp::pow(gomv_semilatus_rectum_p,2)))/gomv_semilatus_rectum_p;

        float_type C = mp::pow(mp::pow(gomv_black_hole_spin_a,2) - gomv_mass_black_hole_M*gomv_semilatus_rectum_p,2);
        
        float_type x = mp::sqrt((- N - mp::sqrt(-4*C*F + N*N)*(gomv_black_hole_spin_a/mp::abs(gomv_black_hole_spin_a)))/2/F); 
        
        gomv_energy = mp::sqrt(1 - ((1 - mp::pow(gomv_eccentricity_e,2))*gomv_mass_black_hole_M*(1 - ((1 - mp::pow(gomv_eccentricity_e,2))*mp::pow(x,2))/mp::pow(gomv_semilatus_rectum_p,2)))/gomv_semilatus_rectum_p);

        gomv_angular_momentum = x + gomv_black_hole_spin_a*gomv_energy; 

        gomv_carter_constant = 0.;
    }
    else
    {
        float_type delta_max; 
        gomm_compute_delta_function(gomv_r1, delta_max);
        float_type d_max = delta_max*(mp::pow(gomv_r1,2)+mp::pow(gomv_black_hole_spin_a,2)*gomv_z1*gomv_z1); 
        float_type f_max = mp::pow(gomv_r1,4)+mp::pow(gomv_black_hole_spin_a, 2)*(gomv_r1*(gomv_r1+2) + gomv_z1*gomv_z1 * delta_max); 
        float_type g_max = 2*gomv_black_hole_spin_a*gomv_r1;  
        float_type h_max = gomv_r1*(gomv_r1-2) + (gomv_z1*gomv_z1* delta_max)/(1 - gomv_z1*gomv_z1); 

        float_type delta_min; 
        gomm_compute_delta_function(gomv_r2, delta_min);
        float_type d_min = delta_min*(mp::pow(gomv_r2,2)+mp::pow(gomv_black_hole_spin_a,2)*gomv_z1*gomv_z1); 
        float_type f_min = mp::pow(gomv_r2,4)+mp::pow(gomv_black_hole_spin_a, 2)*(gomv_r2*(gomv_r2+2) + gomv_z1*gomv_z1 * delta_min); 
        float_type g_min = 2*gomv_black_hole_spin_a*gomv_r2;   
        float_type h_min = gomv_r2*(gomv_r2 - 2) + (gomv_z1*gomv_z1* delta_min)/(1 - gomv_z1*gomv_z1); 
        
        float_type kappa = d_min * h_max - d_max * h_min;   
        float_type epsilon = d_min * g_max - d_max * g_min; 
        float_type rho = f_min * h_max - f_max * h_min;     
        float_type eta = f_min * g_max - f_max * g_min;     
        float_type sigma = g_min * h_max - g_max * h_min;   

        float_type sqrt_argument = sigma * ( sigma * mp::pow(epsilon, 2) + rho * epsilon * kappa - eta * mp::pow(kappa, 2));
        float_type num = kappa * rho + 2*epsilon*sigma - 2.*gomv_prograde_retrograde_orbit*mp::pow(sqrt_argument,0.5);
        float_type den = mp::pow(rho,2) + 4*eta * sigma;
        gomv_energy = mp::pow( num/den, 0.5);

        float_type addend_1 = mp::pow((g_min * gomv_energy)/( h_min),2);
        float_type addend_2 = (f_min * mp::pow( gomv_energy , 2) - d_min)/(h_min) ;
        float_type power_argument = addend_1 + addend_2;
        gomv_angular_momentum = - (g_min * gomv_energy)/(h_min) + gomv_prograde_retrograde_orbit * mp::pow(power_argument ,0.5);  

        gomv_beta = mp::pow(gomv_black_hole_spin_a, 2)*(1 - mp::pow(gomv_energy,2));
        gomv_carter_constant = gomv_z1*gomv_z1* (gomv_beta + (mp::pow(gomv_angular_momentum,2))/(1 - gomv_z1*gomv_z1));
    } 
}

void GeodesicOrbitalMotion::gomm_get_constants_of_motion(float_type& en, float_type& ang, float_type& crt)
{
    en = gomv_energy;
    ang = gomv_angular_momentum;
    crt = gomv_carter_constant;
}

void GeodesicOrbitalMotion::gomm_compute_frequencies_mino_time()
{
    gomv_Upsilon_r = (constants::pi * mp::pow( (1. - mp::pow(gomv_energy, 2))*(gomv_r1-gomv_r3)*(gomv_r2-gomv_r4) , 0.5))*0.5/boost::math::ellint_1(mp::sqrt(gomv_kr));     
    gomv_Upsilon_z = constants::pi*gomv_z2/2./boost::math::ellint_1(mp::sqrt(gomv_kz));
    
    float_type Upsilon_phi_plus = (2.*gomv_energy* gomv_rp - gomv_black_hole_spin_a * gomv_angular_momentum )/ (gomv_r3 - gomv_rp) * (1. - (gomv_r2 - gomv_r3)/(gomv_r2 - gomv_rp)* boost::math::ellint_3( mp::sqrt(gomv_kr), gomv_hp) / boost::math::ellint_1(mp::sqrt(gomv_kr)) ) ;
    float_type Upsilon_phi_minus = (2.*gomv_energy* gomv_rm - gomv_black_hole_spin_a * gomv_angular_momentum )/ (gomv_r3 - gomv_rm) * (1. - (gomv_r2 - gomv_r3)/(gomv_r2 - gomv_rm)* boost::math::ellint_3( mp::sqrt(gomv_kr), gomv_hm) / boost::math::ellint_1(mp::sqrt(gomv_kr)) ) ;
    float_type Upsilon_phi_r = gomv_black_hole_spin_a/(gomv_rp-gomv_rm)*(Upsilon_phi_plus - Upsilon_phi_minus);
    float_type Upsilon_phi_z = gomv_angular_momentum/boost::math::ellint_1(mp::sqrt(gomv_kz))*boost::math::ellint_3( mp::sqrt(gomv_kz), gomv_z1*gomv_z1);
    gomv_Upsilon_phi = Upsilon_phi_r + Upsilon_phi_z;
    
    float_type Gamma_r_plus = ((4. - gomv_black_hole_spin_a * gomv_angular_momentum/gomv_energy)*gomv_rp - 2.*gomv_black_hole_spin_a*gomv_black_hole_spin_a)/(gomv_r3 - gomv_rp)*(1.-(gomv_r2 - gomv_r3)/(gomv_r2 - gomv_rp)*boost::math::ellint_3( mp::sqrt(gomv_kr), gomv_hp)/boost::math::ellint_1(mp::sqrt(gomv_kr)));
    float_type Gamma_r_minus = ((4. - gomv_black_hole_spin_a * gomv_angular_momentum/gomv_energy)*gomv_rm - 2.*gomv_black_hole_spin_a*gomv_black_hole_spin_a)/(gomv_r3 - gomv_rm)*(1.-(gomv_r2 - gomv_r3)/(gomv_r2 - gomv_rm)*boost::math::ellint_3( mp::sqrt(gomv_kr), gomv_hm)/boost::math::ellint_1(mp::sqrt(gomv_kr)));
    float_type Gamma_r = (4. + gomv_black_hole_spin_a*gomv_black_hole_spin_a)*gomv_energy + gomv_energy*( 0.5 *   ((4. + gomv_r1 + gomv_r2 + gomv_r3)*gomv_r3 - gomv_r1*gomv_r2 + (gomv_r1 - gomv_r3)*(gomv_r2-gomv_r4)*boost::math::ellint_2( mp::sqrt(gomv_kr)) /boost::math::ellint_1(mp::sqrt(gomv_kr)) +
                        + (4. + gomv_r1 + gomv_r2 + gomv_r3 + gomv_r4)*(gomv_r2 - gomv_r3)*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hr)/boost::math::ellint_1(mp::sqrt(gomv_kr)) ) + 2./(gomv_rp-gomv_rm) *(Gamma_r_plus - Gamma_r_minus));
    
    float_type Gamma_z;

    if(gomv_carter_constant == 0)
    {
        Gamma_z = - gomv_black_hole_spin_a*gomv_black_hole_spin_a*gomv_energy;
    }
    else
    {
        Gamma_z = - gomv_black_hole_spin_a*gomv_black_hole_spin_a*gomv_energy + (gomv_energy*gomv_carter_constant)/(1. - gomv_energy*gomv_energy)/(gomv_z1*gomv_z1)* (1. - boost::math::ellint_2(mp::sqrt(gomv_kz))/boost::math::ellint_1(mp::sqrt(gomv_kz)));
    }
    gomv_Gamma = Gamma_r + Gamma_z;
}

void GeodesicOrbitalMotion::gomm_compute_frequencies_standard_time()
{
    gomv_Omega_r = gomv_Upsilon_r/gomv_Gamma;
    gomv_Omega_theta = gomv_Upsilon_z/gomv_Gamma;
    gomv_Omega_phi = gomv_Upsilon_phi/gomv_Gamma;
}

void GeodesicOrbitalMotion::gomm_get_frequencies_mino_time(float_type& Ur, float_type& Uz, float_type& Uphi, float_type& G)
{
    gomm_compute_frequencies_mino_time();
    Ur = gomv_Upsilon_r;
    Uz = gomv_Upsilon_z;
    Uphi = gomv_Upsilon_phi;
    G = gomv_Gamma;
}

void GeodesicOrbitalMotion::gomm_get_frequencies_standard_time(float_type& Or, float_type& Otheta, float_type& Ophi)
{
    gomm_compute_frequencies_standard_time();
    Or = gomv_Omega_r;
    Otheta = gomv_Omega_theta;
    Ophi = gomv_Omega_phi;
}

void GeodesicOrbitalMotion::gomm_compute_angular_frequency_omega()
{
    gomv_angular_frequency_omega = gomv_k * gomv_Omega_theta + gomv_n * gomv_Omega_r + gomv_m * gomv_Omega_phi;
}

void GeodesicOrbitalMotion::gomm_get_angular_frequency_omega(float_type& omega)
{
    omega = gomv_angular_frequency_omega;
}

void GeodesicOrbitalMotion::gomm_compute_r_function(const float_type& qr, float_type& r)
{
    float_type u = boost::math::ellint_1(mp::sqrt(gomv_kr))*qr/constants::pi;
    float_type sn = boost::math::jacobi_sn(mp::sqrt(gomv_kr), u);
    r = (gomv_r3*(gomv_r1 - gomv_r2)*mp::pow(sn,2) - gomv_r2*(gomv_r1-gomv_r3))/
        ((gomv_r1 - gomv_r2)*mp::pow(sn,2) - (gomv_r1-gomv_r3));
}

void GeodesicOrbitalMotion::gomm_compute_z_function(const float_type& qz, float_type& z)
{
    float_type u = 2*boost::math::ellint_1(mp::sqrt(gomv_kz))*(qz+constants::pi/2.)/constants::pi;
    z = gomv_z1*boost::math::jacobi_sn(mp::sqrt(gomv_kz), u);
}

void GeodesicOrbitalMotion::gomm_compute_tr_function(const float_type& qr, float_type& tr)
{
    if(mp::abs(qr - constants::pi) < constants::precision_q)
    {
        tr = 0;
    }
    else
    {
        float_type y = constants::pi;
        float_type tilde_tr_plus = (gomv_rp*(4. - gomv_black_hole_spin_a * gomv_angular_momentum / gomv_energy) - 2.*gomv_black_hole_spin_a*gomv_black_hole_spin_a)/(gomv_r2 - gomv_rp)/(gomv_r3-gomv_rp)*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hp, y);
        float_type tilde_tr_minus = (gomv_rm*(4. - gomv_black_hole_spin_a * gomv_angular_momentum / gomv_energy) - 2.*gomv_black_hole_spin_a*gomv_black_hole_spin_a)/(gomv_r2 - gomv_rm)/(gomv_r3-gomv_rm)*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hm, y);
        float_type tilde_tr_pi =  gomv_energy*(gomv_r2 - gomv_r3)/mp::sqrt((1. - gomv_energy*gomv_energy)*(gomv_r1 - gomv_r3)*(gomv_r2 - gomv_r4))*( (4. + gomv_r1 + gomv_r2 + gomv_r3 + gomv_r4)*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hr, y) - 4./(gomv_rp - gomv_rm)*(tilde_tr_plus -tilde_tr_minus) +
                    (gomv_r1 - gomv_r3)*(gomv_r2 - gomv_r4)/(gomv_r2 - gomv_r3)* (2*boost::math::ellint_2(mp::sqrt(gomv_kr)) - gomv_hr * mp::sin(y)*mp::cos(y) * mp::sqrt(1. - gomv_kr*mp::pow(mp::sin(y),2))/(1. - gomv_hr*mp::pow( mp::sin(y),2))  ) );
    
        float_type u = boost::math::ellint_1(mp::sqrt(gomv_kr))*qr/constants::pi;
        
        float_type sn = boost::math::jacobi_sn(mp::sqrt(gomv_kr), u);
        if(sn >= 0)
        {
            float_type cn = boost::math::jacobi_cn(mp::sqrt(gomv_kr), u);
            y = mp::acos(cn);
        }
        else
        {
            float_type cn = boost::math::jacobi_cn(mp::sqrt(gomv_kr), u);
            y = 2.*constants::pi - mp::acos(cn);   
        }

        tilde_tr_plus = (gomv_rp*(4. - gomv_black_hole_spin_a * gomv_angular_momentum / gomv_energy) - 2.*gomv_black_hole_spin_a*gomv_black_hole_spin_a)/(gomv_r2 - gomv_rp)/(gomv_r3-gomv_rp)*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hp, y);
        tilde_tr_minus = (gomv_rm*(4. - gomv_black_hole_spin_a * gomv_angular_momentum / gomv_energy) - 2.*gomv_black_hole_spin_a*gomv_black_hole_spin_a)/(gomv_r2 - gomv_rm)/(gomv_r3-gomv_rm)*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hm, y);
        float_type tilde_tr_y =  gomv_energy*(gomv_r2 - gomv_r3)/mp::sqrt((1. - gomv_energy*gomv_energy)*(gomv_r1 - gomv_r3)*(gomv_r2 - gomv_r4))*( (4. + gomv_r1 + gomv_r2 + gomv_r3 + gomv_r4)*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hr, y) - 4./(gomv_rp - gomv_rm)*(tilde_tr_plus -tilde_tr_minus) +
                    (gomv_r1 - gomv_r3)*(gomv_r2 - gomv_r4)/(gomv_r2 - gomv_r3)* (boost::math::ellint_2(mp::sqrt(gomv_kr),y) - gomv_hr * mp::sin(y)*mp::cos(y) * mp::sqrt(1. - gomv_kr*mp::pow(mp::sin(y),2))/(1. - gomv_hr*mp::pow( mp::sin(y),2))  ) );
    
        tr = tilde_tr_y - tilde_tr_pi/2./constants::pi*qr;
    }

}

void GeodesicOrbitalMotion::gomm_compute_tz_function(const float_type& qz, float_type& tz)
{
    if(mp::abs(qz - constants::pi) < constants::precision_q || mp::abs(qz - constants::pi*0.5) < constants::precision_q || mp::abs(qz - 0.) < constants::precision_q)
    {
        tz = 0.;
    }
    else
    {
        float_type y = constants::pi;
        float_type tilde_tz_pi = - gomv_energy/(1. - gomv_energy*gomv_energy)*gomv_z2*2*boost::math::ellint_2(sqrt(gomv_kz));

       
        float_type u = 2.*boost::math::ellint_1(mp::sqrt(gomv_kz))* (qz+constants::pi/2.)/constants::pi;

        float_type sn = boost::math::jacobi_sn(mp::sqrt(gomv_kz), u);

        if(sn >= 0)
        {
            float_type cn = boost::math::jacobi_cn(mp::sqrt(gomv_kz), u);
            y = mp::acos(cn);

        }
        else
        {
            float_type cn = boost::math::jacobi_cn(mp::sqrt(gomv_kz), u);
            y = 2.*constants::pi - mp::acos(cn);   

        }


        float_type tilde_tz_y = - gomv_energy/(1. - gomv_energy*gomv_energy)*gomv_z2*boost::math::ellint_2(sqrt(gomv_kz), y);

        tz = tilde_tz_y - tilde_tz_pi/constants::pi*(qz+constants::pi/2.);
    }
} 

void GeodesicOrbitalMotion::gomm_compute_t_function(const float_type& qr, const float_type& qz, const float_type& qt, float_type& t)
{
    float_type tr;
    gomm_compute_tr_function(qr, tr);

    float_type tz;
    gomm_compute_tz_function(qz, tz);

    t = qt + tr + tz;
}

void GeodesicOrbitalMotion::gomm_compute_phir_function(const float_type& qr, float_type& phi_r)
{
    if(mp::abs(qr - constants::pi) < constants::precision_q)
    {
        phi_r = 0;
    }
    else
    {
        float_type y = constants::pi;
        float_type tilde_phi_r_plus = (2.*gomv_rp - gomv_black_hole_spin_a*gomv_angular_momentum/gomv_energy)/(gomv_r2 - gomv_rp)/(gomv_r3 - gomv_rp)*2.*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hp);
        float_type tilde_phi_r_minus = (2.*gomv_rm - gomv_black_hole_spin_a*gomv_angular_momentum/gomv_energy)/(gomv_r2 - gomv_rm)/(gomv_r3 - gomv_rm)*2.*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hm);
        float_type tilde_phi_r_pi = -2. * gomv_black_hole_spin_a*gomv_energy*(gomv_r2 - gomv_r3)/(gomv_rp-gomv_rm)/mp::sqrt((1.-gomv_energy*gomv_energy)*(gomv_r1 - gomv_r3)*(gomv_r2-gomv_r4))*(tilde_phi_r_plus -tilde_phi_r_minus);

        float_type u = boost::math::ellint_1(mp::sqrt(gomv_kr))*qr/constants::pi;
        
        float_type sn = boost::math::jacobi_sn(mp::sqrt(gomv_kr), u);
        if(sn >= 0)
        {
            float_type cn = boost::math::jacobi_cn(mp::sqrt(gomv_kr), u);
            y = mp::acos(cn);
        }
        else
        {
            float_type cn = boost::math::jacobi_cn(mp::sqrt(gomv_kr), u);
            y = 2.*constants::pi - mp::acos(cn);   
        }
        

        tilde_phi_r_plus = (2.*gomv_rp - gomv_black_hole_spin_a*gomv_angular_momentum/gomv_energy)/(gomv_r2 - gomv_rp)/(gomv_r3 - gomv_rp)*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hp, y);
        tilde_phi_r_minus = (2.*gomv_rm - gomv_black_hole_spin_a*gomv_angular_momentum/gomv_energy)/(gomv_r2 - gomv_rm)/(gomv_r3 - gomv_rm)*boost::math::ellint_3(mp::sqrt(gomv_kr), gomv_hm, y);
        float_type tilde_phi_r_y = -2. * gomv_black_hole_spin_a*gomv_energy*(gomv_r2 - gomv_r3)/(gomv_rp-gomv_rm)/mp::sqrt((1.-gomv_energy*gomv_energy)*(gomv_r1 - gomv_r3)*(gomv_r2-gomv_r4))*(tilde_phi_r_plus -tilde_phi_r_minus);
        
        phi_r = tilde_phi_r_y - tilde_phi_r_pi/2./constants::pi*qr;
    }
}

void GeodesicOrbitalMotion::gomm_compute_phiz_function(const float_type& qz, float_type& phi_z)
{
    if(mp::abs(qz - constants::pi) < constants::precision_q || mp::abs(qz - constants::pi*0.5) < constants::precision_q || mp::abs(qz - 0.) < constants::precision_q)
    {
        phi_z = 0;
    }
    else
    {
        float_type tilde_phi_z_pi = - gomv_angular_momentum / gomv_z2 * 2. * boost::math::ellint_3(mp::sqrt(gomv_kz), gomv_z1*gomv_z1);
    
        float_type u = 2.*boost::math::ellint_1(mp::sqrt(gomv_kz))*(qz+constants::pi/2.)/constants::pi;
    
        float_type sn = boost::math::jacobi_sn(mp::sqrt(gomv_kz), u);
        float_type y;
        if(sn >= 0)
        {
            float_type cn = boost::math::jacobi_cn(mp::sqrt(gomv_kz), u);
            y = mp::acos(cn);
        }
        else
        {
            float_type cn = boost::math::jacobi_cn(mp::sqrt(gomv_kz), u);
            y = 2.*constants::pi - mp::acos(cn);   
        }

        float_type tilde_phi_z_y = - gomv_angular_momentum / gomv_z2 * boost::math::ellint_3(mp::sqrt(gomv_kz), gomv_z1*gomv_z1, y);
    
        phi_z = tilde_phi_z_y - tilde_phi_z_pi /constants::pi*(qz+constants::pi/2.);
        
        phi_z = - phi_z;
    }

}

void GeodesicOrbitalMotion::gomm_compute_phi_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& phi)
{
    float_type phi_r;
    gomm_compute_phir_function(qr, phi_r);
    
    float_type phi_z;
    gomm_compute_phiz_function(qz, phi_z);

    phi = qphi + phi_r + phi_z; 
}

void GeodesicOrbitalMotion::gomm_compute_four_velocities(const float_type& qr, const float_type& qz, float_type& ur, float_type& uz)
{
    float_type ellint_1_z = boost::math::ellint_1(mp::sqrt(gomv_kz));
    float_type arg_z = 2*ellint_1_z*qz/constants::pi;
    float_type dn_z = boost::math::jacobi_dn(mp::sqrt(gomv_kz), arg_z);
    float_type cn_z = boost::math::jacobi_cn(mp::sqrt(gomv_kz), arg_z);

    float_type ellint_1_r = boost::math::ellint_1(mp::sqrt(gomv_kr));
    float_type arg_r = ellint_1_r*qr/constants::pi;
    float_type dn_r = boost::math::jacobi_dn(mp::sqrt(gomv_kr), arg_r);
    float_type cn_r = boost::math::jacobi_cn(mp::sqrt(gomv_kr), arg_r);
    float_type sn_r = boost::math::jacobi_sn(mp::sqrt(gomv_kr), arg_r);
    


    float_type rprime = (2*(gomv_r1 - gomv_r2)*(gomv_r1 - gomv_r3)*(gomv_r2 - gomv_r3)*ellint_1_r*cn_r*dn_r*sn_r)/
    (constants::pi*mp::pow(-gomv_r1 + gomv_r3 + (gomv_r1 - gomv_r2)*mp::pow(sn_r,2),2));

    float_type zprime = (2*gomv_z1*ellint_1_z*cn_z*dn_z)/constants::pi;    

    float_type var_r = qr;
    float_type var_z = constants::pi/2. + (qr*gomv_Upsilon_z)/gomv_Upsilon_r;
    float_type z;
    float_type r;
    gomm_get_z_function(var_z, z);
    gomm_get_r_function(var_r, r);
    complex_type rho_r = -(1/(r - constants::I*gomv_black_hole_spin_a*z));
    complex_type rhobar_r = -(1/(r + constants::I*gomv_black_hole_spin_a*z));
    complex_type sigma_r_complex = 1/(rho_r*rhobar_r);
    float_type sigma_r = sigma_r_complex.real();

    float_type ur0 = (gomv_Upsilon_r*rprime)/sigma_r;

    var_r = ((-0.5*constants::pi + qz)*gomv_Upsilon_r)/gomv_Upsilon_z;
    var_z = qz;
    gomm_get_z_function(var_z, z);
    gomm_get_r_function(var_r, r);
    complex_type rho_z = -(1/(r - constants::I*gomv_black_hole_spin_a*z));
    complex_type rhobar_z = -(1/(r + constants::I*gomv_black_hole_spin_a*z));
    complex_type sigma_z_complex = 1/(rho_z*rhobar_z);
    float_type sigma_z = sigma_z_complex.real();

    float_type uz0 = -((gomv_Upsilon_z*zprime)/(mp::sqrt(1 - mp::pow(z,2))*sigma_z));

    var_r = (qz - constants::pi*0.5)*gomv_Upsilon_r/gomv_Upsilon_z;
    var_z = qz;
    gomm_get_r_function(var_r, r);
    gomm_get_z_function(var_z, z);

    uz = (r*r + gomv_black_hole_spin_a*gomv_black_hole_spin_a*z*z)*uz0;

    var_r = qr;
    var_z = qr*(gomv_Upsilon_z/gomv_Upsilon_r) + constants::pi*0.5;
    gomm_get_r_function(var_r, r);
    gomm_get_z_function(var_z, z);

    ur = (r*r + gomv_black_hole_spin_a*gomv_black_hole_spin_a*z*z)*ur0;
}

void GeodesicOrbitalMotion::gomm_get_r_function(const float_type& qr, float_type& r)
{
    gomm_compute_r_function(qr, r);
}

void GeodesicOrbitalMotion::gomm_get_z_function(const float_type& qz, float_type& z)
{
    gomm_compute_z_function(qz, z);
}

void GeodesicOrbitalMotion::gomm_get_tr_function(const float_type& qr, float_type& tr)
{
    gomm_compute_tr_function(qr, tr);
}

void GeodesicOrbitalMotion::gomm_get_tz_function(const float_type& qz, float_type& tz)
{
    gomm_compute_tz_function(qz, tz);
}

void GeodesicOrbitalMotion::gomm_get_t_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& t)
{
    gomm_compute_t_function(qr, qz, qphi, t);
}

void GeodesicOrbitalMotion::gomm_get_phir_function(const float_type& qr, float_type& phi_r)
{
    gomm_compute_phir_function(qr, phi_r);
}

void GeodesicOrbitalMotion::gomm_get_phiz_function(const float_type& qz, float_type& phi_z)
{
    gomm_compute_phiz_function(qz, phi_z);
}

void GeodesicOrbitalMotion::gomm_get_phi_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& phi)
{
    gomm_compute_phi_function(qr, qz, qphi, phi);
}

void GeodesicOrbitalMotion::gomm_get_four_velocities(const float_type& qr, const float_type& qz, float_type& ur, float_type& uz)
{
    gomm_compute_four_velocities(qr, qz, ur, uz);
}

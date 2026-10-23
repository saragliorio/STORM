#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include "config.h"
#include "special_functions.h"

#ifndef GeodesicOrbitalMotion_H
#define GeodesicOrbitalMotion_H

typedef std::vector<float_type> state_type;

// Analytic geodesics: M. van de Meent, Class. Quant. Grav. 37, 145007 (2020).
class GeodesicOrbitalMotion
{
    private:
        float_type gomv_mass_black_hole_M;
        float_type gomv_mass_test_particle_m;
        float_type gomv_black_hole_spin_a;
        float_type gomv_semilatus_rectum_p;
        float_type gomv_eccentricity_e;
        float_type gomv_inclination_angle;

        int gomv_l;
        int gomv_m;
        int gomv_k;
        int gomv_n;

        bool gomv_prograde_retrograde_orbit;

        float_type gomv_theta_min;
        float_type gomv_z1;
        float_type gomv_z2;
        float_type gomv_r1;
        float_type gomv_r2;
        float_type gomv_r3;
        float_type gomv_r4;
        float_type gomv_energy;
        float_type gomv_angular_momentum;
        float_type gomv_carter_constant;

        bool gomv_radial_constants_ready = false;
        float_type gomv_radial_K, gomv_radial_tr_pi, gomv_radial_phir_pi;
        float_type gomv_kr;
        float_type gomv_kz;
        float_type gomv_rp;
        float_type gomv_rm;
        float_type gomv_hr;
        float_type gomv_hp;
        float_type gomv_hm;
        float_type gomv_Upsilon_r;
        float_type gomv_Upsilon_z;
        float_type gomv_Upsilon_phi;
        float_type gomv_Gamma;
        float_type gomv_Omega_r;
        float_type gomv_Omega_theta;
        float_type gomv_Omega_phi;
        float_type gomv_beta;
        float_type gomv_angular_frequency_omega;
        void gomm_compute_delta_function(const float_type& r, float_type& delta);
        void gomm_compute_constants_of_motion();
        void gomm_compute_frequencies_mino_time();
        void gomm_compute_frequencies_standard_time();
        void gomm_compute_angular_frequency_omega();
        void gomm_compute_r_function(const float_type& qr, float_type& r);
        void gomm_compute_z_function(const float_type& qz, float_type& z);
        float_type gomm_compute_tr_pi();
        float_type gomm_compute_phir_pi();
        void gomm_compute_tr_function(const float_type& qr, float_type& tr);
        void gomm_compute_tz_function(const float_type& qz, float_type& tz);
        void gomm_compute_t_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& t);
        void gomm_compute_phir_function(const float_type& qr, float_type& phi_r);
        void gomm_compute_phiz_function(const float_type& qz, float_type& phi_z);
        void gomm_compute_phi_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& phi);
        void gomm_compute_four_velocities(const float_type& qr, const float_type& qz, float_type& ur, float_type& uz);

    public:
        GeodesicOrbitalMotion(float_type M, float_type m, float_type a, float_type p, float_type e, float_type theta_inc, int ll, int mm, int kk, int nn);

        void gomm_get_frequencies_mino_time(float_type& Ur, float_type& Utheta, float_type& Uphi, float_type& G);
        void gomm_get_frequencies_standard_time(float_type& Or, float_type& Otheta, float_type& Ophi);
        void gomm_get_constants_of_motion(float_type& en, float_type& ang, float_type& crt);
        void gomm_get_angular_frequency_omega(float_type& omega);
        void gomm_cache_radial_constants();
        void gomm_get_radial_source(const float_type& qr, float_type& r,
            float_type& tr, float_type& phi_r);

        void gomm_get_r_function(const float_type& qr, float_type& r);
        void gomm_get_z_function(const float_type& qz, float_type& z);
        void gomm_get_tr_function(const float_type& qr, float_type& tr);
        void gomm_get_tz_function(const float_type& qz, float_type& tz);
        void gomm_get_t_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& t);
        void gomm_get_phir_function(const float_type& qr, float_type& phi_r);
        void gomm_get_phiz_function(const float_type& qz, float_type& phi_z);
        void gomm_get_phi_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& phi);
        void gomm_get_four_velocities(const float_type& qr, const float_type& qz, float_type& ur, float_type& uz);
};

#endif
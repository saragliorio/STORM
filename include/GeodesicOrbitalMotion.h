#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include "config.h"
#include "special_functions.h"

#ifndef GeodesicOrbitalMotion_H
#define GeodesicOrbitalMotion_H

// state vector type used for orbital evolution
typedef std::vector<float_type> state_type;

//REFERENCE: M. van de Meent, Class. Quant. Grav. 37, 145007 (2020) (for analytical solutions of the geodesic motion)

// Class that computes the geodesic motion of a test particle orbiting a rotating black hole.
// This class evaluates orbital properties of a compact object moving in the Kerr spacetime.
// It provides routines to compute constants of motion, orbital frequencies, trajectory
// functions, and four-velocity components of the orbiting particle.
class GeodesicOrbitalMotion
{
    private: 

        // --- Physical parameters of the system ---

        float_type gomv_mass_black_hole_M;       // mass of the primary black hole
        float_type gomv_mass_test_particle_m;    // mass of the orbiting compact object
        float_type gomv_black_hole_spin_a;       // dimensionless spin parameter of the black hole

        // --- Orbital parameters ---

        float_type gomv_semilatus_rectum_p;      // semilatus rectum of the orbit
        float_type gomv_eccentricity_e;          // orbital eccentricity
        float_type gomv_inclination_angle;       // orbital inclination angle

        int gomv_l;                              // spherical harmonic index l
        int gomv_m;                              // azimuthal harmonic index m
        int gomv_k;                              // polar harmonic index k
        int gomv_n;                              // radial harmonic index n

        bool gomv_prograde_retrograde_orbit;     // flag indicating prograde or retrograde motion

        // --- Angular motion parameters ---

        float_type gomv_theta_min;               // minimum polar angle reached during the motion
        float_type gomv_z1;                      // auxiliary parameter related to polar motion
        float_type gomv_z2;                      // auxiliary parameter related to polar motion

        // --- Radial turning points of the orbit ---

        float_type gomv_r1;                      // first radial root
        float_type gomv_r2;                      // second radial root
        float_type gomv_r3;                      // third radial root
        float_type gomv_r4;                      // fourth radial root
       
        // --- Constants of motion ---

        float_type gomv_energy;                  // conserved orbital energy
        float_type gomv_angular_momentum;        // conserved azimuthal angular momentum
        float_type gomv_carter_constant;         // Carter constant of motion
        
        // --- Auxiliary quantities used in frequency calculations ---

        float_type gomv_kr;                      // radial elliptic modulus
        float_type gomv_kz;                      // polar elliptic modulus
        float_type gomv_rp;                      // radial parameter related to turning points
        float_type gomv_rm;                      // radial parameter related to turning points
        
        float_type gomv_hr;                      // auxiliary radial quantity
        float_type gomv_hp;                      // auxiliary quantity related to rp
        float_type gomv_hm;                      // auxiliary quantity related to rm

        // --- Fundamental frequencies in Mino time ---

        float_type gomv_Upsilon_r;               // radial frequency in Mino time
        float_type gomv_Upsilon_z;               // polar frequency in Mino time
        float_type gomv_Upsilon_phi;             // azimuthal frequency in Mino time
        
        // --- Frequencies in coordinate time ---

        float_type gomv_Gamma;                   // conversion factor between Mino time and coordinate time frequencies
        float_type gomv_Omega_r;                 // radial frequency in coordinate time        
        float_type gomv_Omega_theta;             // polar frequency in coordinate time
        float_type gomv_Omega_phi;               // azimuthal frequency in coordinate time
        
        // --- Additional parameters ---

        float_type gomv_beta;                    // auxiliary parameter related to polar motion
        float_type gomv_angular_frequency_omega; // combined angular frequency 
    
        void gomm_compute_delta_function(const float_type& r, float_type& delta);
        // Computes the Kerr metric delta function at a given radial coordinate.

        void gomm_compute_constants_of_motion();
        // Computes the constants of motion of the geodesic orbit.

        void gomm_compute_frequencies_mino_time();
        // Computes the fundamental orbital frequencies in Mino time.
       
        void gomm_compute_frequencies_standard_time();
        // Converts Mino-time frequencies to coordinate-time frequencies.

        void gomm_compute_angular_frequency_omega();
        //Computes the angular frequency of the orbit.

        // --- Orbital coordinate functions ---

        void gomm_compute_r_function(const float_type& qr, float_type& r);
        // radial coordinate as a function of radial phase

        void gomm_compute_z_function(const float_type& qz, float_type& z);
        // polar coordinate as a function of polar phase

        void gomm_compute_tr_function(const float_type& qr, float_type& tr);
        // radial contribution to coordinate time

        void gomm_compute_tz_function(const float_type& qz, float_type& tz);
        // polar contribution to coordinate time

        void gomm_compute_t_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& t);
        // full coordinate time as a function of orbital phases

        void gomm_compute_phir_function(const float_type& qr, float_type& phi_r);
        // radial contribution to the azimuthal angle

        void gomm_compute_phiz_function(const float_type& qz, float_type& phi_z);
        // polar contribution to the azimuthal angle

        void gomm_compute_phi_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& phi);
        // total azimuthal angle as a function of orbital phases

        void gomm_compute_four_velocities(const float_type& qr, const float_type& qz, float_type& ur, float_type& uz);
        // Computes the radial and polar components of the four-velocity.
    
        public:

        GeodesicOrbitalMotion(float_type M, float_type m, float_type a, float_type p, float_type e, float_type theta_inc, int ll, int mm, int kk, int nn);
        // Constructor. Initializes the orbital system with physical and orbital parameters.
        
        void gomm_get_frequencies_mino_time(float_type& Ur, float_type& Utheta, float_type& Uphi, float_type& G);
        // Returns the fundamental frequencies in Mino time.
        
        void gomm_get_frequencies_standard_time(float_type& Or, float_type& Otheta, float_type& Ophi);
        // Returns the orbital frequencies in coordinate time.
    
        void gomm_get_constants_of_motion(float_type& en, float_type& ang, float_type& crt);
        // Returns the constants of motion of the orbit.

        void gomm_get_angular_frequency_omega(float_type& omega);
        // Returns the angular frequency of the orbit.

        // --- Public interface to orbital coordinate functions ---

        void gomm_get_r_function(const float_type& qr, float_type& r);
        void gomm_get_z_function(const float_type& qz, float_type& z);
        
        void gomm_get_tr_function(const float_type& qr, float_type& tr);
        void gomm_get_tz_function(const float_type& qz, float_type& tz);
        void gomm_get_t_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& t);
        
        void gomm_get_phir_function(const float_type& qr, float_type& phi_r);
        void gomm_get_phiz_function(const float_type& qz, float_type& phi_z);
        void gomm_get_phi_function(const float_type& qr, const float_type& qz, const float_type& qphi, float_type& phi);

        void gomm_get_four_velocities(const float_type& qr, const float_type& qz, float_type& ur, float_type& uz);
        // Returns the radial and polar components of the four-velocity.
        
        // destructor
        ~GeodesicOrbitalMotion() 
        {
        }
};

#endif
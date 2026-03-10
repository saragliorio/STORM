#pragma once

#include <iostream>
#include <cmath>
#include "config.h"
#include "special_functions.h"
#include "RadialHomogeneousSolution.h"
#include "SpinWeightedSpheroidalHarmonics.h"
#include "GeodesicOrbitalMotion.h"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <memory>  

#ifndef SCALARFLUX_H
#define SCALARFLUX_H

// REFERENCES: our paper 

  

  // Class that computes scalar energy fluxes for a particle orbiting a Kerr black hole.

  // The class handles:
  //   - orbital motion using GeodesicOrbitalMotion
  //   - angular functions using SpinWeightedSpheroidalHarmonics
  //   - radial solutions using RadialHomogeneousSolution
  //   - computation of energy flux at the horizon and at infinity

  // Inputs include:
  //   - spin weight of the field
  //   - black hole parameters (mass, spin)
  //   - orbit parameters (semilatus rectum, eccentricity, inclination)
  //   - mode numbers (l, m, k, n)

  class ScalarFluxMode
{
    private: 
        int sfmv_spin_weight_s;                  // spin weight s
        float_type sfmv_mass_black_hole_M;       // mass of the primary black hole
        float_type sfmv_mass_test_particle_m;    // mass of the orbiting compact object
        float_type sfmv_black_hole_spin_a;       // black hole spin parameter
        float_type sfmv_semilatus_rectum_p;      // semilatus rectum of the orbit
        float_type sfmv_eccentricity_e;          // orbital eccentricity
        float_type sfmv_inclination_angle;       // orbital inclination angle
        int sfmv_l;                              // spherical harmonic index l
        int sfmv_m;                              // azimuthal harmonic index m
        int sfmv_k;                              // polar harmonic index k
        int sfmv_n;                              // radial harmonic index n

        std::unique_ptr<GeodesicOrbitalMotion> sfmv_orbmot; // orbit object for trajectory and frequencies
        float_type sfmv_energy;                  // energy of the particle
        float_type sfmv_angular_momentum;        // angular momentum of the particle
        float_type sfmv_carter_constant;         // Carter constant of the orbit
        float_type sfmv_Upsilon_r;               // radial Mino frequency
        float_type sfmv_Upsilon_theta;           // polar Mino frequency
        float_type sfmv_Upsilon_phi;             // azimuthal Mino frequency
        float_type sfmv_Gamma;                   // scaling factor for Mino to coordinate time
        float_type sfmv_Omega_r;                 // radial frequency in coordinate time
        float_type sfmv_Omega_theta;             // polar frequency in coordinate time
        float_type sfmv_Omega_phi;               // azimuthal frequency in coordinate time
        float_type sfmv_angular_frequency_omega; // combined orbital angular frequency

        std::unique_ptr<SpinWeightedSpheroidalHarmonics> sfmv_swh;   // spin-weighted spheroidal harmonic object
        std::unique_ptr<RadialHomogeneousSolution> sfmv_homsol;      // radial homogeneous solution object

        float_type sfmv_lambda;                                     // angular eigenvalue

        complex_type sfmv_W;                  // Wronskian for flux computation

        float_type sfmv_scalar_energy_flux_infinity; // scalar flux at infinity
        float_type sfmv_scalar_energy_flux_horizon;  // scalar flux at the horizon
        
        float_type sfmv_r_min;                // minimum radial coordinate of the orbit
        float_type sfmv_r_max;                // maximum radial coordinate of the orbit
        
        // --- Internal integration routines for flux computation ---
        void sfmm_compute_I_plus(const float_type& qr, complex_type& I_plus);    // radial integral for outgoing mode
        void sfmm_compute_I_minus(const float_type& qr, complex_type& I_minus);  // radial integral for ingoing mode
        void sfmm_compute_I_z(const float_type& qz, float_type& Iz);             // polar integral over theta
        
        void sfmm_compute_Iz1(float_type& Iz1);                                  // subcomponent of polar integral
        void sfmm_compute_Iz2(float_type& Iz2);                                  // subcomponent of polar integral
        void sfmm_compute_Ir1_plus(complex_type& Ir1_plus);                      // radial subcomponent for outgoing flux
        void sfmm_compute_Ir2_plus(complex_type& Ir2_plus);                      // radial subcomponent for outgoing flux
        void sfmm_compute_Ir1_minus(complex_type& Ir1_minus);                    // radial subcomponent for ingoing flux
        void sfmm_compute_Ir2_minus(complex_type& Ir2_minus);                    // radial subcomponent for ingoing flux

        // --- Energy flux computation for generic orbits ---
        void sfmm_compute_energy_flux_horizon_infinity_generic_orbit();          // compute both horizon and infinity fluxes on generic orbits
        void sfmm_compute_energy_flux_infinity_generic_orbit();                  // compute flux at infinity on generic orbits
        void sfmm_compute_energy_flux_horizon_generic_orbit();                   // compute flux at horizon on generic orbits 

        // --- Energy flux computation for circular orbits ---
        void sfmm_compute_energy_flux_horizon_infinity_circular_orbit();        // compute both horizon and infinity fluxes on circular orbits
        void sfmm_compute_energy_flux_infinity_circular_orbit();                // compute flux at infinity on circular orbits
        void sfmm_compute_energy_flux_horizon_circular_orbit();                 // compute flux at horizon on circular orbits
        
        // --- Energy flux computation for eccentric orbits ---
        void sfmm_compute_energy_flux_horizon_infinity_eccentric_orbit();       // compute both horizon and infinity fluxes on eccnetric orbits
        void sfmm_compute_energy_flux_horizon_eccentric_orbit();                 // compute flux at horizon on eccnetric orbits
        void sfmm_compute_energy_flux_infinity_eccentric_orbit();                // compute flux at infinity on eccnetric orbits 

    public:
        // constructor: initialize flux computation for a specific mode
        ScalarFluxMode(int s, float_type M, float_type m, float_type a, float_type p, float_type e, float_type theta_inc, int ll, int mm, int kk, int nn);
        
        // --- Public interface for generic orbits ---
        void sfmm_get_energy_flux_horizon_infinity_generic_orbit(float_type& flux_horizon, float_type& flux_infinity);
        void sfmm_get_energy_flux_infinity_generic_orbit(float_type& flux_infinity);
        void sfmm_get_energy_flux_horizon_generic_orbit(float_type& flux_horizon);
        
        // --- Public interface for eccentric orbits ---
        void sfmm_get_energy_flux_horizon_infinity_eccentric_orbit(float_type& flux_horizon, float_type& flux_infinity);
        void sfmm_get_energy_flux_horizon_eccentric_orbit(float_type& flux_horizon);
        void sfmm_get_energy_flux_infinity_eccentric_orbit(float_type& flux_infinity);

        // --- Public interface for circular orbits ---
        void sfmm_get_energy_flux_horizon_infinity_circular_orbit(float_type& flux_horizon, float_type& flux_infinity);
        void sfmm_get_energy_flux_infinity_circular_orbit(float_type& flux_infinity);
        void sfmm_get_energy_flux_horizon_circular_orbit(float_type& flux_horizon);

        // destructor
        ~ScalarFluxMode() 
        {
        }
    };

#endif
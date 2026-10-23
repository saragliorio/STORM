#pragma once

#include <iostream>
#include <cmath>
#include "config.h"
#include "special_functions.h"
#include "ScalarRadialBackend.h"
#include "SpinWeightedSpheroidalHarmonics.h"
#include "GeodesicOrbitalMotion.h"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <memory>
#include <map>
#include <utility>

#ifndef SCALARFLUX_H
#define SCALARFLUX_H

struct ScalarFluxOptions {
    ScalarRadialBackendKind backend = ScalarRadialBackendKind::Original;
    CompactifiedCharacteristicRadialSolution::Options endpoint_options;
    float_type quadrature_tolerance = constants::precision_trapezoidal_quadrature;

    // Evaluate nonzero frequencies below the legacy zero-flux threshold.
    bool bypass_small_frequency_guard = false;
};

class ScalarFluxMode
{
    private:
        int sfmv_spin_weight_s;
        float_type sfmv_mass_black_hole_M;
        float_type sfmv_mass_test_particle_m;
        float_type sfmv_black_hole_spin_a;
        float_type sfmv_semilatus_rectum_p;
        float_type sfmv_eccentricity_e;
        float_type sfmv_inclination_angle;
        int sfmv_l;
        int sfmv_m;
        int sfmv_k;
        int sfmv_n;

        std::unique_ptr<GeodesicOrbitalMotion> sfmv_orbmot;
        float_type sfmv_Gamma;
        float_type sfmv_angular_frequency_omega;

        std::unique_ptr<SpinWeightedSpheroidalHarmonics> sfmv_swh;
        std::unique_ptr<ScalarRadialBackend> sfmv_homsol;
        ScalarFluxOptions sfmv_options;
        bool sfmv_explicit_options = false;
        complex_type sfmv_amplitude_horizon = 0;
        complex_type sfmv_amplitude_infinity = 0;
        bool sfmv_amplitude_horizon_ready = false;
        bool sfmv_amplitude_infinity_ready = false;
        std::map<float_type, std::pair<float_type, float_type>> sfmv_radial_source_samples;

        float_type sfmv_lambda;

        complex_type sfmv_W;

        float_type sfmv_scalar_energy_flux_infinity;
        float_type sfmv_scalar_energy_flux_horizon;

        float_type sfmv_r_min;
        float_type sfmv_r_max;

        void sfmm_get_radial_source_sample(const float_type& qr, float_type& r, float_type& phase);
        void sfmm_compute_I_plus(const float_type& qr, complex_type& I_plus, float_type* radius = nullptr);
        void sfmm_compute_I_minus(const float_type& qr, complex_type& I_minus, float_type* radius = nullptr);
        void sfmm_compute_I_z(const float_type& qz, float_type& Iz);

        void sfmm_compute_Iz1(float_type& Iz1);
        void sfmm_compute_Iz2(float_type& Iz2);
        void sfmm_compute_Ir1_plus(complex_type& Ir1_plus);
        void sfmm_compute_Ir2_plus(complex_type& Ir2_plus);
        void sfmm_compute_Ir1_minus(complex_type& Ir1_minus);
        void sfmm_compute_Ir2_minus(complex_type& Ir2_minus);

        void sfmm_compute_energy_flux_horizon_infinity_generic_orbit();
        void sfmm_compute_energy_flux_infinity_generic_orbit();
        void sfmm_compute_energy_flux_horizon_generic_orbit();

        void sfmm_compute_energy_flux_horizon_infinity_circular_orbit();
        void sfmm_compute_energy_flux_infinity_circular_orbit();
        void sfmm_compute_energy_flux_horizon_circular_orbit();

        void sfmm_compute_energy_flux_horizon_infinity_eccentric_orbit();
        void sfmm_compute_energy_flux_horizon_eccentric_orbit();
        void sfmm_compute_energy_flux_infinity_eccentric_orbit();

        bool sfmm_frequency_guard() const;
        float_type sfmm_quadrature_tolerance() const;
        void sfmm_record_amplitude(bool horizon, complex_type& amplitude);
        void sfmm_zero_flux_and_amplitudes();

    public:
        ScalarFluxMode(int s, float_type M, float_type m, float_type a, float_type p, float_type e, float_type theta_inc, int ll, int mm, int kk, int nn);
        ScalarFluxMode(int s, float_type M, float_type particle_mass, float_type a,
            float_type p, float_type e, float_type theta_inc, int ll,
            int azimuthal_m, int kk, int nn, const ScalarFluxOptions& options);

        // Available after both fluxes have been computed.
        void sfmm_get_complex_source_amplitudes(complex_type& A_H,
            complex_type& A_inf) const;
        void sfmm_rescale_radial_basis(const complex_type& in_scale,
            const complex_type& up_scale);
        float_type sfmm_angular_frequency() const { return sfmv_angular_frequency_omega; }
        float_type sfmm_lambda() const { return sfmv_lambda; }
        const ScalarRadialBackend& sfmm_radial_backend() const { return *sfmv_homsol; }

        void sfmm_get_energy_flux_horizon_infinity_generic_orbit(float_type& flux_horizon, float_type& flux_infinity);
        void sfmm_get_energy_flux_infinity_generic_orbit(float_type& flux_infinity);
        void sfmm_get_energy_flux_horizon_generic_orbit(float_type& flux_horizon);

        void sfmm_get_energy_flux_horizon_infinity_eccentric_orbit(float_type& flux_horizon, float_type& flux_infinity);
        void sfmm_get_energy_flux_horizon_eccentric_orbit(float_type& flux_horizon);
        void sfmm_get_energy_flux_infinity_eccentric_orbit(float_type& flux_infinity);

        void sfmm_get_energy_flux_horizon_infinity_circular_orbit(float_type& flux_horizon, float_type& flux_infinity);
        void sfmm_get_energy_flux_infinity_circular_orbit(float_type& flux_infinity);
        void sfmm_get_energy_flux_horizon_circular_orbit(float_type& flux_horizon);
};

#endif

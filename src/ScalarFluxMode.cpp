#include "ScalarFluxMode.h"
#include <stdexcept>

namespace {
template<class F>
auto checked_trapezoidal(F f, const float_type& a, const float_type& b,
    const float_type& tolerance) -> decltype(f(a))
{
    float_type error = 0, L1 = 0;
    auto integral = boost::math::quadrature::trapezoidal(f, a, b,
        tolerance, 12, &error, &L1);
    if (!(error <= tolerance*L1))
        throw std::runtime_error("Scalar flux quadrature did not converge");
    return integral;
}
}

ScalarFluxMode::ScalarFluxMode(int s, float_type M, float_type m, float_type a,
    float_type p, float_type e, float_type theta_inc, int ll, int mm, int kk, int nn)
    : ScalarFluxMode(s,M,m,a,p,e,theta_inc,ll,mm,kk,nn,ScalarFluxOptions{}) {
    sfmv_explicit_options = false;
}

ScalarFluxMode::ScalarFluxMode(int s, float_type M, float_type m, float_type a,
    float_type p, float_type e, float_type theta_inc, int ll, int mm, int kk, int nn,
    const ScalarFluxOptions& options) :
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
        sfmv_n(nn),
        sfmv_options(options),
        sfmv_explicit_options(true)
    {

        if (sfmv_options.backend != ScalarRadialBackendKind::Original &&
            (!mp::isfinite(M) || !mp::isfinite(m) || !mp::isfinite(a) ||
             !mp::isfinite(p) || !mp::isfinite(e) || !mp::isfinite(theta_inc)))
            throw std::invalid_argument("Compactified scalar source requires finite masses, spin and orbital parameters");
        if (sfmv_options.backend != ScalarRadialBackendKind::Original &&
            (s != 0 || M != 1))
            throw std::invalid_argument("Compactified scalar flux supports only s=0 and M=1");
        if (sfmv_options.backend != ScalarRadialBackendKind::Original && e >= 1)
            throw std::invalid_argument("Compactified scalar source requires a bound orbit with e<1");
        if (sfmv_options.backend != ScalarRadialBackendKind::Original &&
            !mp::isfinite(sfmv_options.quadrature_tolerance))
            throw std::invalid_argument("Source quadrature tolerance must be finite and positive");
        if (sfmv_options.quadrature_tolerance <= 0)
            throw std::invalid_argument("Source quadrature tolerance must be positive");

        if (mp::abs(sfmv_black_hole_spin_a) > 0.999)
        {
            std::cerr << "Error in ScalarFluxMode class: black hole spin parameter a must be <= 0.999" << std::endl;
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

        float_type Upsilon_r, Upsilon_theta, Upsilon_phi;
        sfmv_orbmot->gomm_get_frequencies_mino_time(Upsilon_r, Upsilon_theta, Upsilon_phi, sfmv_Gamma);

        float_type Omega_r, Omega_theta, Omega_phi;
        sfmv_orbmot->gomm_get_frequencies_standard_time(Omega_r, Omega_theta, Omega_phi);

        sfmv_orbmot->gomm_get_angular_frequency_omega(sfmv_angular_frequency_omega);
        if (sfmv_options.backend != ScalarRadialBackendKind::Original && e > 0)
            sfmv_orbmot->gomm_cache_radial_constants();

        sfmv_swh = std::make_unique<SpinWeightedSpheroidalHarmonics>(sfmv_spin_weight_s, sfmv_l, sfmv_m, sfmv_black_hole_spin_a*sfmv_angular_frequency_omega);
        sfmv_swh->swshm_get_lambda(sfmv_lambda);

        sfmv_r_min = sfmv_semilatus_rectum_p/(1 + sfmv_eccentricity_e);
        sfmv_r_max = sfmv_semilatus_rectum_p/(1 - sfmv_eccentricity_e);

        sfmv_homsol = std::make_unique<ScalarRadialBackend>(sfmv_spin_weight_s,
            sfmv_mass_black_hole_M,sfmv_black_hole_spin_a,sfmv_mass_test_particle_m,
            sfmv_lambda,sfmv_angular_frequency_omega,sfmv_l,sfmv_m,
            sfmv_r_min,sfmv_r_max,sfmv_options.backend,sfmv_options.endpoint_options,
            sfmv_options.bypass_small_frequency_guard);
        sfmv_W = sfmv_homsol->wronskian((sfmv_r_min+sfmv_r_max)/2);

    }

bool ScalarFluxMode::sfmm_frequency_guard() const {
    return sfmv_angular_frequency_omega == 0 ||
        (!sfmv_options.bypass_small_frequency_guard &&
         mp::abs(sfmv_angular_frequency_omega) <= constants::precision_threshold);
}

float_type ScalarFluxMode::sfmm_quadrature_tolerance() const {
    return sfmv_explicit_options ? sfmv_options.quadrature_tolerance
        : constants::precision_trapezoidal_quadrature;
}

void ScalarFluxMode::sfmm_record_amplitude(bool horizon, complex_type& amplitude) {
    if (sfmv_options.backend != ScalarRadialBackendKind::Original)
        amplitude *= horizon ? sfmv_homsol->horizon_amplitude()
            : sfmv_homsol->infinity_amplitude();
    if (horizon) {
        sfmv_amplitude_horizon = amplitude;
        sfmv_amplitude_horizon_ready = true;
    } else {
        sfmv_amplitude_infinity = amplitude;
        sfmv_amplitude_infinity_ready = true;
    }
}

void ScalarFluxMode::sfmm_zero_flux_and_amplitudes() {
    sfmv_scalar_energy_flux_horizon = 0;
    sfmv_scalar_energy_flux_infinity = 0;
    sfmv_amplitude_horizon = 0;
    sfmv_amplitude_infinity = 0;
    sfmv_amplitude_horizon_ready = true;
    sfmv_amplitude_infinity_ready = true;
}

void ScalarFluxMode::sfmm_get_complex_source_amplitudes(complex_type& A_H,
    complex_type& A_inf) const {
    if (!sfmv_amplitude_horizon_ready || !sfmv_amplitude_infinity_ready)
        throw std::logic_error("Compute both fluxes before reading cached source amplitudes");
    A_H = sfmv_amplitude_horizon;
    A_inf = sfmv_amplitude_infinity;
}

void ScalarFluxMode::sfmm_rescale_radial_basis(const complex_type& in_scale,
    const complex_type& up_scale) {
    sfmv_homsol->rescale_basis(in_scale,up_scale);
    sfmv_W = sfmv_homsol->wronskian((sfmv_r_min+sfmv_r_max)/2);
    sfmv_amplitude_horizon_ready = false;
    sfmv_amplitude_infinity_ready = false;
}

void ScalarFluxMode::sfmm_get_radial_source_sample(const float_type& qr,
    float_type& r, float_type& phase)
{
    const auto found = sfmv_radial_source_samples.find(qr);
    if (found != sfmv_radial_source_samples.end()) {
        r = found->second.first;
        phase = found->second.second;
        return;
    }
    float_type tr, phi_r;
    sfmv_orbmot->gomm_get_radial_source(qr, r, tr, phi_r);
    phase = mp::cos(sfmv_n*qr + sfmv_angular_frequency_omega*tr-sfmv_m*phi_r);
    sfmv_radial_source_samples.emplace(qr, std::make_pair(r, phase));
}

void ScalarFluxMode::sfmm_compute_I_plus(const float_type& qr, complex_type& I_plus, float_type* radius)
{
    float_type tr, phi_r, r, phase;
    const bool reuse = sfmv_options.backend != ScalarRadialBackendKind::Original &&
        sfmv_eccentricity_e > 0 && sfmv_inclination_angle == 0;
    if (reuse) {
        sfmm_get_radial_source_sample(qr, r, phase);
    } else if (sfmv_options.backend == ScalarRadialBackendKind::Original) {
        sfmv_orbmot->gomm_get_tr_function(qr, tr);
        sfmv_orbmot->gomm_get_phir_function(qr, phi_r);
        sfmv_orbmot->gomm_get_r_function(qr, r);
    } else {
        sfmv_orbmot->gomm_get_radial_source(qr, r, tr, phi_r);
    }
    if (radius) *radius = r;
    complex_type R_in;

    sfmv_homsol->evaluate_in(r, R_in);

    if (reuse) I_plus = phase*R_in;
    else I_plus = mp::cos(sfmv_n*qr + sfmv_angular_frequency_omega*tr-sfmv_m*phi_r)*R_in;
}

void ScalarFluxMode::sfmm_compute_I_minus(const float_type& qr, complex_type& I_minus, float_type* radius)
{
    float_type tr, phi_r, r, phase;
    const bool reuse = sfmv_options.backend != ScalarRadialBackendKind::Original &&
        sfmv_eccentricity_e > 0 && sfmv_inclination_angle == 0;
    if (reuse) {
        sfmm_get_radial_source_sample(qr, r, phase);
    } else if (sfmv_options.backend == ScalarRadialBackendKind::Original) {
        sfmv_orbmot->gomm_get_tr_function(qr, tr);
        sfmv_orbmot->gomm_get_phir_function(qr, phi_r);
        sfmv_orbmot->gomm_get_r_function(qr, r);
    } else {
        sfmv_orbmot->gomm_get_radial_source(qr, r, tr, phi_r);
    }
    if (radius) *radius = r;

    complex_type R_up;

    sfmv_homsol->evaluate_up(r, R_up);

    if (reuse) I_minus = phase*R_up;
    else I_minus = mp::cos(sfmv_n*qr + sfmv_angular_frequency_omega*tr-sfmv_m*phi_r)*R_up;
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

    Iz1 = checked_trapezoidal(f_Iz1, qz_start, qz_end, sfmm_quadrature_tolerance());
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

    Iz2 = sfmv_black_hole_spin_a*sfmv_black_hole_spin_a*checked_trapezoidal(f_Iz2, qz_start, qz_end, sfmm_quadrature_tolerance());
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

    Ir1_plus = checked_trapezoidal(f_Ir1_plus, qr_start, qr_end, sfmm_quadrature_tolerance());
}

void ScalarFluxMode::sfmm_compute_Ir2_plus(complex_type& Ir2_plus)
{
    float_type qr_start = 0.;
    float_type qr_end = constants::pi;

    auto f_Ir2_plus = [this](float_type qr)
    {
        complex_type I_plus;
        float_type r;
        if (sfmv_options.backend == ScalarRadialBackendKind::Original) {
            sfmm_compute_I_plus(qr, I_plus);
            sfmv_orbmot->gomm_get_r_function(qr, r);
        } else {
            sfmm_compute_I_plus(qr, I_plus, &r);
        }

        return I_plus*r*r;
    };

    Ir2_plus = checked_trapezoidal(f_Ir2_plus, qr_start, qr_end, sfmm_quadrature_tolerance());
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

    Ir1_minus = checked_trapezoidal(f_Ir1_minus, qr_start, qr_end, sfmm_quadrature_tolerance());

}

void ScalarFluxMode::sfmm_compute_Ir2_minus(complex_type& Ir2_minus)
{
    float_type qr_start = 0.;
    float_type qr_end = constants::pi;

    auto f_Ir2_minus = [this](float_type qr)
    {
        complex_type I_minus;
        float_type r;
        if (sfmv_options.backend == ScalarRadialBackendKind::Original) {
            sfmm_compute_I_minus(qr, I_minus);
            sfmv_orbmot->gomm_get_r_function(qr, r);
        } else {
            sfmm_compute_I_minus(qr, I_minus, &r);
        }

        return I_minus*r*r;
    };

    Ir2_minus = checked_trapezoidal(f_Ir2_minus, qr_start, qr_end, sfmm_quadrature_tolerance());
}

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
    sfmm_record_amplitude(true, delta_phi_minus);

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
    sfmm_record_amplitude(false, delta_phi_plus);

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
    sfmm_record_amplitude(false, delta_phi_plus);
    complex_type delta_phi_minus = coeff * (Iz1*Ir2_minus + Iz2*Ir1_minus);
    sfmm_record_amplitude(true, delta_phi_minus);

    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;

    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );

}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_generic_orbit(float_type& flux_horizon)
{
    if(sfmm_frequency_guard())
    {
        sfmm_zero_flux_and_amplitudes();
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
    if(sfmm_frequency_guard())
    {
        sfmm_zero_flux_and_amplitudes();
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
    if(sfmm_frequency_guard())
    {
        sfmm_zero_flux_and_amplitudes();
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

void ScalarFluxMode::sfmm_compute_energy_flux_horizon_eccentric_orbit()
{
    complex_type Ir2_minus;
    sfmm_compute_Ir2_minus(Ir2_minus);

    complex_type coeff = (4./constants::pi/sfmv_Gamma/sfmv_W);

    complex_type S{0};
    sfmv_swh->swshm_get_S(0., 0., S);
    float_type S_real = S.real();

    complex_type delta_phi_minus = coeff * (Ir2_minus*S_real*constants::pi);
    sfmm_record_amplitude(true, delta_phi_minus);

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
    sfmm_record_amplitude(false, delta_phi_plus);

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
    sfmm_record_amplitude(false, delta_phi_plus);
    complex_type delta_phi_minus = coeff * (Ir2_minus*S_real*constants::pi);
    sfmm_record_amplitude(true, delta_phi_minus);

    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;

    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );
}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_eccentric_orbit(float_type& flux_horizon)
{
    if (sfmv_options.backend != ScalarRadialBackendKind::Original &&
        sfmv_angular_frequency_omega == 0) {
        sfmm_zero_flux_and_amplitudes();
        flux_horizon = 0;
        return;
    }

    sfmm_compute_energy_flux_horizon_eccentric_orbit();
    flux_horizon = sfmv_scalar_energy_flux_horizon;
}

void ScalarFluxMode::sfmm_get_energy_flux_infinity_eccentric_orbit(float_type& flux_infinity)
{
    if (sfmv_options.backend != ScalarRadialBackendKind::Original &&
        sfmv_angular_frequency_omega == 0) {
        sfmm_zero_flux_and_amplitudes();
        flux_infinity = 0;
        return;
    }

    sfmm_compute_energy_flux_infinity_eccentric_orbit();
    flux_infinity = sfmv_scalar_energy_flux_infinity;
}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_infinity_eccentric_orbit(float_type& flux_horizon, float_type& flux_infinity)
{
    if (sfmv_options.backend != ScalarRadialBackendKind::Original &&
        sfmv_angular_frequency_omega == 0) {
        sfmm_zero_flux_and_amplitudes();
        flux_horizon = 0; flux_infinity = 0;
        return;
    }

    sfmm_compute_energy_flux_horizon_infinity_eccentric_orbit();
    flux_horizon = sfmv_scalar_energy_flux_horizon;
    flux_infinity = sfmv_scalar_energy_flux_infinity;
}

void ScalarFluxMode::sfmm_compute_energy_flux_horizon_circular_orbit()
{
    float_type r_radius = sfmv_semilatus_rectum_p;

    complex_type R_up;
    sfmv_homsol->evaluate_up(r_radius, R_up, true);

    float_type phi = 0;
    float_type x = 0.;
    complex_type S{0};
    sfmv_swh->swshm_get_S(x, phi, S);

    complex_type alpha = - 4. * constants::pi * r_radius * r_radius * S/sfmv_Gamma/sfmv_W;

    complex_type delta_phi_minus = alpha*R_up;
    sfmm_record_amplitude(true, delta_phi_minus);

    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;

    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );
}

void ScalarFluxMode::sfmm_compute_energy_flux_infinity_circular_orbit()
{
    float_type r_radius = sfmv_semilatus_rectum_p;
    complex_type R_in;
    sfmv_homsol->evaluate_in(r_radius, R_in, true);

    float_type phi = 0;
    float_type x = 0.;
    complex_type S{0};
    sfmv_swh->swshm_get_S(x, phi, S);

    complex_type alpha = - 4. * constants::pi * r_radius * r_radius * S/sfmv_Gamma/sfmv_W;

    complex_type delta_phi_plus = alpha*R_in;
    sfmm_record_amplitude(false, delta_phi_plus);

    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
}

void ScalarFluxMode::sfmm_compute_energy_flux_horizon_infinity_circular_orbit()
{
    float_type r_radius = sfmv_semilatus_rectum_p;
    complex_type R_in;
    sfmv_homsol->evaluate_in(r_radius, R_in, true);

    complex_type R_up;
    sfmv_homsol->evaluate_up(r_radius, R_up, true);

    float_type phi = 0;
    float_type x = 0.;
    complex_type S{0};
    sfmv_swh->swshm_get_S(x, phi, S);

    complex_type alpha = - 4. * constants::pi * r_radius * r_radius * S/sfmv_Gamma/sfmv_W;

    complex_type delta_phi_plus = alpha*R_in;
    sfmm_record_amplitude(false, delta_phi_plus);
    complex_type delta_phi_minus = alpha*R_up;
    sfmm_record_amplitude(true, delta_phi_minus);

    float_type r_plus = sfmv_mass_black_hole_M + mp::pow(mp::pow(sfmv_mass_black_hole_M, 2) - mp::pow(sfmv_black_hole_spin_a, 2), 0.5);
    float_type omega_plus = sfmv_black_hole_spin_a/(2*sfmv_mass_black_hole_M*r_plus);
    float_type k_mkn = sfmv_angular_frequency_omega - sfmv_m * omega_plus;

    sfmv_scalar_energy_flux_infinity = (1./16./constants::pi) * mp::pow(sfmv_angular_frequency_omega , 2) * mp::pow(mp::abs(delta_phi_plus), 2);
    sfmv_scalar_energy_flux_horizon =  (1./16./constants::pi) * sfmv_angular_frequency_omega * k_mkn *      mp::pow(mp::abs(delta_phi_minus),2)*(r_plus*r_plus + sfmv_black_hole_spin_a*sfmv_black_hole_spin_a );
}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_circular_orbit(float_type& flux_horizon)
{
    if (sfmv_options.backend != ScalarRadialBackendKind::Original &&
        sfmv_angular_frequency_omega == 0) {
        sfmm_zero_flux_and_amplitudes();
        flux_horizon = 0;
        return;
    }

    sfmm_compute_energy_flux_horizon_circular_orbit();
    flux_horizon = sfmv_scalar_energy_flux_horizon;
}

void ScalarFluxMode::sfmm_get_energy_flux_infinity_circular_orbit(float_type& flux_infinity)
{
    if (sfmv_options.backend != ScalarRadialBackendKind::Original &&
        sfmv_angular_frequency_omega == 0) {
        sfmm_zero_flux_and_amplitudes();
        flux_infinity = 0;
        return;
    }

    sfmm_compute_energy_flux_infinity_circular_orbit();
    flux_infinity = sfmv_scalar_energy_flux_infinity;
}

void ScalarFluxMode::sfmm_get_energy_flux_horizon_infinity_circular_orbit(float_type& flux_horizon, float_type& flux_infinity)
{
    if (sfmv_options.backend != ScalarRadialBackendKind::Original &&
        sfmv_angular_frequency_omega == 0) {
        sfmm_zero_flux_and_amplitudes();
        flux_horizon = 0; flux_infinity = 0;
        return;
    }

    sfmm_compute_energy_flux_horizon_infinity_circular_orbit();
    flux_horizon = sfmv_scalar_energy_flux_horizon;
    flux_infinity = sfmv_scalar_energy_flux_infinity;
}

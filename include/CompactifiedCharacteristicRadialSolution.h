#pragma once

#include "config.h"
#include "HorizonStaticReference.h"
#include <array>
#include <cstddef>
#include <vector>
#include <boost/numeric/odeint.hpp>

// M=1, s=0; R_H = exp(i H (omega r_star - m chi)) y_H / r.
class CompactifiedCharacteristicRadialSolution {
public:
    struct Options {
        float_type endpoint_tolerance = float_type("1e-34");
        float_type integration_absolute_tolerance = float_type("1e-26");
        float_type integration_relative_tolerance = float_type("1e-26");
        float_type interpolation_tolerance = float_type("1e-26");
        float_type outgoing_tolerance = 0;
        unsigned horizon_order = 64;
        unsigned infinity_max_order = 160;
        unsigned max_endpoint_retries = 10;
        std::size_t max_steps = 200000;
        bool direct_static_startup = false;
    };

    CompactifiedCharacteristicRadialSolution(float_type a, float_type lambda,
        float_type omega, int l, int m, float_type r_min, float_type r_max);
    CompactifiedCharacteristicRadialSolution(float_type a, float_type lambda,
        float_type omega, int l, int m, float_type r_min, float_type r_max,
        const Options& options);

    void evaluate_in(const float_type& r, complex_type& R, complex_type& Rr,
        complex_type& Rrr) const;
    void evaluate_up(const float_type& r, complex_type& R, complex_type& Rr,
        complex_type& Rrr) const;
    void characteristic(int H, const float_type& r, complex_type& y,
        complex_type& y_sigma) const;
    void reconstruct(int H, const float_type& r, const complex_type& y,
        const complex_type& y_sigma, complex_type& R, complex_type& Rr,
        complex_type& Rrr) const;
    static void operator_coefficients(const float_type& a, const float_type& lambda,
        const float_type& omega, int m, int H, const float_type& sigma,
        complex_type& A, complex_type& B, complex_type& C);
    static float_type tortoise(const float_type& a, const float_type& r);
    static float_type azimuthal_phase(const float_type& a, const float_type& r);
    complex_type horizon_amplitude() const { return T_H_; }
    complex_type infinity_amplitude() const { return T_inf_; }
    complex_type wronskian(const float_type& r) const;
    complex_type wronskian_characteristic(const float_type& r) const;
    void rescale_basis(const complex_type& in_scale, const complex_type& up_scale);

private:
    using State = std::array<complex_type,2>;
    using Stepper = boost::numeric::odeint::bulirsch_stoer_dense_out<State,
        float_type, State, float_type>;
    struct DenseStep {
        float_type left, right;
        Stepper snapshot;
        DenseStep(float_type l, float_type r, const Stepper& s)
            : left(std::move(l)), right(std::move(r)), snapshot(s) {}
    };
    struct TaylorStep {
        float_type left, right;
        std::vector<complex_type> coefficients;
    };
    float_type a_,lambda_,omega_,rp_,rm_,b_,j_,w_,r_min_,r_max_,sigma_min_,sigma_max_;
    int l_,m_;
    HorizonStaticReference horizon_reference_;
    Options options_;
    complex_type T_H_, T_inf_, in_scale_=1, up_scale_=1;
    std::vector<DenseStep> in_dense_, up_dense_;
    std::vector<TaylorStep> in_taylor_, up_taylor_;
    State in_seed_{}, up_seed_{};
    float_type in_seed_sigma_=0, up_seed_sigma_=0;
    void validate();
    void initialize_endpoints();
    void integrate_static(int H, float_type start, const State& initial);
    void integrate_polynomial_taylor(int H, float_type start, const State& initial);
    void static_horizon_polynomials(std::vector<complex_type>& A,
        std::vector<complex_type>& B, std::vector<complex_type>& C) const;
    void coefficients(int H, const float_type& sigma, complex_type& A,
        complex_type& B, complex_type& C) const;
    void ingoing_coefficients(const float_type& sigma, complex_type& A,
        complex_type& B, complex_type& C) const;
    float_type checked_source_radius(const float_type& r) const;
    State state_at(int H, const float_type& sigma, bool characteristic_state=true) const;
    State horizon_series(float_type& start) const;
    State static_horizon_series(float_type& start) const;
    State infinity_series(float_type& start) const;
    State static_series(float_type& start) const;
};

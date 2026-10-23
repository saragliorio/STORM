#include "ScalarRadialBackend.h"
#include <stdexcept>

ScalarRadialBackend::ScalarRadialBackend(int s, float_type M, float_type a,
    float_type particle_mass, float_type lambda, float_type omega,
    int l, int azimuthal_m, float_type r_min, float_type r_max,
    ScalarRadialBackendKind kind)
    : ScalarRadialBackend(s,M,a,particle_mass,lambda,omega,l,azimuthal_m,
        r_min,r_max,kind,CompactifiedCharacteristicRadialSolution::Options{}) {}

ScalarRadialBackend::ScalarRadialBackend(int s, float_type M, float_type a,
    float_type particle_mass, float_type lambda, float_type omega,
    int l, int azimuthal_m, float_type r_min, float_type r_max,
    ScalarRadialBackendKind kind,
    const CompactifiedCharacteristicRadialSolution::Options& options,
    bool bypass_original_small_frequency_guard)
    : kind_(kind), M_(M), a_(a), r_min_(r_min), r_max_(r_max) {
    if (kind_ != ScalarRadialBackendKind::Original && (s != 0 || M != 1))
        throw std::invalid_argument("Compactified scalar solver supports only s=0 and M=1");
    if (kind_ == ScalarRadialBackendKind::Endpoint) {
        compactified_ = std::make_unique<CompactifiedCharacteristicRadialSolution>(
            a,lambda,omega,l,azimuthal_m,r_min,r_max,options);
        return;
    }
    original_ = std::make_unique<RadialHomogeneousSolution>(
        s,M,a,particle_mass,lambda,omega,l,azimuthal_m);

    complex_type nu, B_inc;
    original_->rhsm_get_nu(nu);
    original_->rhsm_get_B_inc(B_inc);
    if (!bypass_original_small_frequency_guard && mp::abs(omega) <= constants::precision_threshold)
        original_W_ = -(2.*l+1);
    else
        original_W_ = 2.*constants::I*omega*B_inc;
}

void ScalarRadialBackend::evaluate_in(const float_type& r, complex_type& R,
    bool circular) {
    if (kind_ == ScalarRadialBackendKind::Original) {
        if (circular) original_->rhsm_get_R_in_mst(r,R);
        else original_->rhsm_get_R_in_numerical_integration(r_min_,r_max_,r,R);
    } else {
        complex_type Rr,Rrr;
        compactified_->evaluate_in(r,R,Rr,Rrr);
    }
}
void ScalarRadialBackend::evaluate_up(const float_type& r, complex_type& R,
    bool circular) {
    if (kind_ == ScalarRadialBackendKind::Original) {
        if (circular) original_->rhsm_get_R_up_mst(r,R);
        else original_->rhsm_get_R_up_numerical_integration(r_min_,r_max_,r,R);
    } else {
        complex_type Rr,Rrr;
        compactified_->evaluate_up(r,R,Rr,Rrr);
    }
}
void ScalarRadialBackend::evaluate_in(const float_type& r, complex_type& R,
    complex_type& Rr, complex_type& Rrr, bool circular) {
    if (kind_ == ScalarRadialBackendKind::Original) {
        if (circular) original_->rhsm_get_R_in_and_1der_and_2der_mst(r,R,Rr,Rrr);
        else original_->rhsm_get_R_in_numerical_integration_and_1der_and_2der(
            r_min_,r_max_,r,R,Rr,Rrr);
    } else compactified_->evaluate_in(r,R,Rr,Rrr);
}
void ScalarRadialBackend::evaluate_up(const float_type& r, complex_type& R,
    complex_type& Rr, complex_type& Rrr, bool circular) {
    if (kind_ == ScalarRadialBackendKind::Original) {
        if (circular) original_->rhsm_get_R_up_and_1der_and_2der_mst(r,R,Rr,Rrr);
        else original_->rhsm_get_R_up_numerical_integration_and_1der_and_2der(
            r_min_,r_max_,r,R,Rr,Rrr);
    } else compactified_->evaluate_up(r,R,Rr,Rrr);
}
complex_type ScalarRadialBackend::wronskian(const float_type& r) const {
    return kind_ == ScalarRadialBackendKind::Original
        ? original_W_ : compactified_->wronskian(r);
}
complex_type ScalarRadialBackend::local_wronskian(const float_type& r, bool circular) {
    complex_type Ri,Rir,Rirr,Ru,Rur,Rurr;
    evaluate_in(r,Ri,Rir,Rirr,circular);
    evaluate_up(r,Ru,Rur,Rurr,circular);
    return (r*r-2*M_*r+a_*a_)*(Ri*Rur-Ru*Rir);
}
complex_type ScalarRadialBackend::horizon_amplitude() const {
    return kind_ == ScalarRadialBackendKind::Original
        ? complex_type(1) : compactified_->horizon_amplitude();
}
complex_type ScalarRadialBackend::infinity_amplitude() const {
    return kind_ == ScalarRadialBackendKind::Original
        ? complex_type(1) : compactified_->infinity_amplitude();
}
void ScalarRadialBackend::rescale_basis(const complex_type& in_scale,
    const complex_type& up_scale) {
    if (!compactified_)
        throw std::logic_error("Basis rescaling requires the compactified backend");
    compactified_->rescale_basis(in_scale,up_scale);
}

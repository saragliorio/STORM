#pragma once

#include "CompactifiedCharacteristicRadialSolution.h"
#include "RadialHomogeneousSolution.h"
#include <memory>

enum class ScalarRadialBackendKind {
    Original,
    Endpoint
};

class ScalarRadialBackend {
public:
    ScalarRadialBackend(int s, float_type M, float_type a,
        float_type particle_mass, float_type lambda, float_type omega,
        int l, int azimuthal_m, float_type r_min, float_type r_max,
        ScalarRadialBackendKind kind = ScalarRadialBackendKind::Original);
    ScalarRadialBackend(int s, float_type M, float_type a,
        float_type particle_mass, float_type lambda, float_type omega,
        int l, int azimuthal_m, float_type r_min, float_type r_max,
        ScalarRadialBackendKind kind,
        const CompactifiedCharacteristicRadialSolution::Options& options,
        bool bypass_original_small_frequency_guard = false);

    void evaluate_in(const float_type& r, complex_type& R,
        bool circular = false);
    void evaluate_up(const float_type& r, complex_type& R,
        bool circular = false);
    void evaluate_in(const float_type& r, complex_type& R, complex_type& Rr,
        complex_type& Rrr, bool circular = false);
    void evaluate_up(const float_type& r, complex_type& R, complex_type& Rr,
        complex_type& Rrr, bool circular = false);
    complex_type wronskian(const float_type& r) const;
    complex_type local_wronskian(const float_type& r, bool circular = false);
    complex_type horizon_amplitude() const;
    complex_type infinity_amplitude() const;
    ScalarRadialBackendKind kind() const { return kind_; }
    const CompactifiedCharacteristicRadialSolution* compactified() const {
        return compactified_.get();
    }
    void rescale_basis(const complex_type& in_scale, const complex_type& up_scale);

private:
    ScalarRadialBackendKind kind_;
    float_type M_, a_, r_min_, r_max_;
    complex_type original_W_;
    std::unique_ptr<RadialHomogeneousSolution> original_;
    std::unique_ptr<CompactifiedCharacteristicRadialSolution> compactified_;
};

#include "CompactifiedCharacteristicRadialSolution.h"
#include "ScalarFluxMode.h"
#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
using Solver = CompactifiedCharacteristicRadialSolution;
const complex_type I(0, 1);
const float_type algebra_tolerance = mp::pow(float_type(10), 20-std::numeric_limits<float_type>::digits10);
const float_type numerical_tolerance("1e-16");
int checks = 0;

void require(bool condition, const std::string& label) {
    ++checks;
    if (!condition) throw std::runtime_error(label);
}

void near(const complex_type& value, const complex_type& reference,
          const std::string& label, const float_type& tolerance = algebra_tolerance) {
    const float_type error = reference == 0 ? mp::abs(value) : mp::abs(value-reference)/mp::abs(reference);
    require(error <= tolerance, label+": error "+error.str(12));
}

template<class Exception, class Operation>
void rejects(Operation operation, const std::string& label) {
    bool caught = false;
    try { operation(); } catch (const Exception&) { caught = true; }
    require(caught, label);
}

void operator_and_reconstruction() {
    for (const char* spin : {"0", "0.9"})
      for (int m : {-2, 0, 2})
        for (const char* frequency : {"-0.07", "0", "0.07"})
          for (int H : {-1, 1})
            for (const char* coordinate : {"0.00001", "0.03", "0.8", "0.999"}) {
                const float_type a(spin), omega(frequency), lambda("6.013"), sigma(coordinate);
                const float_type root = mp::sqrt(1-a*a), rp = 1+root, rm = 1-root, r = rp/sigma;
                const float_type delta = (r-rp)*(r-rm), delta_r = 2*r-2;
                const float_type q = -rp/(r*r), qr = 2*rp/(r*r*r);
                const float_type K = omega*(r*r+a*a)-a*m;
                const complex_type L = I*H*K/delta-1/r;
                const complex_type Lr = I*H*(2*omega*r*delta-K*delta_r)/(delta*delta)+1/(r*r);
                complex_type A, B, C;
                Solver::operator_coefficients(a, lambda, omega, m, H, sigma, A, B, C);
                near(A, delta*q*q, "transformed physical operator A");
                near(B, delta*(qr+2*L*q)+delta_r*q, "transformed physical operator B");
                near(C, delta*(Lr+L*L)+delta_r*L+K*K/delta-lambda, "transformed physical operator C");
                const float_type rstar = r+rp/root*mp::log((r-rp)/2)-rm/root*mp::log((r-rm)/2);
                near(Solver::tortoise(a, r), rstar, "tortoise additive constant");
                near(Solver::azimuthal_phase(a, r), a/(2*root)*mp::log((r-rp)/(r-rm)), "azimuthal phase");
            }
    const float_type a("0.9"), omega("0.07"), lambda("6.013"), r("9.137");
    const float_type root = mp::sqrt(1-a*a), rp = 1+root, rm = 1-root;
    Solver solver(a, lambda, omega, 2, 2, 4, 20);
    for (int H : {-1, 1}) {
        const float_type delta = (r-rp)*(r-rm), K = omega*(r*r+a*a)-2*a;
        const float_type rstar = r+rp/root*mp::log((r-rp)/2)-rm/root*mp::log((r-rm)/2);
        const float_type chi = a/(2*root)*mp::log((r-rp)/(r-rm));
        const complex_type G = mp::exp(I*H*(omega*rstar-2*chi))/r;
        const complex_type y(float_type("0.31"), float_type("-0.27"));
        const complex_type ys(float_type("-0.51"), float_type("0.92"));
        complex_type R, Rr, Rrr;
        solver.reconstruct(H, r, y, ys, R, Rr, Rrr);
        near(R, G*y, "physical field phase");
        near(Rr, G*(-rp/(r*r)*ys+(I*H*K/delta-1/r)*y), "physical derivative phase");
        near(delta*Rrr+(2*r-2)*Rr, -(K*K/delta-lambda)*R, "physical radial equation");
    }
    near(solver.wronskian(r), solver.wronskian_characteristic(r), "Wronskian reconstruction");
    near(solver.wronskian(4), solver.wronskian(20), "conserved Wronskian", numerical_tolerance);
}

void horizon_polynomial() {
    for (int l : {0, 1, 2, 4}) {
        HorizonStaticReference reference(l);
        for (const char* coordinate : {"1", "0.8", "0.25", "0.05"}) {
            const float_type s(coordinate), x = 2/s-1, xs = -2/(s*s), xss = 4/(s*s*s);
            float_type P = 1, Px = 0, Pxx = 0;
            if (l == 1) { P = x; Px = 1; }
            if (l == 2) { P = (3*x*x-1)/2; Px = 3*x; Pxx = 3; }
            if (l == 4) { P = (35*x*x*x*x-30*x*x+3)/8; Px = (140*x*x*x-60*x)/8; Pxx = (420*x*x-60)/8; }
            complex_type S, Ss, Sss;
            reference.evaluate(s, S, Ss, Sss);
            near(S, P, "static horizon Legendre polynomial");
            near(Ss, Px*xs, "static horizon first derivative");
            near(Sss, Pxx*xs*xs+Px*xss, "static horizon second derivative");
        }
    }
}

void static_solutions() {
    for (int l : {0, 1, 2}) {
        Solver solver(0, l*(l+1), 0, l, 0, 4, 20);
        for (const char* radius : {"4", "6.137", "12.231", "20"}) {
            const float_type r(radius), x = r-1, Q0 = mp::log((x+1)/(x-1))/2, Q0r = -1/(x*x-1);
            float_type P = 1, Pr = 0, Q = Q0, Qr = Q0r, in_scale = 1, up_scale = 1;
            if (l == 1) { P = x; Pr = 1; Q = x*Q0-1; Qr = Q0+x*Q0r; in_scale = float_type(1)/2; up_scale = 6; }
            if (l == 2) { P = (3*x*x-1)/2; Pr = 3*x; Q = P*Q0-3*x/2; Qr = Pr*Q0+P*Q0r-float_type(3)/2; in_scale = float_type(1)/2; up_scale = 30; }
            complex_type Rin, Rrin, Rrrin, Rup, Rrup, Rrrup;
            solver.evaluate_in(r, Rin, Rrin, Rrrin);
            solver.evaluate_up(r, Rup, Rrup, Rrrup);
            near(Rin, in_scale*P, "static ingoing field", numerical_tolerance);
            near(Rup, up_scale*Q, "static outgoing field", numerical_tolerance);
            near(r*Rrin, r*in_scale*Pr, "static ingoing derivative", numerical_tolerance);
            near(Rrup, up_scale*Qr, "static outgoing derivative", numerical_tolerance);
            near(solver.wronskian(r), -in_scale*up_scale, "static exact Wronskian", numerical_tolerance);
        }
    }
}

void radius_guards() {
    for (bool singleton : {false, true}) {
        const float_type lo = 8, hi = singleton ? 8 : 12;
        Solver solver(0, 2, float_type("0.05"), 1, 1, lo, hi);
        for (int side : {-1, 1}) {
            const float_type boundary = side == -1 ? lo : hi;
            const float_type roundoff = boundary+side*2*std::numeric_limits<float_type>::epsilon()*boundary;
            require(roundoff != boundary, "boundary perturbation is representable");
            for (int H : {-1, 1}) {
                std::array<complex_type, 3> exact, perturbed;
                auto evaluate = [&](const float_type& radius, std::array<complex_type, 3>& output) {
                    if (H == -1) solver.evaluate_in(radius, output[0], output[1], output[2]);
                    else solver.evaluate_up(radius, output[0], output[1], output[2]);
                };
                evaluate(boundary, exact);
                evaluate(roundoff, perturbed);
                require(exact == perturbed, "roundoff radius canonicalizes before reconstruction");
                rejects<std::out_of_range>([&] { evaluate(boundary*(1+side*float_type("1e-20")), perturbed); }, "out-of-range radius rejected");
                rejects<std::invalid_argument>([&] { evaluate(std::numeric_limits<float_type>::quiet_NaN(), perturbed); }, "NaN radius rejected");
            }
            require(solver.wronskian(boundary) == solver.wronskian(roundoff), "Wronskian radius canonicalization");
        }
    }
    rejects<std::invalid_argument>([] { Solver bad(1, 2, float_type("0.05"), 1, 1, 8, 12); }, "extremal spin rejected");
    rejects<std::invalid_argument>([] { Solver bad(0, 2, float_type("0.05"), 1, 2, 8, 12); }, "invalid azimuthal index rejected");
    rejects<std::invalid_argument>([] { Solver bad(0, 2, float_type("0.05"), 1, 1, 2, 12); }, "horizon source radius rejected");
    rejects<std::invalid_argument>([] { Solver bad(0, float_type("2.001"), 0, 1, 0, 8, 12); }, "invalid static separation constant rejected");
}

void source_quadrature_failure() {
    ScalarFluxOptions options;
    options.backend = ScalarRadialBackendKind::Endpoint;
    options.quadrature_tolerance = float_type("1e-100");
    ScalarFluxMode mode(0, 1, 1, float_type("0.9"), 10, float_type("0.3"), 45, 1, 1, 0, 0, options);
    float_type horizon = 1, infinity = 2;
    bool caught = false;
    try {
        mode.sfmm_get_energy_flux_horizon_infinity_generic_orbit(horizon, infinity);
    } catch (const std::runtime_error& error) {
        caught = std::string(error.what()) == "Scalar flux quadrature did not converge";
    }
    require(caught, "unconverged compactified source quadrature rejected");
    require(horizon == 1 && infinity == 2, "failed quadrature leaves flux outputs unchanged");
    complex_type Ah, Ai;
    rejects<std::logic_error>([&] { mode.sfmm_get_complex_source_amplitudes(Ah, Ai); }, "failed quadrature leaves amplitudes unavailable");
}

void eccentric_flux() {
    ScalarFluxOptions options;
    options.backend = ScalarRadialBackendKind::Endpoint;
    options.quadrature_tolerance = float_type("1e-18");
    ScalarFluxMode mode(0, 1, 1, float_type("0.9"), 8, float_type("0.3"), 0, 2, 2, 0, 1, options);
    float_type horizon, infinity, repeated_horizon, repeated_infinity;
    complex_type Ah, Ai, repeated_Ah, repeated_Ai;
    rejects<std::logic_error>([&] { mode.sfmm_get_complex_source_amplitudes(Ah, Ai); }, "amplitudes unavailable before integration");
    mode.sfmm_get_energy_flux_horizon_infinity_eccentric_orbit(horizon, infinity);
    mode.sfmm_get_complex_source_amplitudes(Ah, Ai);
    require(horizon < 0 && infinity > 0, "superradiant signed fluxes");
    mode.sfmm_get_energy_flux_horizon_eccentric_orbit(repeated_horizon);
    mode.sfmm_get_energy_flux_infinity_eccentric_orbit(repeated_infinity);
    near(repeated_horizon, horizon, "separate horizon source quadrature", numerical_tolerance);
    near(repeated_infinity, infinity, "separate infinity source quadrature", numerical_tolerance);
    mode.sfmm_rescale_radial_basis(complex_type(float_type("-1.7"), float_type("0.8")),
                                  complex_type(float_type("0.4"), float_type("-1.3")));
    rejects<std::logic_error>([&] { mode.sfmm_get_complex_source_amplitudes(repeated_Ah, repeated_Ai); }, "rescaling invalidates cached amplitudes");
    for (int repeat = 0; repeat < 2; ++repeat) {
        mode.sfmm_get_energy_flux_horizon_infinity_eccentric_orbit(repeated_horizon, repeated_infinity);
        mode.sfmm_get_complex_source_amplitudes(repeated_Ah, repeated_Ai);
        near(repeated_Ah, Ah, "horizon amplitude rescaling invariance", numerical_tolerance);
        near(repeated_Ai, Ai, "infinity amplitude rescaling invariance", numerical_tolerance);
        near(repeated_horizon, horizon, "horizon flux rescaling invariance", numerical_tolerance);
        near(repeated_infinity, infinity, "infinity flux rescaling invariance", numerical_tolerance);
    }
}
}

int main() {
    const float_type quadrature_tolerance = constants::precision_trapezoidal_quadrature;
    try {
        operator_and_reconstruction();
        horizon_polynomial();
        static_solutions();
        radius_guards();
        source_quadrature_failure();
        constants::precision_trapezoidal_quadrature = float_type("1e-100");
        eccentric_flux();
        constants::precision_trapezoidal_quadrature = quadrature_tolerance;
        std::cout << "PASS compactified_test checks=" << checks << '\n';
    } catch (const std::exception& error) {
        constants::precision_trapezoidal_quadrature = quadrature_tolerance;
        std::cerr << "FAIL compactified_test: " << error.what() << '\n';
        return 1;
    }
}

#include "ScalarFluxMode.h"
#include <boost/math/special_functions/next.hpp>
#include <stdexcept>

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

template<class Exception, class F>
void require_throws(F f, const char* message, const char* expected = nullptr)
{
    try { f(); }
    catch (const Exception& error) {
        require(!expected || std::string(error.what()) == expected, message);
        return;
    }
    throw std::runtime_error(message);
}

int main()
{
    try {
        for (int l : {0, 1}) {
            SpinWeightedSpheroidalHarmonics angular(0, l, 0, float_type("0.1"));
            float_type lambda;
            angular.swshm_get_lambda(lambda);
            require(mp::isfinite(lambda) && mp::abs(lambda-l*(l+1)) < float_type("0.02"),
                    "Nonzero-spheroidicity scalar m=0 eigenvalue failed");
        }

        // Independent Gamma values computed by mpmath at 185 decimal digits.
        for (bool reflected : {false, true}) {
            const complex_type z(float_type(reflected ? "-2.5" : "0.5"), float_type("0.7"));
            const complex_type reference(
                float_type(reflected ? "-0.1598187163629329021126466705447361386373136396019200000853745085087394743725827390567081"
                                     : "0.5274198573183713735855009636998388629041407905745274712989309341887886363889873365608353"),
                float_type(reflected ? "-0.1575665490815152719398898329656713168534639575348074609796590228127845363747799132128053"
                                     : "-0.6404494748452189853189286858243899970242319270007559024376364178573976416435270165670250"));
            complex_type gamma, next, ratio;
            special_functions::sf_compute_gamma_function_z(z, gamma);
            special_functions::sf_compute_gamma_function_z(z+1, next);
            require(mp::abs((gamma-reference)/reference) < float_type("1e-75"), "Complex Gamma accuracy failed");
            require(mp::abs((next-z*gamma)/(z*gamma)) < float_type("1e-75"), "Gamma recurrence failed");
            complex_type aliased = z;
            special_functions::sf_compute_gamma_function_z(aliased, aliased);
            require(mp::abs((aliased-gamma)/gamma) < float_type("1e-75"), "Aliased Gamma evaluation failed");
            for (int n : {-1, 0, 1}) {
                special_functions::sf_compute_gamma_ratio(z, n, ratio);
                const complex_type expected = n < 0 ? 1/(z-1) : (n > 0 ? z : complex_type(1));
                require(mp::abs(ratio-expected) < float_type("1e-75"), "Integer Gamma ratio failed");
                aliased = z;
                special_functions::sf_compute_gamma_ratio(aliased, n, aliased);
                require(mp::abs(aliased-expected) < float_type("1e-75"), "Aliased Gamma ratio failed");
            }
        }

        constants::digits_precision_newton_algorithm = 180;
        constants::accuracy_continued_fraction_computation = float_type("1e-80");
        constants::accuracy_sum_computation = float_type("1e-75");
        constants::cons_counts = 12;
        RadialHomogeneousSolution tiny(0, 1, 0, 1, 0, float_type("5e-11"), 0, 0);
        complex_type analytic_value, derivative, value_only;
        tiny.rhsm_get_R_in_and_1der_mst(4, analytic_value, derivative);
        tiny.rhsm_get_R_in_mst(4, value_only);
        require(mp::abs((analytic_value-value_only)/value_only) < float_type("1e-50"), "MST value getters disagree");
        auto value = [&](const float_type& r) {
            complex_type result;
            tiny.rhsm_get_R_in_mst(r, result);
            return result;
        };
        for (const float_type& h : {float_type("1e-8"), float_type("1e-10")}) {
            const complex_type difference = (-value(4-3*h)+9*value(4-2*h)-45*value(4-h)
                                             +45*value(4+h)-9*value(4+2*h)+value(4+3*h))/(60*h);
            require(mp::abs((derivative-difference)/difference) < float_type("1e-35"),
                    "Tiny-frequency MST derivative disagrees with independent finite difference");
        }
        constants::imax_coefficients = 1;
        require_throws<std::out_of_range>([] {
            RadialHomogeneousSolution exhausted(0, 1, 0, 1, 0, float_type("5e-11"), 0, 0);
        }, "Exhausted MST coefficients did not throw");
        constants::imax_coefficients = 500;
        constants::digits_precision_newton_algorithm = 30;
        constants::accuracy_continued_fraction_computation = float_type("1e-40");
        constants::accuracy_sum_computation = float_type("1e-30");
        constants::cons_counts = 5;

        RadialHomogeneousSolution radial(0, 1, 0, 1, 2, float_type("0.1"), 1, 1);
        for (bool up : {false, true}) {
            for (bool derivatives : {false, true}) {
                auto sample = [&](const float_type& r) {
                    std::array<complex_type, 3> result{};
                    if (derivatives && up) radial.rhsm_get_R_up_numerical_integration_and_1der_and_2der(4, 10, r, result[0], result[1], result[2]);
                    else if (derivatives) radial.rhsm_get_R_in_numerical_integration_and_1der_and_2der(4, 10, r, result[0], result[1], result[2]);
                    else if (up) radial.rhsm_get_R_up_numerical_integration(4, 10, r, result[0]);
                    else radial.rhsm_get_R_in_numerical_integration(4, 10, r, result[0]);
                    return result;
                };
                const auto lower = sample(4), upper = sample(10);
                require(lower[0] != 0 && upper[0] != 0, "Dense radial endpoint returned zero");
                require(sample(bm::float_prior(float_type(4))) == lower, "Lower one-ULP endpoint failed");
                require(sample(bm::float_next(float_type(10))) == upper, "Upper one-ULP endpoint failed");
                require_throws<std::out_of_range>([&] { sample(float_type("3.99")); }, "Uncovered lower radial query did not throw");
                require_throws<std::out_of_range>([&] { sample(float_type("10.01")); }, "Uncovered upper radial query did not throw");
            }
        }

        RadialHomogeneousSolution uncovered(0, 1, 0, 1, 2, float_type("0.1"), 1, 1);
        require_throws<std::out_of_range>([&] {
            uncovered.rhsm_get_R_in_numerical_integration(4, 4, 4, value_only);
        }, "Missing horizon dense data did not throw");
        require_throws<std::out_of_range>([&] {
            uncovered.rhsm_get_R_up_numerical_integration(4, 4, 4, value_only);
        }, "Missing infinity dense data did not throw");

        constants::precision_trapezoidal_quadrature = float_type("1e-100");
        require_throws<std::runtime_error>([] {
            ScalarFluxMode flux(0, 1, 1, float_type("0.9"), 10, float_type("0.3"), 45, 1, 1, 0, 0);
            float_type horizon, infinity;
            flux.sfmm_get_energy_flux_horizon_infinity_generic_orbit(horizon, infinity);
        }, "Unconverged scalar quadrature did not throw", "Scalar flux quadrature did not converge");
        std::cout << "MST regressions passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

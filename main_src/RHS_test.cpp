#include "RadialHomogeneousSolution.h"
#include "GeodesicOrbitalMotion.h"
#include "SpinWeightedSpheroidalHarmonics.h"

// Example demonstrating how to compute homogeneous solutions
// of the radial Teukolsky equation using the RadialHomogeneousSolution class.
// The code computes:
// 1) Spin-weighted spheroidal eigenvalues
// 2) Renormalized angular momentum ν
// 3) Asymptotic amplitudes
// 4) Radial homogeneous solutions using the MST method
// 5) Radial homogeneous solutions using numerical integration

int main() 
{

    // Set numerical precision for output
    std::cout << std::setprecision(20);
    
    // ==================================
    // Example with arbitrary r and omega
    // ==================================

    int s = 0;                                   // Spin weight s
    std::cout << "s = " << s << std::endl;

    float_type a{"0.9"};                          // Kerr spin parameter
    float_type M{"1"};                            // Black hole mass
    float_type mu{"1"};                           // Test particle mass
    std::cout << "a = " << a << std::endl;

    int l = 1;                                    // Spheroidal harmonic degree
    int m = 1;                                    // Azimuthal number
    std::cout << "l = " << l << std::endl;
    std::cout << "m = " << m << std::endl;

    float_type r{"8"};                            // Radius of the test particle
    float_type omega{"4.1"};                      // Angular frequency
    std::cout << "r = " << r << std::endl;
    std::cout << "omega = " << omega << std::endl;


    // Construct SpinWeightedSpheroidalHarmonics object
    SpinWeightedSpheroidalHarmonics* swsh_obj = new SpinWeightedSpheroidalHarmonics(s, l, m, a*omega);

    // Retrieve angular eigenvalue λ
    float_type lambda;
    swsh_obj->swshm_get_lambda(lambda);
    std::cout << "lambda = " << lambda << std::endl;

    // Construct RadialHomogeneousSolution object
    RadialHomogeneousSolution* rhs_obj = new RadialHomogeneousSolution(s, M, a, mu, lambda, omega, l, m);

    // Retrieve renormalized angular momentum ν
    complex_type nu{0};
    rhs_obj->rhsm_get_nu(nu);
    std::cout << "nu = " << nu << std::endl;
    
    // Retrieve asymptotic amplitudes
    complex_type B_trans;                         // Transmission amplitude at horizon
    complex_type C_trans;                         // Up-mode transmission amplitude
    complex_type B_inc;                           // Incoming amplitude at infinity
    complex_type B_ref;                           // Reflected amplitude at infinity

    rhs_obj->rhsm_get_B_trans(B_trans);
    rhs_obj->rhsm_get_C_trans(C_trans);
    rhs_obj->rhsm_get_B_inc(B_inc);
    rhs_obj->rhsm_get_B_ref(B_ref);

    std::cout << "B_trans = " << B_trans << std::endl;
    std::cout << "C_trans = " << C_trans << std::endl;
    std::cout << "B_inc = " << B_inc << std::endl;
    std::cout << "B_ref = " << B_ref << std::endl;

    // Compute radial solutions using MST formalism
    complex_type R_in;                            // Ingoing homogeneous solution
    complex_type R_up;                            // Outgoing homogeneous solution

    complex_type R_in_1der;                       // First derivative of R_in
    complex_type R_in_2der;                       // Second derivative of R_in
    complex_type R_up_1der;                       // First derivative of R_up
    complex_type R_up_2der;                       // Second derivative of R_up

    // Evaluate MST radial solutions at radius r
    rhs_obj->rhsm_get_R_in_and_1der_and_2der_mst(r, R_in, R_in_1der, R_in_2der);
    rhs_obj->rhsm_get_R_up_and_1der_and_2der_mst(r, R_up, R_up_1der, R_up_2der);

    std::cout << "R_in (MST method) = " << R_in << std::endl;
    std::cout << "R_up (MST method) = " << R_up << std::endl;
    std::cout << "R_in_1der (MST method) = " << R_in_1der << std::endl;
    std::cout << "R_up_1der (MST method) = " << R_up_1der << std::endl;
    std::cout << "R_in_2der (MST method) = " << R_in_2der << std::endl;
    std::cout << "R_up_2der (MST method) = " << R_up_2der << std::endl;
    
    // =============================
    // Example using orbital motion
    // =============================

    int n = 0;                                    // Radial harmonic index
    int k = 0;                                    // Polar harmonic index
    float_type e{"0.3"};                           // Orbital eccentricity
    float_type theta{"1.6"};                       // Inclination parameter

    // Construct geodesic orbital motion object
    GeodesicOrbitalMotion* gom_obj = new GeodesicOrbitalMotion(M, m, a, r, e, theta, l, m, k, n);

    // Retrieve orbital frequency
    float_type angular_frequency; 
    gom_obj->gomm_get_angular_frequency_omega(angular_frequency);
    
    // Compute spheroidal eigenvalue corresponding to orbital frequency
    SpinWeightedSpheroidalHarmonics* swsh_obj_1 = new SpinWeightedSpheroidalHarmonics(s, l, m, a*angular_frequency);

    swsh_obj_1->swshm_get_lambda(lambda);
    std::cout << "lambda = " << lambda << std::endl;

    // Construct new radial solution object
    RadialHomogeneousSolution* rhs_obj_1 = new RadialHomogeneousSolution(s, M, a, mu, lambda, omega, l, m);

    // Radial turning points of eccentric orbit
    float_type r_min = r/(1 + e);                  // Periastron radius
    float_type r_max = r/(1 - e);                  // Apastron radius

    std::cout << "r_min = " << r_min << std::endl;
    std::cout << "r_max = " << r_max << std::endl;

    // Compute radial solutions via numerical integration
    rhs_obj_1->rhsm_get_R_in_numerical_integration_and_1der_and_2der(
        r_min, r_max, r, R_in, R_in_1der, R_in_2der);

    rhs_obj_1->rhsm_get_R_up_numerical_integration_and_1der_and_2der(
        r_min, r_max, r, R_up, R_up_1der, R_up_2der);

    std::cout << "R_in (numerical integration method) = " << R_in << std::endl;
    std::cout << "R_up (numerical integration method) = " << R_up << std::endl;

    std::cout << "R_in_1der (numerical integration method) = " << R_in_1der << std::endl;
    std::cout << "R_in_2der (numerical integration method) = " << R_in_2der << std::endl;

    std::cout << "R_up_1der (numerical integration method) = " << R_up_1der << std::endl;
    std::cout << "R_up_2der (numerical integration method) = " << R_up_2der << std::endl;

    // Evaluate MST radial solutions at radius r
    rhs_obj_1->rhsm_get_R_in_and_1der_and_2der_mst(r, R_in, R_in_1der, R_in_2der);
    rhs_obj_1->rhsm_get_R_up_and_1der_and_2der_mst(r, R_up, R_up_1der, R_up_2der);

    std::cout << "R_in (MST method) = " << R_in << std::endl;
    std::cout << "R_up (MST method) = " << R_up << std::endl;
    std::cout << "R_in_1der (MST method) = " << R_in_1der << std::endl;
    std::cout << "R_up_1der (MST method) = " << R_up_1der << std::endl;
    std::cout << "R_in_2der (MST method) = " << R_in_2der << std::endl;
    std::cout << "R_up_2der (MST method) = " << R_up_2der << std::endl;


    delete swsh_obj;
    delete rhs_obj;
    delete gom_obj;
    delete swsh_obj_1;
    delete rhs_obj_1;
}
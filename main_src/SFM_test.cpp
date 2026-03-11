#include "ScalarFluxMode.h"

/*
 Example demonstrating the usage of the ScalarFluxMode class.

 This program initializes a Kerr black hole system with a test particle orbit
 and a specific scalar mode (s, l, m, k, n). It then computes and prints:
  - scalar energy flux at the black hole horizon
  - scalar energy flux at infinity

 The program uses the ScalarFluxMode class, which internally handles:
  - orbital motion via GeodesicOrbitalMotion
  - angular functions via SpinWeightedSpheroidalHarmonics
  - radial homogeneous solutions via RadialHomogeneousSolution
  - integration over the source terms to obtain energy fluxes
*/

int main()
{
    // Set numerical precision for output
    std::cout << std::setprecision(20) << std::endl;

    // --- Orbital and black hole parameters ---
    int s = 0;                       // spin weight of the scalar field
    float_type M("1.0");             // mass of the primary black hole
    float_type mu("1.0");            // mass of the orbiting particle
    float_type a("0.9");             // black hole spin
    float_type p("8.0");             // semilatus rectum of orbit
    float_type e("0.3");             // orbital eccentricity
    float_type theta_inc("20.0");    // inclination angle of orbit

    std::cout << "a = " << a << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "e = " << e << std::endl;
    std::cout << "theta_inc = " << theta_inc << std::endl;

    // --- Harmonic indices for the scalar mode ---
    int l = 1;    // spherical harmonic index
    int m = 1;    // azimuthal harmonic index
    int k = 0;    // polar harmonic index
    int n = 0;    // radial harmonic index

    std::cout << "(l, m, n, k) = (" << l  << ", " << m << ", " << n <<", " << k  << ") "<< std::endl;

    // --- Initialize ScalarFluxMode object for this mode and orbit ---
    ScalarFluxMode sfm_obj = ScalarFluxMode(s, M, m, a, p, e, theta_inc, l, m, k, n);
    
    // --- Compute scalar energy flux for a generic orbit ---
    float_type flux_horizon;     // flux at black hole horizon
    float_type flux_infinity;    // flux at infinity
    sfm_obj.sfmm_get_energy_flux_horizon_infinity_generic_orbit(flux_horizon, flux_infinity);  
    
    // --- Output results ---
    std::cout << "Energy flux at the horizon: " << flux_horizon << std::endl;
    std::cout << "Energy flux at infinity: " << flux_infinity << std::endl;
}
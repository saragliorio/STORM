#include "GeodesicOrbitalMotion.h"

/*
 Example program demonstrating the usage of the GeodesicOrbitalMotion class.
 
 This program initializes a Kerr black hole system and a test particle orbit.
 It then computes and prints:
  - constants of motion
  - orbital frequencies (Mino time and coordinate time)
  - combined angular frequency
  - orbital coordinate functions
  - contributions to coordinate time and azimuthal motion
 */

int main() 
{
    // Set numerical precision for output
    std::cout << std::setprecision(20) << std::endl;
    
    // --- Physical parameters of the system ---
    float_type M{"1"};           // mass of the primary black hole
    float_type mu{"1"};          // mass of the orbiting particle
    float_type a{"0.9"};         // black hole spin parameter

    // --- Orbital parameters ---
    float_type p{"8"};           // semilatus rectum
    float_type e {"0.3"};        // orbital eccentricity
    float_type theta_inc{"5.3"}; // orbital inclination angle

    std::cout << "a = " << a << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "e = " << e << std::endl;
    std::cout << "theta_inc = " << theta_inc << std::endl;
    
    // --- Harmonic indices associated with the orbit ---
    int l = 1;   // spherical harmonic index
    int m = 1;   // azimuthal index
    int k = 1;   // polar harmonic index
    int n = 1;   // radial harmonic index
    
    std::cout << "(l, m, n, k) = (" << l  << ", " << m << ", " << n <<", " << k  << ") "<< std::endl;

    // create an instance of the geodesic orbital motion class
    GeodesicOrbitalMotion gom_obj = GeodesicOrbitalMotion (M, mu, a, p, e, theta_inc, l, m, k, n);

    // --- Get constants of motion ---
    float_type energy;
    float_type angular_momentum;
    float_type carter_constant;
    gom_obj.gomm_get_constants_of_motion(energy, angular_momentum, carter_constant);
    std::cout << "Energy = " << energy <<std::endl;
    std::cout << "Angular momentum = " << angular_momentum <<std::endl;
    std::cout << "Carter constant = " << carter_constant <<std::endl;
    
    // --- Get fundamental frequencies in Mino time ---
    float_type Upsilon_r;
    float_type Upsilon_theta;
    float_type Upsilon_phi;
    float_type Gamma;
    gom_obj.gomm_get_frequencies_mino_time(Upsilon_r, Upsilon_theta, Upsilon_phi, Gamma);
    std::cout << "Upsilon_r = " << Upsilon_r << std::endl;
    std::cout << "Upsilon_theta = " << Upsilon_theta << std::endl;
    std::cout << "Upsilon_phi = " << Upsilon_phi <<std::endl;
    std::cout << "Gamma = " << Gamma <<std::endl;

    // --- Get orbital frequencies in coordinate time ---
    float_type Omega_r;
    float_type Omega_theta;
    float_type Omega_phi;
    gom_obj.gomm_get_frequencies_standard_time(Omega_r, Omega_theta, Omega_phi);
    std::cout << "Omega_r = " << Omega_r << std::endl;
    std::cout << "Omega_theta = " << Omega_theta << std::endl;
    std::cout << "Omega_phi = " << Omega_phi <<std::endl;
    
    // --- Get the angular frequency ---
    float_type omega;
    gom_obj.gomm_get_angular_frequency_omega(omega);
    std::cout << "angular_frequency_omega = " << omega <<std::endl;

    
    // --- Orbital phase variables ---
    float_type qr{"2.4"};     // radial phase
    float_type qz{"1.8"};     // polar phase
    float_type qphi{"3.2"};   // azimuthal phase

    // --- Evaluate orbital coordinates ---
    float_type r, z;
    gom_obj.gomm_get_r_function(qr, r);
    gom_obj.gomm_get_z_function(qz, z);
    std::cout << "r = " << r <<std::endl;
    std::cout << "z = " << z <<std::endl;

    // --- Evaluate coordinate time contributions ---
    float_type tr, tz, t;
    gom_obj.gomm_get_tr_function(qr, tr);
    gom_obj.gomm_get_tz_function(qz, tz);
    gom_obj.gomm_get_t_function(qr, qz, qphi, t);
    std::cout << "tr = " << tr <<std::endl;
    std::cout << "tz = " << tz <<std::endl;
    std::cout << "t = " << t <<std::endl;

    // --- Evaluate azimuthal motion contributions ---
    float_type phi_r, phi_z, phi;        
    gom_obj.gomm_get_phir_function(qr, phi_r);
    gom_obj.gomm_get_phiz_function(qz, phi_z);
    gom_obj.gomm_get_phi_function(qr, qz, qphi, phi);
    std::cout << "phi_r = " << phi_r <<std::endl;
    std::cout << "phi_z = " << phi_z <<std::endl;
    std::cout << "phi = " << phi <<std::endl;

    return 0;
}
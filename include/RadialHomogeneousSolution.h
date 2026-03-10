#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/math/tools/fraction.hpp>
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/numeric/odeint.hpp>
#include "config.h"
#include "special_functions.h"

#ifndef HOMOGENEOUSESOLUTION_H
#define HOMOGENEOUSSOLUTION_H

//REFERENCES:
//  -  S. Mano, H. Suzuki, and E. Takasugi, Prog. Theor.Phys. 95, 1079 (1996) (for MST method)
//  -  M. Sasaki and H. Tagoshi, Living Rev. Rel. 6, 6 (2003) (for MST method)
//  -  Z. Nasipak, Classical and Quantum Gravity 42, 165001 (2025) (for monodromy method)
//  -  A. Zenginoğlu's 2011 Physical Review D 83 (127502) (for numerical integration method)


// Class that computes homogeneous solutions of the radial Teukolsky equation
// using the MST (Mano–Suzuki–Takasugi) formalism and numerical integration.
class RadialHomogeneousSolution 
{
    // rhsv stands for "Radial Homogeneous Solution Variable"

    // ===== Physical parameters of the Teukolsky equation =====
    int rhsv_spin_weight_s;                         // spin weight s
    float_type rhsv_mass_black_hole_M;              // mass of the primary black hole
    float_type rhsv_mass_test_particle_m;           // mass of the orbiting compact object
    float_type rhsv_black_hole_spin_a;              // spin parameter of the black hole
    float_type rhsv_angular_eingenvalue_lambda;     // Eigenvalue of the angular Teukolsky equation
    float_type rhsv_black_hole_spin_per_unit_mass_q;// Dimensionless spin q = a/M
    float_type rhsv_angular_frequency_omega;        // Orbital frequency of the test particle
    int rhsv_l;                                     // spherical harmonic index l
    int rhsv_m;                                     // azimuthal harmonic index m
    bool rhsv_negative_freq;                        // Flag indicating if the frequency is negative
    float_type rhsv_epsilon;                        // MST parameter 
    float_type rhsv_kappa;                          // MST parameter
    float_type rhsv_tau;                            // MST parameter

    // ===== Renormalized angular momentum =====
    complex_type rhsv_nu_root;                      // Root of the continued fraction equation defining ν
    complex_type rhsv_minus_nu_root_minus_1;        // Second root −ν−1 

    // ===== Minimal solutions of the MST recurrence relation =====
    std::vector<complex_type> rhsv_minimal_coefficient_n_positive_nu;  // Minimal solution coefficients for n>0 (ν branch)
    std::vector<complex_type> rhsv_minimal_coefficient_n_negative_nu;  // Minimal solution coefficients for n<0 (ν branch)
    std::vector<complex_type> rhsv_minimal_coefficient_n_positive_nu1; // Minimal solution coefficients for n>0 (−ν−1 branch)
    std::vector<complex_type> rhsv_minimal_coefficient_n_negative_nu1; // Minimal solution coefficients for n<0 (−ν−1 branch)

    // ===== Asymptotic amplitudes of the radial solutions =====
    complex_type rhsv_B_trans;                      // Transmission coefficient at the horizon
    complex_type rhsv_B_inc;                        // Incoming wave amplitude at infinity
    complex_type rhsv_B_ref;                        // Reflected wave amplitude at infinity
    complex_type rhsv_C_trans;                      // Transmission amplitude at infinity

    // ===== Flags controlling numerical integrations =====
    bool rhsv_numerical_integration_in;             // Flag indicating if the "in" solution was integrated numerically
    bool rhsv_numerical_integration_up;             // Flag indicating if the "up" solution was integrated numerically

    // Structures implementing continued fractions and numerical integrations
    struct continued_fraction_Rn;                   // Continued fraction representation for R_n
    struct continued_fraction_Ln;                   // Continued fraction representation for L_n
    struct numerical_integration;                   // Struct used for numerical integration from the horizon
    struct numerical_integration_infinity;          // Struct used for numerical integration from infinity
    
    // ===== Data structure storing dense output from ODE integration =====
    struct DenseStepData 
    {
        float_type t_start;
        float_type dt;
        std::array<complex_type, 2> x;

        boost::numeric::odeint::bulirsch_stoer_dense_out<std::array<complex_type, 2>, float_type,
        std::array<complex_type, 2>, float_type> stepper_snapshot;
    };

    std::vector<DenseStepData> rhsv_dense_data_horizon;   // Stored integration steps from the horizon
    std::vector<DenseStepData> rhsv_dense_data_infinity;  // Stored integration steps from infinity

    // ===== Internal methods =====

    void rhsm_compute_alpha_n (const int n, complex_type& alpha_n, const complex_type& nu); 
    // Computes α_n coefficient of the MST three-term recurrence relation

    void rhsm_compute_beta_n (const int n, complex_type& beta_n, const complex_type& nu); 
    // Computes β_n coefficient of the MST recurrence relation

    void rhsm_compute_gamma_n (const int n, complex_type& gamma_n, const complex_type& nu); 
    // Computes γ_n coefficient of the MST recurrence relation

    void rhsm_compute_Rn (const int n, complex_type& Rn, const complex_type& nu); 
    // Computes R_n ratio used in the continued fraction for minimal solutions
                                                            
    void rhsm_compute_Ln (const int n, complex_type& Ln, const complex_type& nu); 
    // Computes L_n ratio used in the continued fraction for minimal solutions

    void rhsm_function_renormalized_angular_momentum (const complex_type& nu, complex_type& f_nu);
    // Function whose root determines the renormalized angular momentum ν

    void rhsm_function_renormalized_angular_momentum_derivative(const complex_type& nu, complex_type& result); 
    // Derivative of the ν root equation used in Newton iteration

    void rhsm_newton_algorithm_complex();
    // Newton method to find the complex root ν

    void rhsm_monodromy_method_for_nu();
    // Alternative way to compute ν using monodromy method

    void rhsm_compute_minimal_coefficient_n_positive(const int n, const complex_type& nu, std::vector<complex_type>& minimal_coefficient_n_positive);
    // Computes minimal MST coefficients for n>0

    void rhsm_compute_minimal_coefficient_n_negative(const int n, const complex_type& nu, std::vector<complex_type>& minimal_coefficient_n_negative);
    // Computes minimal MST coefficients for n<0

    void rhsm_sum_minimal_coefficients(const complex_type& nu, const std::vector<complex_type>& minimal_coefficient_n_positive, const std::vector<complex_type>& minimal_coefficient_n_negative , complex_type& result);
    // Computes the sum of minimal coefficients

    void rhsm_compute_A_minus(const complex_type& nu, const std::vector<complex_type>& minimal_coefficient_n_positive, const std::vector<complex_type>& minimal_coefficient_n_negative, complex_type& A_minus);
    // Computes the MST coefficient A⁻

    void rhsm_compute_K(const complex_type& nu, const std::vector<complex_type>& minimal_coefficient_n_positive, const std::vector<complex_type>& minimal_coefficient_n_negative, complex_type& K);
    // Computes the normalization constant K of the MST formalism

    void rhsm_compute_A_plus(const complex_type& nu, const std::vector<complex_type>& minimal_coefficient_n_positive, const std::vector<complex_type>& minimal_coefficient_n_negative, const complex_type& sum_minimal_coefficients, complex_type& A_plus);
    // Computes the MST coefficient A⁺

    void rhsm_compute_constants_quantities();
    // Computes all normalization constants and asymptotic amplitudes

    void rhsm_compute_R0(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, complex_type& R0, const float_type& x );
    // Computes the MST hypergeometric radial solution R₀

    void rhsm_compute_R0_and_1der(const complex_type& nu, const std::vector<complex_type>& minimal_coefficient_n_positive, const std::vector<complex_type>& minimal_coefficient_n_negative, complex_type& R0, complex_type& R0_der, const float_type& x);
    // Computes R₀ and its first radial derivative

    void rhsm_compute_R0_and_1der_and_2der(const complex_type& nu, const std::vector<complex_type>& minimal_coefficient_n_positive, const std::vector<complex_type>& minimal_coefficient_n_negative, complex_type& R0, complex_type& R0_der, complex_type& R0_2der, const float_type& x);
    // Computes R₀ and its first and second derivatives

    void rhsm_compute_Rup(const complex_type& nu, const std::vector<complex_type>& rhsv_minimal_coefficient_n_positive, const std::vector<complex_type>& rhsv_minimal_coefficient_n_negative, const float_type& z, complex_type& R_up);
    // Computes the up MST radial solution (not renormalized)
    void rhsm_compute_Rup_and_1der(const complex_type& nu, const std::vector<complex_type>& minimal_coefficient_n_positive, const std::vector<complex_type>& minimal_coefficient_n_negative, const float_type& z, complex_type& R_up, complex_type& R_up_der);
    // Computes the up solution (not renormalized) and its first derivative

    void rhsm_compute_Rup_and_1der_and_2der(const complex_type& nu, const std::vector<complex_type>& minimal_coefficient_n_positive, const std::vector<complex_type>& minimal_coefficient_n_negative, const float_type& z, complex_type& R_up, complex_type& R_up_der, complex_type& R_up_2der);
    // Computes the up solution (not renormalized) and its first and second derivatives

    void rhsm_compute_R_in_mst(const float_type& r, complex_type& R_in);
    // Computes the in homogeneous solution (renormalized) using MST expansion 

    void rhsm_compute_R_in_and_1der_mst(const float_type& r, complex_type& R_in,complex_type& R_in_der);
    // Computes the in solution (renormalized) and its first derivative (MST)

    void rhsm_compute_R_in_and_1der_and_2der_mst(const float_type& r, complex_type& R_in,complex_type& R_in_der, complex_type& R_in_2der);
    // Computes the in solution (renormalized) and its first and second derivatives (MST)

    void rhsm_compute_R_up_mst(const float_type& r, complex_type& R_up);
    // Computes the up homogeneous solution (renormalized) using MST

    void rhsm_compute_R_up_and_1der_mst(const float_type& r, complex_type& R_up, complex_type& R_up_der);
    // Computes the up solution (renormalized) and its first derivative

    void rhsm_compute_R_up_and_1der_and_2der_mst(const float_type& r, complex_type& R_up, complex_type& R_up_der, complex_type& R_up_2der);
    // Computes the up solution (renormalized) and its first and second derivatives
    
    void rhsm_compute_coefficients_of_eq_numerical_integration(const int H, const float_type& r, complex_type& y2der, complex_type& y2, complex_type& y1);
    // Computes coefficients of the radial ODE used in numerical integration
    
    void rhsm_compute_R_in_numerical_integration(const float_type& r_start_integration, const float_type& r_end_integration);
    // Numerically integrates the in solution between two radii
    
    void rhsm_compute_R_up_numerical_integration(const float_type& r_start_integration, const float_type& r_end_integration);
    // Numerically integrates the up solution between two radii
    
public:

    RadialHomogeneousSolution(int s, float_type M, float_type a, float_type m, float_type lambda, float_type omega, int ll, int mm);
    // Constructor: initializes parameters and computes ν and normalization constants

    void rhsm_get_R_in_numerical_integration(const float_type& r_start_integration, const float_type& r_end_integration, const float_type&r, complex_type& R_in);
    // Returns numerically integrated in solution at radius r

    void rhsm_get_R_up_numerical_integration(const float_type& r_start_integration, const float_type& r_end_integration, const float_type&r, complex_type& R_up);
    // Returns numerically integrated up solution at radius r

    void rhsm_get_R_in_numerical_integration_and_1der_and_2der(const float_type& r_start_integration, const float_type& r_end_integration, const float_type&r, complex_type& R_in, complex_type& R_in_der, complex_type& R_in_2der);
    // Returns numerically integrated in solution and its derivatives

    void rhsm_get_R_up_numerical_integration_and_1der_and_2der(const float_type& r_start_integration, const float_type& r_end_integration, const float_type&r, complex_type& R_up, complex_type& R_up_der, complex_type& R_up_2der);
    // Returns numerically integrated up solution and its derivatives

    void rhsm_get_nu(complex_type& res);
    // Returns the renormalized angular momentum ν

    void rhsm_get_R_in_mst(const float_type& r,complex_type& res);
    // Returns the in solution computed via MST expansion

    void rhsm_get_R_up_mst(const float_type& r,complex_type& res);
    // Returns the up solution computed via MST expansion

    void rhsm_get_R_in_and_1der_mst(const float_type&r, complex_type& R_in, complex_type& R_in_der);
    // Returns the in solution and its derivative via MST

    void rhsm_get_R_up_and_1der_mst(const float_type&r, complex_type& R_up, complex_type& R_up_der);
    // Returns the up solution and its derivative via MST

    void rhsm_get_R_in_and_1der_and_2der_mst(const float_type&r, complex_type& R_in, complex_type& R_in_der, complex_type& R_in_2der);
    // Returns the in solution and first two derivatives via MST

    void rhsm_get_R_up_and_1der_and_2der_mst(const float_type&r, complex_type& R_up, complex_type& R_up_der, complex_type& R_up_2der);
    // Returns the up solution and first two derivatives via MST

    void rhsm_get_B_trans(complex_type& res);
    // Returns transmission amplitude at the horizon

    void rhsm_get_C_trans(complex_type& res);
    // Returns transmission amplitude of the up solution

    void rhsm_get_B_inc(complex_type& res);
    // Returns incoming amplitude at infinity

    void rhsm_get_B_ref(complex_type& res);
    // Returns reflected amplitude at infinity
    
    ~RadialHomogeneousSolution() 
    {
    }
    
};

#endif
#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/math/tools/fraction.hpp>
#include <boost/math/differentiation/finite_difference.hpp>
#include <omp.h>
#include "config.h"
#include "special_functions.h"

#ifndef SPINWEIGHTEDHARMONICS_H
#define SPINWEIGHTEDHARMONICS_H

// REFERENCES: 
//  - Fujita-Tagoshi, Prog.Theor.Phys. 112 (2004) 415-450 (for low-frequency expansion of eigenvalue and eigenvalue solution via continued fraction method)
//  - Leaver, Proc.Roy.Soc.Lond.A 402 (1985) 285-298 (for eigenfunction solution via Leaver series method)
//  - Hughes, Phys.Rev.D 61 (2000) 8, 084004 (for spectral method solution of eigenvalue)
//  -  M. Casals, A. C. Ottewill, and N. Warburton, Proc. Roy. Soc. Lond. A 475, 20180701 (2019) (for high-frequency expansion of eigenvalue)

// Class that computes spin-weighted spheroidal harmonics and their eigenvalues
class SpinWeightedSpheroidalHarmonics 
{

  // ===== Input parameters =====

  int swshv_s;           // spin weight s
  int swshv_l;           // spherical harmonic index l
  int swshv_m;           // azimuthal harmonic index m
  float_type swshv_zeta; // oblateness parameter ζ = aω

  // ===== Auxiliary coefficients used in angular equation =====

  int swshv_alpha;       // auxiliary parameter α
  int swshv_beta;        // auxiliary parameter β
  float_type swshv_k1;   // auxiliary coefficient k1
  float_type swshv_k2;   // auxiliary coefficient k2

  // ===== Eigenvalue quantities =====

  float_type swshv_eigenvalue; // angular eigenvalue E
  float_type swshv_lambda;     // separation constant λ

  // ===== Leaver series variables =====

  std::vector<float_type> swshv_an_leaver; // reries coefficients a_n of Leaver expansion
  complex_type swshv_normalization_leaver; // normalization constant of the eigenfunction
  int swshv_sign_leaver;                   // sign convention used for normalization

  // ===== Harmonic values =====

  complex_type swshv_S;              // spin Weighted Spheroidal Harmonic S(θ,φ)
  complex_type swshv_S_prime;        // first derivative dS/dθ
  complex_type swshv_S_double_prime; // second derivative d²S/dθ²

  // Continued fraction helper structure
  struct swhs_continued_fraction_Rn;

  // ===== Matrix utilities for spectral method =====

  void swshm_find_off_diagonal_max(const std::vector<std::vector<float_type>>& A, int& index_row, int& index_column, float_type& maxVal);
  // Find largest off-diagonal matrix element

  void swshm_construct_identity_matrix(std::vector<std::vector<float_type>>& J);
  // Construct identity matrix

  void swshm_multiply_matrices(const std::vector<std::vector<float_type>>& A, const std::vector<std::vector<float_type>>& B, std::vector<std::vector<float_type>>& result);
  // Matrix multiplication

  void swshm_transpose_matrix(const std::vector<std::vector<float_type>>& A, std::vector<std::vector<float_type>>& result); 
  // Matrix transpose

  void swshm_jacobi_eigenvalues(std::vector<std::vector<float_type>>& A, std::vector<float_type>& eigenvalues);
  // Compute eigenvalues using Jacobi diagonalization


  // ===== Spectral method eigenvalue computation =====

  void kHat(const float_type& l, float_type& result);      // Compute coefficient k̂(l)
  void k2(const float_type& l, float_type& result);        // Compute coefficient k²(l)
  void kTilde2(const float_type& l, float_type& result);   // Compute coefficient k̃²(l)

  void swshm_compute_eigenvalue_with_spectral_method(float_type& eigenvalue);
  // Compute eigenvalue using spectral expansion


  // ===== Analytical approximations for eigenvalue =====

  void swshm_compute_eigenvalue_low_frequency_expansion(float_type& eigenvalue);
  // Low-frequency expansion (|ζ| << 1)

  void swshm_compute_eigenvalue_high_frequency_expansion(float_type& eigenvalue);
  // High-frequency expansion (|ζ| >> 1)


  // ===== Continued fraction method =====

  void swshm_compute_alpha_n(const int n, float_type& alpha_n);
  // Recurrence coefficient α_n

  void swshm_compute_beta_n(const int n, const float_type& e, float_type& beta_n);
  // Recurrence coefficient β_n

  void swshm_compute_gamma_n(const int n, float_type& gamma_n);
  // Recurrence coefficient γ_n

  void swshm_compute_Rn(const int n, const float_type& E, float_type& Rn); 
  // Compute continued fraction ratio R_n

  void swshm_compute_Ln(const int n, const float_type& E, float_type& Ln); 
  // Compute continued fraction ratio L_n

  void swshm_compute_function_eigenvalue(const float_type& e, float_type& f_e); 
  // Eigenvalue equation f(E)

  void swshm_compute_derivative_function_eigenvalue(const float_type& e, float_type& f_primec);
  // Derivative of eigenvalue equation

  void swshm_find_root_for_eigenvalue();
  // Newton method to find eigenvalue root


  // ===== Leaver series method for eigenfunction =====

  void shwm_compute_alpha_n_leaver(const int n, float_type& a_n);
  // Compute α_n coefficient of Leaver recurrence

  void shwm_compute_beta_n_leaver(const int n, float_type& beta_n);
  // Compute β_n coefficient of Leaver recurrence

  void shwm_compute_beta_n_leaver_for_eigenvalue(const int n, const float_type& e, float_type& beta_n);
  // β_n used when eigenvalue is unknown

  void shwm_compute_gamma_n_leaver(const int n, float_type& gamma_n);
  // Compute γ_n coefficient of Leaver recurrence

  void shwm_compute_a_n_vector_leaver();
  // Compute vector of Leaver expansion coefficients

  void shwm_compute_normalization_vector_leaver();
  // Compute normalization constant

  void shwm_compute_sign_leaver();
  // Determine sign convention of the harmonic

  void shwm_compute_S(const float_type& theta, const float_type& phi);
  // Compute harmonic S(θ,φ)

  void shwm_compute_S_dS_ddS(const float_type& theta, const float_type& phi);
  // Compute S and its first and second derivatives


public:

  // Constructor initializing (s,l,m,ζ)
  SpinWeightedSpheroidalHarmonics(int ss, int ll, int mm, float_type obliteness);

  // Return angular eigenvalue
  void swshm_get_eigenvalue(float_type& eig);

  // Return separation constant λ
  void swshm_get_lambda(float_type& lambda);

  // Compute harmonic S(θ,φ)
  void swshm_get_S(const float_type& x, const float_type& phi, complex_type& S);

  // Compute S, dS, and d²S
  void swshm_get_S_dS_ddS(const float_type& x, const float_type& phi, complex_type& S, complex_type& dS, complex_type& ddS);

  // Destructor
  ~SpinWeightedSpheroidalHarmonics() 
  {
  }
};

#endif
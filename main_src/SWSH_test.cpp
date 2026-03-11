#include "SpinWeightedSpheroidalHarmonics.h"

int main()
{ 
  // Set numerical precision for output
  std::cout << std::setprecision(20) << std::endl;
 
  int s = 0;   // Spin weight
  int l = 1;   // Harmonic degree
  int m = 1;   // Azimuthal number

  std::cout << "s = " << s << std::endl;
  std::cout << "l = " << l << std::endl;
  std::cout << "m = " << m << std::endl;
  
  // ===== Physical parameters =====
  float_type a("0.8");      // Kerr spin parameter
  float_type omega("6.7");  // Mode frequency

  // Oblateness parameter ζ = aω
  float_type oblateness = a*omega;
  std::cout << "oblateness = " << oblateness << std::endl;

  // Construct SpinWeightedSpheroidalHarmonics object
  SpinWeightedSpheroidalHarmonics swsh_obj = SpinWeightedSpheroidalHarmonics(s, l, m, oblateness);
  
  // ===== Compute eigenvalue =====
  float_type eigenvalue;                      // Angular eigenvalue
  swsh_obj.swshm_get_eigenvalue(eigenvalue); // Retrieve eigenvalue
  std::cout << "eigenvalue = " << eigenvalue << std::endl;

  // ===== Compute separation constant λ =====
  float_type lambda;                          // Separation constant
  swsh_obj.swshm_get_lambda(lambda);         // Retrieve λ
  std::cout << "lambda = " << lambda << std::endl;
  
  // ===== Evaluate SWSH at given angular position =====
  float_type theta = 0.3*constants::pi;       // Polar angle θ
  float_type x = mp::cos(theta);              // x = cos(θ)
  float_type phi = 0;                         // Azimuthal angle φ
  
  // ===== Compute SWSH and derivatives =====
  complex_type S;                             // SWSH S(θ,φ)
  complex_type dS;                            // First derivative dS/dθ
  complex_type ddS;                           // Second derivative d²S/dθ²

  swsh_obj.swshm_get_S_dS_ddS(x, phi, S, dS, ddS); // Compute S, dS, d²S

  std::cout << "S = " << S << std::endl;
  std::cout << "dS = " << dS << std::endl;
  std::cout << "ddS = " << ddS << std::endl;

  return(0);
}
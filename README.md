<p align="center">
  <img src="STORM.png" width="500">
</p>

# STORM

This repository contains a C++ library for computing scalar fluxes on generic orbits around Kerr black holes. 
It is designed for high-precision computations relevant to Extreme Mass-Ratio Inspirals (EMRIs) and uses Boost Multiprecision to ensure numerical accuracy.

The code provides:

- **SpinWeightedSpheroidalHarmonics**: class that computes spin-weighted spheroidal harmonics.
- **RadialHomogeneousSolution**: class that computes homogeneous solutions of the radial Teukolsky equation. 
- **GeodesicOrbitalMotion**: class that computes quantities related to the geodesic orbital motion of a test particle around Kerr black holes.
- **ScalarFluxMode**: class that computes scalar energy fluxes from a particle in a generic orbit around a Kerr black hole.
- **special_functions**: Namespace providing high-precision implementations of special functions (gamma functions, hypergeometric functions) for use with Boost multiprecision types.
- **config.h**: Global configuration file where numerical precision, tolerances, and algorithm parameters can be set for the entire code.

---

## Dependencies

The project requires:

- **C++17 compiler** (verified with g++, should work with any C++17-compliant compiler)
- **Boost library (version 1.85 or higher)**, including the following components:
  - **Multiprecision types** (`cpp_complex`, `cpp_bin_float`) for high-precision floating-point and complex arithmetic.
  - **Mathematical constants** (`boost::math::constants`) for mathematical constants. 
  - **Continued fraction** (`boost::math::tools::fraction`) for continued fraction computation. 
  - **Quadrature routines** (`boost::math::quadrature`) for numerical integration.
  - **Special functions** (`boost::math::special_functions`) for special functions evaluation.
  - **ODE integration and root-finding utilities** (`boost::numeric::odeint`) for ODE integration.
  - **Finite difference utilities** (`boost::math::differentiation`) for numerical derivative computation.

The Boost library is publicly available [here](https://www.boost.org/).
Make sure Boost is installed and update the `INCLUDES` and `BOOST_LIB_DIR` paths in the Makefile.

---

## Installation / Compilation

This project uses a Makefile to compile the library and example programs.
To Compile all classes and examples:

```bash
make
```
This will generate:

- `RHS_test` – executable associated to the example of RadialHomogeneousSolution

- `SWSH_test`  → executable associated to the example of SpinWeightedSpheroidalHarmonics

- `GOM_test`  → executable associated to the example of GeodesicOrbitalMotion

- `SFM_test`  → executable associated to the example of ScalarFluxMode

To compile individual targets:

```bash
make RHSTEST
make GOMTEST
make SWSHTEST
make SFMTEST
```

Clean generated files:
```bash
make clean
```

## Directory Structure

```
.
├── include/                  # Header files for all classes
├── src/                      # Source files for all classes
├── main_src/                 # Example programs
├── obj/                      # Object files (generated)
├── exe/                      # Executables (generated)
├── Makefile                  # Build instructions
└── README.md                 # This file
```

---

## Classes and Examples

1. **RadialHomogeneousSolution**  
    - *Header file*: `include/RadialHomogeneousSolution.h`  
    - *Source file*: `src/RadialHomogeneousSolution.cpp`  
    - *Example usage*: `main_src/RHS_test.cpp`  

    - *Purpose:* Computes radial homogeneous solutions of the Teukolsky equation for arbitrary spin weight `s`.  
      The class provides:  
        - Renormalized angular momentum (computed using both the monodromy method and continued fraction method from Fujita & Tagoshi)  
        - Asymptotic amplitudes at infinity and at the horizon (MST method)
        - Radial solutions at the horizon and at infinity (MST and numerical integration methods)
        - First and second derivatives of the radial solution (MST and numerical integration methods)

2. **SpinWeightedSpheroidalHarmonics**  
    - *Header file*: `include/SpinWeightedSpheroidalHarmonics.h`  
    - *Source file*: `src/SpinWeightedSpheroidalHarmonics.cpp`  
    - *Example usage*: `main_src/SWSH_test.cpp`  

    - *Purpose:* Computes spin-weighted spheroidal harmonics and their associated eigenvalues for arbitrary spin weight `s`.  

      The class provides:  
        - Eigenvalues of the spin-weighted spheroidal harmonics:
          - Computed using the continued fraction method (Fujita & Tagoshi)  
          - Initial guess obtained from the spectral method 
          - Optionally, low- and high-oblateness expansions are also available  
        - Eigenfunctions:
          - Computed using Leaver’s method for high accuracy  



3. **GeodesicOrbitalMotion**  
    - *Header file*: `include/GeodesicOrbitalMotion.h`  
    - *Source file*: `src/GeodesicOrbitalMotion.cpp`  
    - *Example usage*: `main_src/GOM_test.cpp`  

    - *Purpose:* Computes the orbital motion of a particle around a Kerr black hole on generic orbits (inclined and eccentric).  

      The class provides:  
        - Constants of motion: energy, angular momentum, and Carter constant.
        - Orbital frequencies: Mino time frequencies, coordinate time frequencies, angular frequency.
        - Geodesic solution: radial, polar, coordinate time, and azimuthal angle functions.

4. **ScalarFluxMode**  
    - *Header file*: `include/ScalarFluxMode.h`  
    - *Source file*: `src/ScalarFluxMode.cpp`  
    - *Example usage*: `main_src/SFM_test.cpp`  

    - *Purpose:* Computes scalar energy fluxes emitted by a particle orbiting a Kerr black hole, both at the horizon and at infinity.  

      The class supports:  
        - Orbit types:
          - Circular orbits  
          - Eccentric orbits  
          - Generic orbits  
          - Energy fluxes: 
          - Horizon flux  
          - Infinity flux  

    - *Notes on computation:*
    - Flux computations are currently performed on a single CPU.  
    - Parallelization over multiple flux calculations (e.g., for different modes or orbital parameters) is planned for future versions.

All example programs demonstrate how to:

- **Initialize physical parameters of the system:**
  - Black hole mass and spin
  - Particle mass
  - Orbital parameters (semi-latus rectum, eccentricity, inclination)
  -  Mode indices

- **Create class objects and compute quantities**

- **Retrieve results using the provided getter functions**

- **Print or further process the outputs**


## Compilation Notes
- Numerical tolerances, precision, and algorithm parameters are configurable in `config.h`.

## TODO / Future Work
- Add a class for gravitational fluxes (almost done, currently in testing phase).  
- Implement static modes (i.e., ω = 0) in the raRadialHomogeneousSolutiondial class (almost done, currently in testing phase).  
- Implement the computation of the separatrix.  
- Add routines to compute angular momentum and Carter constant fluxes.  
- Create a utility/file to perform the sum over modes for circular, eccentric, and generic orbits.

## Authors / Contributors

- **Sara Gliorio** – Main developer and maintainer. For questions or issues, contact: <sara.gliorio@gssi.it>
- **Matteo Della Rocca** – Contributions to the `makefile` and the `SpinWeightedSpheroidalHarmonics` class
- **Next collaborator:** maybe you!


## Citation

If you use this code in your research, please cite:

- S. Gliorio et al., *Adiabatic evolution of asymmetric binaries on generic orbits with new fundamental fields I: characterization of gravitational wave fluxes*, [arXiv:2603.2603.10116](https://arxiv.org/abs/2603.10116)
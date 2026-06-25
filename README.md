# NAIAD

This is a one-dimensional SN neutron transport solver for Cartesian geometry. The main idea of this work was to investigate anisotropic scatteirng in calculations with high leakage.

# Usage

To build.

```
git clone <url>
cmake -B build
cmake --build build
```

Naiad uses CMake to build the code. You can use typical CMake configuration options.

Note that Naiad uses OpenMP (shared-memory) parallelism with C++. I have developed this on a Mac laptop, and the clang++ compiler installed with MacOS does not natively support OpenMP.

If you too are using a Mac, my recommended installation is to use [homebrew](https://brew.sh/) to install a version of `llvm`.

# Features

- Ready now:
    - Many energy groups.
    - High-order SN (up to S64).
    - Parallel calculations with OpenMP.
    - Anisotropic scattering.
- Work in progress:
    - Many spatial discretization schemes (discrete ordinates, step characteristics, linear characteristics, quadratic characteristics, and discontinuous Galerkin).

# Name

This project is a successor to my other project [Siren](https://github.com/wcdawn/siren). The naiad is another ancient seductress in the waters of ancient Greece.

# Support

I have recently learned that some people actually read these README pages. If you have read this and have any questions about this project at all, please email me. My email address is not hard to find. My full C.V. is on my [website](https://wcdawn.github.io/).

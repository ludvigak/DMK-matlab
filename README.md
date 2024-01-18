# DMK Matlab
Matlab implementation of the Dual-space Multilevel Kernel-splitting (DMK) method

## Setup

To initialize submodules, run
```
git submodule update --init --recursive
```

## Usage

* Start Matlab in root directory
* `init` initializes paths
* `runtests` runs the test suite
* `demo_dmk` runs DMK and compares it to a direct summation
* `demo_scaling` shows linear scaling
*  `show_interp_error` analyzes the interpolation error in the difference kernel
* `dev_*` files are bits and pieces of code used when developing functionality
  * `dev_ewald_sum` has code for developing and checking Ewald splits

## Code structure

* `external/` contains third-party packages in the form of Git submodules
* `src/` contains all routines needed for DMKL
* `src/+kernels/` contains kernel splittings
* `src/+approx/` is a package with the routines for     approximation using Chebyshev functions and tensor product grids. Holds its own tests in `tests.m`
* `test/` contains all unit tests

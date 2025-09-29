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
* `demo_dmk` runs DMK and compares it to a direct summation. Showcases use for all kernels and free space / periodic.
* `demo_scaling` shows linear scaling

## Code structure

* `dev_scripts/`: bits and pieces of code used when developing functionality
* `experiments/`: experiments reported in paper
* `external/`: third-party packages in the form of Git submodules
* `src/`: all routines needed for DMK
* `src/+kernels/`: kernel splittings
* `src/+approx/` is a package with the routines for approximation using Chebyshev functions and tensor product grids. Holds its own tests in `tests.m`
* `test/`: unit tests


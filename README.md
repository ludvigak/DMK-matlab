# DMK Matlab
Matlab implementation of the Dual-space Multilevel Kernel-splitting (DMK) method for Stokes and Laplace potentials.

Companion code to the preprint [Fast summation of Stokes potentials using a new kernel-splitting in the DMK framework](https://arxiv.org/abs/2509.21471) by Ludvig af Klinteberg, Leslie Greengard, Shidong Jiang, and Anna-Karin Tornberg.

This is a prototype implementation of the DMK algorithm, designed to aid in understanding the algorithm and in the design of new kernels or other modifications. It is *not* a fast HPC code. For that, please refer to the official implementation hosted by the Flatiron Institute at https://github.com/flatironinstitute/dmk

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
* Scripts in `experiments/` reproduce results from paper.
  
## Code structure

* `dev_scripts/`: bits and pieces of code used when developing functionality
* `external/`: third-party packages in the form of Git submodules
* `src/`: all routines needed for DMK
* `src/+kernels/`: kernel splittings
* `src/+approx/` is a package with the routines for approximation using Chebyshev functions and tensor product grids. Holds its own tests in `tests.m`
* `test/`: unit tests


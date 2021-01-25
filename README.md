# Adaptive Stagger-Tilted Grid

Implementation of paper **An adaptive staggered-tilted grid for incompressible flow simulation** ([doi](https://dl.acm.org/doi/10.1145/3414685.3417837))

NOTE: This repo only contains 2d implementation and is still under development.

## Build

1. download dependency
    ```
    git submodule update --init --recursive
    ```
2. build
   ```
   mkdir build & cd build
   cmake -DCMAKE_BUILD_TYPE=RELEASE ../
   make -j4
   ```

The above process will generated a shared library under directory *./python*. 

## Run

Under the *python* directory, run `jupyter-lab` (requires installation of jupyter and related python packages).

The default **ast_smoke_2d** notebook will run a simple 2d smoke simulation and visualize (and save) the simulation result.
# stochasticQP

## Building with CMake on ManeFrame II (M2)

```
module purge
module load cmake/3.19.1 gcc-9.2 gurobi/8.0.1
git clone git@github.com:SMU-SODA/stochasticQP.git
cd stochasticQP
mkdir build
cd build
cmake ..
make
```


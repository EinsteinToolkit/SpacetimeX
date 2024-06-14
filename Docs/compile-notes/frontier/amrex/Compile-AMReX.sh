#!/bin/bash

cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_INSTALL_PREFIX=${HOME}/local/amrex-24.06-cray \
      -DCMAKE_PREFIX_PATH='/opt/rocm-6.0.0/lib/cmake/AMDDeviceLibs;/opt/rocm-6.0.0/lib/cmake/amd_comgr;/opt/rocm-6.0.0/lib/cmake/hip;/opt/rocm-6.0.0/lib/cmake/hiprand;/opt/rocm-6.0.0/lib/cmake/hsa-runtime64;/opt/rocm-6.0.0/lib/cmake/rocprim;/opt/rocm/lib/cmake/rocrand' \
      -DAMReX_GPU_BACKEND=HIP \
      -DAMReX_AMD_ARCH=gfx90a \
      -DAMReX_FORTRAN=OFF \
      -DAMReX_FORTRAN_INTERFACES=OFF \
      -DAMReX_OMP=OFF \
      -DAMReX_PARTICLES=ON \
      -DAMReX_PRECISION=DOUBLE \
      ..
      # -DCMAKE_CXX_COMPILER=${ROCM_PATH}/bin/amdclang++ \
      # -DCMAKE_C_COMPILER=${ROCM_PATH}/bin/amdclang \

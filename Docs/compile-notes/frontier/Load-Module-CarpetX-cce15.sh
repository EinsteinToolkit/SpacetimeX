#!/bin/bash

module load craype-accel-amd-gfx90a
module load rocm/5.3.0

export MPICH_GPU_SUPPORT_ENABLED=1
export PE_MPICH_GTL_DIR_amd_gfx90a="-L/opt/cray/pe/mpich/8.1.23/gtl/lib"
export PE_MPICH_GTL_LIBS_amd_gfx90a="-lmpi_gtl_hsa"

module load adios2/2.8.3
module load amrex/22.11
module load boost/1.79.0-cxx17
module load cray-fftw/3.3.10.3
module load cray-hdf5-parallel/1.12.2.1
module load gsl/2.7.1
module load hwloc/2.5.0
module load libjpeg-turbo/2.1.0
module load openblas/0.3.17
module load openpmd-api/0.14.4
module load zlib/1.2.11

#!/bin/bash

module purge
module load PrgEnv-cray/8.5.0
module load amd-mixed/6.0.0
module load cce/17.0.0
module load cpe/23.12
module load rocm/6.0.0
module load craype-accel-amd-gfx90a
module load craype-x86-trento
module unload darshan-runtime

################################################################################
# from AMReX-Astro
################################################################################
#module load PrgEnv-gnu
#module load cray-mpich/8.1.27
#module load craype-accel-amd-gfx90a
#module load amd-mixed/6.0.0
#module unload darshan-runtime

################################################################################
# from vass
################################################################################
#module load cce/17.0.0
#module load cpe/23.12
#module load rocm/6.0.0
#module load craype-accel-amd-gfx90a
#module load cray-hdf5-parallel
#module load cray-fftw
#module load boost
#module load adios2
#module load hsi
#module load gsl
#module load lfs-wrapper
#module load cray-dsmml/0.2.2


################################################################################
# from simfactory
################################################################################
#module load craype-accel-amd-gfx90a
#module load rocm/5.3.0
#
#export MPICH_GPU_SUPPORT_ENABLED=1
#export PE_MPICH_GTL_DIR_amd_gfx90a="-L/opt/cray/pe/mpich/8.1.23/gtl/lib"
#export PE_MPICH_GTL_LIBS_amd_gfx90a="-lmpi_gtl_hsa"
#
#module load adios2/2.8.3
#module load amrex/22.11
#module load boost/1.79.0-cxx17
#module load cray-fftw/3.3.10.3
#module load cray-hdf5-parallel/1.12.2.1
#module load gsl/2.7.1
#module load hwloc/2.5.0
#module load libjpeg-turbo/2.1.0
#module load openblas/0.3.17
#module load openpmd-api/0.14.4
#module load zlib/1.2.11


################################################################################
# Old module list
################################################################################

#module load craype-x86-trento
#module load libfabric/1.15.2.0
#module load craype-network-ofi
#module load perftools-base/22.12.0
#module load xpmem/2.6.2-2.5_2.22__gd067c3f.shasta
#module load cray-pmi/6.1.8
#module load DefApps/default
#module load ccache/4.5.1
#module load cmake/3.23.2
#module load craype-accel-amd-gfx90a
##module load rocm/5.2.0
#module load rocm/5.3.0
#module load craype/2.7.19
#module load cray-dsmml/0.2.2
#module load cray-libsci/22.12.1.1
#module load PrgEnv-gnu/8.3.3

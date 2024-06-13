#!/bin/bash

module purge
module load PrgEnv-cray/8.5.0
module load cce/17.0.0
module load cpe/23.12
module load rocm/6.0.0
module load craype-accel-amd-gfx90a
module load craype-x86-trento
module load openpmd-api
module load boost
module load cray-fftw
module load cray-hdf5-parallel
module load gsl
module load hwloc
module load openblas
module load zlib
module load adios2
module load libjpeg-turbo
module unload darshan-runtime


################################################################################
# old 
################################################################################

#module purge
#module load DefApps/default
#module load libfabric/1.15.2.0
#module load craype-x86-trento
#module load craype-network-ofi
#module load cray-pmi/6.1.8
#module load craype-accel-amd-gfx90a
#module load cce/15.0.0
#module load craype/2.7.19
#module load cray-dsmml/0.2.2
#module load cray-mpich/8.1.23
#module load cray-libsci/22.12.1.1
#
#module load perftools-base/22.12.0
#module load xpmem/2.6.2-2.5_2.22__gd067c3f.shasta

#!/bin/bash

module purge
module load PrgEnv-cray/8.5.0
module load amd-mixed/6.0.0
module load cpe/23.12
module load craype-accel-amd-gfx90a
module load craype-x86-trento
module unload darshan-runtime

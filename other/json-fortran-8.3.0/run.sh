#!/bin/bash

#export OMP_PROC_BIND=true
#export OMP_PLACES=cores
#export OMP_NUM_THREADS=1
#export SCOREP_METRIC_PAPI=PAPI_L2_DCM

# Load modules
#ml papi

# CMake
cd /build/
#cmake ..
#scorep-wrapper --create gfortran .
#export PATH="/cfs/klemming/home/e/ejacobss/Private/Project_Code/build:$PATH"
SCOREP_WRAPPER=off cmake .. -DCMAKE_Fortran_COMPILER=scorep-gfortran

# Make
make

# Run code
cd /cfs/klemming/home/e/ejacobss/Private/Project_Code/build/src/main
scan srun -n1 ./mainout
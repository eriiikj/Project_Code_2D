#!/bin/bash

#export OMP_PROC_BIND=true
#export OMP_PLACES=cores
#export OMP_NUM_THREADS=1
#export SCOREP_METRIC_PAPI=PAPI_L2_DCM

# Load modules
# ml PDC/21.11
# ml CMake/3.21.2
# ml PrgEnv-gnu
# ml Scalasca/2.6-cpeGNU-21.11
# ml papi


export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_NUM_THREADS=1

# -------- Normal make and execution ----


# # CMake
# cd /home/er7128ja/Nextcloud/Projekt/Project_Code
# # cmake .. -DCMAKE_Fortran_COMPILER=gfortran
# cmake .. -DCMAKE_Fortran_COMPILER=ifort

# Make
cd /home/er7128ja/Nextcloud/Projekt/Project_Code_2D/build
make

# Run code
cd /home/er7128ja/Nextcloud/Projekt/Project_Code_2D/build/src/main
./mainout
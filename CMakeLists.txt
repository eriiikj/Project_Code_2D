# ----- General setup -----
cmake_minimum_required(VERSION 3.13)

# Specify compiler
set(CMAKE_Fortran_COMPILER "ifort")
# set(CMAKE_Fortran_COMPILER "gfortran")

# OpenMP
set(omp_run FALSE)

# Get compiler name
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# Make sure that the default build type is RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

message("Build type:")
message("${CMAKE_BUILD_TYPE}")

message("Fortran compiler:")
message("${CMAKE_Fortran_COMPILER}")

# ----- Project Setup -----

# Set the project name
project(mainout)

# Add executable
# add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/main/main.f90) 

# Enable Fortran language
enable_language(Fortran)

# Add OpenMP
find_package(OpenMP REQUIRED)
message("CMAKE_FC_FLAGS:")
message("${CMAKE_FC_FLAGS}")
message("${OPENMP_FOUND}")

# ----- Add static libs -----
if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  set(jsonlib_path ${PROJECT_SOURCE_DIR}/src/staticlib/json_gfortran)
elseif (Fortran_COMPILER_NAME STREQUAL "scorep-gfortran")
  set(jsonlib_path ${PROJECT_SOURCE_DIR}/src/staticlib/json_gfortran_scorep)
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  set(jsonlib_path ${PROJECT_SOURCE_DIR}/src/staticlib/json_intel)
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

add_library(jsonlib STATIC IMPORTED)
set_target_properties(jsonlib PROPERTIES IMPORTED_LOCATION ${jsonlib_path}/lib/libjsonfortran.a)

message("jsonlib_path:")
message("${jsonlib_path}")

# Add static lib mods
set(jsonlib_mods ${jsonlib_path}/include)

# ----- Add own libraries -----
# Specify Somelib version, somelib_ifort or somelib_gfortran. Note that ifort can run both version, whereas gfortran only runs with somelib_gfortran
if ((Fortran_COMPILER_NAME STREQUAL "gfortran")  OR (Fortran_COMPILER_NAME STREQUAL "scorep-gfortran"))
  set(somelib_variant somelib_gfortran)
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  set(somelib_variant somelib_ifort)
endif((Fortran_COMPILER_NAME STREQUAL "gfortran") OR (Fortran_COMPILER_NAME STREQUAL "scorep-gfortran"))
message("Somelib version:")
message("${somelib_variant}")

# Set paths to mod files being built in binary tree
set(mflib_mods        ${PROJECT_BINARY_DIR}/src/lib/mflib)
set(somelib_mods      ${PROJECT_BINARY_DIR}/src/lib/${somelib_variant})
set(whiskerlib_mods   ${PROJECT_BINARY_DIR}/src/modules/whiskermodules)


# Add libs in subfolders, as well as executable
add_subdirectory(src)

# ----- Add imported libraries -----

# Set mkl lib paths
set(MKLroot /software/Intel/XE2018/mkl/lib/intel64)

# Add mkl lib
if (Fortran_COMPILER_NAME STREQUAL "ifort")
    set(mkllib1_path ${MKLroot}/libmkl_intel_lp64.a)
    if (omp_run)
      set(mkllib2_path ${MKLroot}/libmkl_intel_thread.a)
    else()
      set(mkllib2_path ${MKLroot}/libmkl_sequential.a)
    endif (omp_run)
  elseif ((Fortran_COMPILER_NAME STREQUAL "gfortran") OR (Fortran_COMPILER_NAME STREQUAL "scorep-gfortran"))
    set(mkllib1_path ${MKLroot}/libmkl_gf_lp64.a)
    if (omp_run)
      set(mkllib2_path ${MKLroot}/libmkl_gnu_thread.a)
    else()
      set(mkllib2_path ${MKLroot}/libmkl_sequential.a)
    endif (omp_run)
endif (Fortran_COMPILER_NAME STREQUAL "ifort")

set(mkllib3_path ${MKLroot}/libmkl_core.a)

# Add imported mkl libs to workspace
add_library(mkllib1 STATIC IMPORTED)
add_library(mkllib2 STATIC IMPORTED)
add_library(mkllib3 STATIC IMPORTED)
set_target_properties(mkllib1 PROPERTIES IMPORTED_LOCATION ${mkllib1_path})
set_target_properties(mkllib2 PROPERTIES IMPORTED_LOCATION ${mkllib2_path})
set_target_properties(mkllib3 PROPERTIES IMPORTED_LOCATION ${mkllib3_path})

if (Fortran_COMPILER_NAME STREQUAL "ifort")
  if (omp_run)
    set(mkllib_link -Wl,--start-group mkllib1 mkllib2 mkllib3 -Wl,--end-group -liomp5 -lpthread -lm -ldl)
  else()
    set(mkllib_link -Wl,--start-group mkllib1 mkllib2 mkllib3 -Wl,--end-group -lpthread -lm -ldl)
  endif (omp_run)
elseif ((Fortran_COMPILER_NAME STREQUAL "gfortran") OR (Fortran_COMPILER_NAME STREQUAL "scorep-gfortran"))
  if (omp_run)
  set(mkllib_link -Wl,--start-group mkllib1 mkllib2 mkllib3 -Wl,--end-group -lgomp -lpthread -lm -ldl)
  else()
  set(mkllib_link -Wl,--start-group mkllib1 mkllib2 mkllib3 -Wl,--end-group -lpthread -lm -ldl)
  endif (omp_run)
endif(Fortran_COMPILER_NAME STREQUAL "ifort")


# Add mkl include root
#set(mkl_include_root /software/Intel/XE2018/compilers_and_libraries_2018.1.163/linux/mkl/include)


# Include matlab link
set(matlib_link -L/software/matlab/R2019b/bin/glnxa64  -leng -lmat -lmx -lut)
# set(matlib_link -L/usr/local/MATLAB/R2022b/bin/glnxa64 -leng -lmat -lmx -lut)


# ----- Link libraries to target -----

# Include location of mod files to target
target_include_directories(${PROJECT_NAME} PUBLIC ${mflib_mods} ${somelib_mods} ${whiskerlib_mods} ${jsonlib_mods})

# Link libraries to target
if (omp_run)
  target_link_libraries(${PROJECT_NAME} PUBLIC mflib whiskerlib somelib jsonlib ${matlib_link} ${mkllib_link} OpenMP::OpenMP_Fortran)
else()
  target_link_libraries(${PROJECT_NAME} PUBLIC mflib whiskerlib somelib jsonlib ${matlib_link} ${mkllib_link})
endif (omp_run)
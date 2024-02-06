add_library(mflib mf_datatypes.f90 mf_utils.f90)
set_target_properties(mflib PROPERTIES PREFIX "")
if (omp_run)
    target_link_libraries(mflib PUBLIC OpenMP::OpenMP_Fortran)
endif(omp_run)


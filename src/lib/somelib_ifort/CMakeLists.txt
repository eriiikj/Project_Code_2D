add_library(somelib memory_util.f90 abaqus_util.f90 elem_large_cont_2d.f90 fem_system.f90 fem_util.f90 mater_hyperel.f90 mater_J2iso.f90 mater_J2iso_Cu.f90 mater_J2iso_Sn.f90 mater_large.f90 matlab_util.f90 matrix_util.f90 memory_util.f90 some_constants.f90 sparse_util.f90 wrt2vtk.f90)
set_target_properties(somelib PROPERTIES PREFIX "")
if (omp_run)
    target_link_libraries(somelib PUBLIC OpenMP::OpenMP_Fortran)
endif(omp_run)

add_library(whiskerlib mesh_module.f90 grain_module.f90 imc_transf.f90 mater_hyperel_IMC.f90 mfc_util.f90 rhapsloop.f90 json_read.f90 level_set.f90 ls_utils.f90 ls_reconstruction_routines.f90 ls_sorting_routines.f90 ls_vp_routines.f90 ls_update_routines.f90 ls_types.f90 diffusion_types.f90 diffusion.f90 gp_util.f90 output_data.f90)
set_target_properties(whiskerlib PROPERTIES PREFIX "")

add_dependencies(whiskerlib mflib somelib)

# ----- Link whiskerlib to dependent libs ------
target_include_directories(whiskerlib PUBLIC ${mflib_mods} ${somelib_mods} ${jsonlib_mods})
if (omp_run)
    target_link_libraries(whiskerlib PUBLIC OpenMP::OpenMP_Fortran)
endif(omp_run)

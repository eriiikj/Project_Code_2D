module diffusion_types
    ! Last modified 
    ! E. Jacobsson 2023-11-28
      
    ! mf_datatypes
    use mf_datatypes

    ! Somelib
    use fem_system

    implicit none

    ! grain mesh system
    type grain_mesh_system
        integer               :: nalloc_nelm, nalloc_nnod, nalloc_nbc
        integer               :: nelm, nnod
        integer               :: dofnod, nelmx, nelmy, nrgp, nodel, ndof, dofel
        integer,  allocatable :: enod(:,:), bcnod(:,:), bcdof(:)
        real(dp), allocatable :: coord(:,:), bcval(:), bcval_idx(:)
        real(dp), allocatable :: a(:), r(:), ed(:,:), j_flux(:,:), jint(:)
        real(dp), allocatable :: p(:)
        
        ! Diffusion coefficient in grain
        real(dp)              :: M                
    end type grain_mesh_system


    ! Diffusion system
    type diffusion_system
        ! Save each grain_mesh_system object in an array called grain_meshes
        type(grain_mesh_system), allocatable :: grain_meshes(:)

        ! Thermodynamic parameters 3 d-20
        real(dp)              :: D_diff_material(3) = [(2.877d-36)*(1d3)**2d0*3600d0, (4.575d-20)*(1d3)**2d0*3600d0, &
        (2.452d-17)*(1d3)**2d0*3600d0]    
        real(dp)              :: thermo_parameterA(3)  = [1.0133d5, 4d5, 4.2059d6]         ! Parameter A [J/mol]
        real(dp)              :: thermo_parameterB(3)  = [-2.1146d4, -6.9892d3, 7.1680d3]  ! Parameter B [J/mol]
        real(dp)              :: thermo_parameterxb(3) = [0.10569d0, 0.41753d0, 0.99941d0] ! Parameter xbar [J/mol]
        real(dp)              :: molar_volumes(3)      = [7.09d0, 10.7d0, 16.3d0]*1d3      ! Molar volume [mm3/mol]

        ! Time
        real(dp)              :: generate_mesh_time=0d0, read_mesh_time=0d0, solve_time=0d0

        ! Global diffusion mesh
        integer,  allocatable :: enod(:,:)
        real(dp), allocatable :: coord(:,:)
        integer               :: nelm, nnod, nodel
        real(dp), allocatable :: ls(:,:), ls_ed(:,:,:)

    end type diffusion_system


end module diffusion_types
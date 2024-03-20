module output_data

  ! OpenMP
  use omp_lib

  ! mflib
  use mf_datatypes

  ! somelib
  use mater_J2iso_Cu
  use mater_J2iso_Sn
  use wrt2vtk
  use matlab_util

  ! whiskerlib
  use mesh_module
  use grain_module
  use gp_util
  use json_read
  use ls_types

  implicit none

  type plot_system
    real(dp), allocatable :: cau(:,:,:), sxx(:,:), syy(:,:), szz(:,:), sxy(:,:)    
    real(dp), allocatable :: biax(:,:), biax_Cu(:,:), biax_IMC(:,:), biax_Sn(:,:), biax_Sn2(:,:)
    real(dp), allocatable :: vm(:,:), vm_Cu(:,:), vm_IMC(:,:), vm_Sn(:,:)
    real(dp), allocatable :: alpha(:,:), alpha_Cu(:,:), alpha_Sn(:,:)
    real(dp), allocatable :: plastic(:,:), plastic_Cu(:,:), plastic_Sn(:,:), edd(:,:)
    real(dp), allocatable :: nodplot(:), nodplot2(:), nodplot_gp(:), E_save_gp(:,:)
    real(dp)              :: biax_Sn_avg
    real(dp), allocatable :: pgp(:,:), p_nod(:), p_ed(:,:) ! Hydrostatic pressure

    ! Plot gauss point mesh
    real(dp), allocatable :: gpx(:,:), gpy(:,:), coordgp(:,:)
    integer,  allocatable :: enodgp(:,:) 

    ! GP plot
    real(dp), allocatable :: gp_mat_label(:,:), gp_mat_label_nodplot(:), gp_phi_label(:,:), gp_phi_label_nodplot(:)
  end type plot_system

contains

  subroutine init_plot_system(pq, input_location, mesh)
    implicit none

    ! Intent inout
    type(plot_system), intent(inout)         :: pq
    
    ! Intent in
    character(len=:),allocatable, intent(in) :: input_location
    type(mesh_system), intent(in)            :: mesh

    ! Subroutine variables
    character(len=40)                        :: matfilename
    integer                                  :: ierr, nrgp, nelm, nnod, nnodgp

    ! Short form
    nrgp   = mesh%nrgp
    nelm   = mesh%nelm
    nnod   = mesh%nnod
    nnodgp = mesh%nnodgp

    ! Cauchy
    allocate(pq%cau(4,mesh%nrgp,mesh%nelm), stat=ierr)
    allocate(pq%sxx(nrgp,nelm), stat=ierr)
    allocate(pq%syy(nrgp,nelm), stat=ierr)
    allocate(pq%szz(nrgp,nelm), stat=ierr)
    allocate(pq%sxy(nrgp,nelm), stat=ierr)

    ! Hydrostatic pressure
    allocate(pq%pgp(nrgp,nelm), stat=ierr)
    allocate(pq%p_nod(nnod),stat=ierr)
    allocate(pq%p_ed(4,nelm),stat=ierr)

    ! Vm
    allocate(pq%vm(nrgp,nelm), stat=ierr)
    allocate(pq%vm_Cu(nrgp,nelm), stat=ierr)
    allocate(pq%vm_IMC(nrgp,nelm), stat=ierr)
    allocate(pq%vm_Sn(nrgp,nelm), stat=ierr)

    ! Biaxial
    allocate(pq%biax(nrgp,nelm), stat=ierr)
    allocate(pq%biax_Cu(nrgp,nelm), stat=ierr)
    allocate(pq%biax_IMC(nrgp,nelm), stat=ierr)
    allocate(pq%biax_Sn(nrgp,nelm), stat=ierr)
    allocate(pq%biax_Sn2(nrgp,nelm), stat=ierr)

    ! Plastic
    allocate(pq%alpha(nrgp,nelm), stat=ierr)
    allocate(pq%alpha_Cu(nrgp,nelm), stat=ierr)
    allocate(pq%alpha_Sn(nrgp,nelm), stat=ierr)
    allocate(pq%plastic(nrgp,nelm), stat=ierr)
    allocate(pq%plastic_Cu(nrgp,nelm), stat=ierr)
    allocate(pq%plastic_Sn(nrgp,nelm), stat=ierr)

    ! D avg
    allocate(pq%E_save_gp(nrgp,nelm), stat=ierr)

    ! Gp mesh
    allocate(pq%gpx(4,nelm), stat=ierr)
    allocate(pq%gpy(4,nelm), stat=ierr)
    allocate(pq%enodgp(4,nelm), stat=ierr)
    allocate(pq%coordgp(2,nnodgp))
    allocate(pq%gp_mat_label_nodplot(nnodgp), stat=ierr)    
    allocate(pq%gp_phi_label_nodplot(nnodgp), stat=ierr)

    ! Plot
    allocate(pq%nodplot(nnod), stat=ierr)
    allocate(pq%nodplot2(nnod), stat=ierr)
    allocate(pq%edd(nrgp,nelm), stat=ierr)
    allocate(pq%nodplot_gp(nnodgp), stat=ierr)    

    ! Gp plot    
    allocate(pq%gp_mat_label(nrgp,nelm), stat=ierr)
    allocate(pq%gp_phi_label(nrgp,nelm), stat=ierr)

    ! Change to input location
    call chdir(input_location)
    call chdir('mat_files')

    ! Filename
    write (matfilename, "(A,I0,A)") 'execution_data.mat'

    ! --- Write to file ---
    call matWrt2f(trim(matfilename), mesh%ex , 'ex', 'w')
    call matWrt2f(trim(matfilename), mesh%ey , 'ey', 'u')

    

    return
  end subroutine init_plot_system



  subroutine compute_plot_quantities(pq, mesh, lssys)
    implicit none

    ! Intent inout
    type(plot_system), intent(inout) :: pq
    type(mesh_system), intent(inout) :: mesh

    ! Intent in
    type(ls_system), intent(in)      :: lssys

    ! Subroutine variables
    integer :: i

    ! Cauchy stress
    pq%sxx = pq%cau(1,:,:)
    pq%syy = pq%cau(2,:,:)
    pq%szz = pq%cau(3,:,:)
    pq%sxy = pq%cau(4,:,:)
    
    ! Calculate von Mises stress
    pq%vm     = sqrt(1d0/2d0*((pq%sxx-pq%syy)**2d0+(pq%sxx-pq%szz)**2d0+(pq%szz-pq%syy)**2d0)+3d0*pq%sxy**2d0)
    pq%vm_Cu  = 0d0
    pq%vm_IMC = 0d0
    pq%vm_Sn  = 0d0
    where (lssys%hphi_gp(:,:,1).gt.0.99d0)                                     pq%vm_Cu  = pq%vm
    where (lssys%hphi_gp(:,:,1).lt.0.99d0 .or. lssys%hphi_gp(:,:,5).lt.0.99d0) pq%vm_IMC = pq%vm
    ! where (lssys%hphi_gp(:,:,5).gt.0.99d0)                                     pq%vm_Sn  = pq%vm
    ! pq%vm_Sn  = pq%vm*lssys%hphi_gp_plot(:,:,6) + pq%vm*lssys%hphi_gp_plot(:,:,7) + pq%vm*lssys%hphi_gp_plot(:,:,8)
    pq%vm_Sn  = pq%vm*lssys%sn_hphi_gp_plot


    ! Calculate biaxial stress
    pq%biax     = (pq%sxx + pq%szz)/2d0
    pq%biax_Cu  = 0d0
    pq%biax_IMC = 0d0
    pq%biax_Sn  = 0d0
    where (lssys%hphi_gp(:,:,1).gt.0.99d0)                                     pq%biax_Cu  = pq%biax
    where (lssys%hphi_gp(:,:,1).lt.0.99d0 .or. lssys%hphi_gp(:,:,5).lt.0.99d0) pq%biax_IMC = pq%biax
    ! where (lssys%hphi_gp(:,:,5).gt.0.99d0)                                     pq%biax_Sn  = pq%biax
    ! pq%biax_Sn  = pq%biax*lssys%hphi_gp_plot(:,:,6) + pq%biax*lssys%hphi_gp_plot(:,:,7) + pq%biax*lssys%hphi_gp_plot(:,:,8)
    pq%biax_Sn  = pq%biax*lssys%sn_hphi_gp_plot
    pq%biax_Sn2 = pq%biax*lssys%sn_hphi_gp_plot2

    ! Compute volume averaged biaxial stress in Sn layer
    call compute_biax_Sn_avg2(pq%biax_Sn_avg, mesh%newcoord, mesh%enod, mesh%nelm, pq%biax_Sn2, lssys)

    ! Compute plastic response
    call J2iso_Sn_getVal('plastic',pq%plastic)

    ! Compute hydrostatic pressure (only in Sn phase)
    pq%pgp = 0d0
    ! where (lssys%hphi_gp(:,:,4).gt.0.99d0) pq%pgp = pq%biax
    pq%pgp = pq%biax
    call integ2nod(pq%p_nod,pq%pgp,'qu4',mesh%enod)
    call extract(pq%p_ed, pq%p_nod, mesh%enod)

    return
  end subroutine compute_plot_quantities


  subroutine plot_vtk(pq, input_location, nnodgp, mesh, lssys, i_IMC)
    implicit none

    ! Intent inout
    type(plot_system), intent(inout)         :: pq
    type(mesh_system), intent(inout)         :: mesh
    integer, intent(inout)                   :: nnodgp

    ! Intent in
    type(ls_system), intent(in)              :: lssys
    character(len=:),allocatable, intent(in) :: input_location
    integer, intent(in)                      :: i_IMC

    ! Subroutine variables
    character(len=45)                        :: vtk_filename, fieldname
    integer                                  :: g, i
    

    ! --- Enter VTK folder -----
    call chdir(input_location)
    call chdir('VTK')

    ! --- Open main VTK ---
    if (i_IMC.lt.10) then
      write (vtk_filename, "(A,I1,A)") 'main_', i_IMC,'.vtk'
    elseif (i_IMC.lt.100) then
      write (vtk_filename, "(A,I2,A)") 'main_', i_IMC,'.vtk'
    elseif (i_IMC.lt.1000) then
      write (vtk_filename, "(A,I3,A)") 'main_', i_IMC,'.vtk'
    endif

    
    call VTKopen(trim(vtk_filename))

    ! Write mesh in current configuration to main.vtk
    call VTKPlotMesh('qu4', mesh%enod, mesh%newcoord)

    ! Write material data to main.vtk
    ! call VTKPlotNodeVal(real(lssys%material_nod, real(dp)),'material_at_nods')

    ! Write a_glob to main.vtk
    write(fieldname,'(A8,I1)') 'a_glob'
    call VTKPlotNodeVal(-minval(lssys%a,2),trim(fieldname))    

    ! Write phi data to main.vtk
    do g=1,lssys%ngrains
      write(fieldname,'(A8,I1)') 'phi_', g
      call VTKPlotNodeVal(lssys%phi(:,g),trim(fieldname))
    enddo

    ! Write hphi data to main.vtk
    do g=1,lssys%ngrains
      write(fieldname,'(A8,I1)') 'hphi_', g
      call VTKPlotNodeVal(lssys%hphi(:,g),trim(fieldname))
    enddo

    ! Write sn_hphi_plot data to main.vtk
    write(fieldname,'(A12,I1)') 'sn_hphi_plot'
    call VTKPlotNodeVal(lssys%sn_hphi_plot,trim(fieldname))

    ! Write all vm data to main.vtk
    call integ2nod(pq%nodplot,pq%vm,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'vm_avaraged_to_nodes')

    ! Write vm_Cu to main.vtk
    call integ2nod(pq%nodplot,pq%vm_Cu,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'vm_Cu_avaraged_to_nodes')

    ! Write vm_IMC to main.vtk
    call integ2nod(pq%nodplot,pq%vm_IMC,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'vm_IMC_avaraged_to_nodes')

    ! Write vm_Sn to main.vtk
    call integ2nod(pq%nodplot,pq%vm_Sn,'qu4',mesh%enod)    
    call VTKPlotNodeVal(pq%nodplot,'vm_Sn_avaraged_to_nodes')


    
    ! Write all biax data
    call integ2nod(pq%nodplot,pq%biax,'qu4',mesh%enod)
    pq%nodplot2 = pq%nodplot
    call VTKPlotNodeVal(pq%nodplot,'biax_avaraged_to_nodes')

    ! Write biax data for only Cu
    call integ2nod(pq%nodplot,pq%biax_Cu,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'biax_Cu_avaraged_to_nodes')

    ! Write biax data for only IMC
    call integ2nod(pq%nodplot,pq%biax_IMC,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'biax_IMC_avaraged_to_nodes')

    ! Write biax data for only Sn
    call integ2nod(pq%nodplot,pq%biax_Sn,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'biax_Sn_avaraged_to_nodes')

    ! Write hydrostatic pressure    
    call VTKPlotNodeVal(pq%p_nod,'p_avaraged_to_nodes')

    ! Write data for E
    call integ2nod(pq%nodplot,pq%E_save_gp,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'E_avaraged_to_nodes')

    ! Write plastic zone to main.vtk    
    call integ2nod(pq%nodplot,pq%plastic,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'plastic_avaraged_to_nodes')

    ! Write hardening parameter alpha to main.vtk
    call J2iso_Sn_getVal('alpha',pq%alpha)
    call integ2nod(pq%nodplot,pq%alpha,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'alpha_avaraged_to_nodes')

    ! Write initial IMC to main.vtk
    pq%alpha = 0d0
    where (lssys%old_IMC_gp_flag) pq%alpha = 1d0
    call integ2nod(pq%nodplot,pq%alpha,'qu4',mesh%enod)
    call VTKPlotNodeVal(pq%nodplot,'initial_IMC')

    ! ! Write element sets
    ! call VTKPlotCellVal(mesh%element_set,'element_set')

    ! Close main.vtk
    call VTKclose()




    ! --- Calculate gp material ---

    ! 1) Define enodgp = reshape(1:nnodgp,(4,nelm)). Note completely different from enod.
    pq%enodgp = reshape([(i,i=1,mesh%nnodgp)],[4,mesh%nelm])

    ! 2) Get coords of gauss points for all nods in gp order 1,2,4,3 (counterclockwise)
    call gpxtr(mesh%gpx, mesh%gpy, mesh%newcoord, mesh%enod)
    call getcoords(pq%coordgp, mesh%gpx([1,2,4,3],:), mesh%gpy([1,2,4,3],:), pq%enodgp)

    ! ! 4) Label all gauss points with material
    ! pq%gp_mat_label = 0d0
    ! where (grains%Cu_gps)      pq%gp_mat_label = 1d0
    ! where (grains%IMC_gps)     pq%gp_mat_label = 2d0
    ! where (grains%IMCSn_gps)   pq%gp_mat_label = 3d0
    ! where (grains%Sn_gps)      pq%gp_mat_label = 4d0
    ! pq%gp_mat_label_nodplot = reshape(pq%gp_mat_label([1,2,4,3],:),[nnodgp])

    ! ! Phi
    ! pq%gp_phi_label = 0d0
    ! pq%gp_phi_label_nodplot = reshape(phi_gp([1,2,4,3],:),[nnodgp])



    ! --- Open mesh_with_gpmat.vtk and store gp material data ---
    if (i_IMC.lt.10) then
      write (vtk_filename, "(A,I1,A)") 'mesh_with_gpmat_', i_IMC,'.vtk'
    elseif (i_IMC.lt.100) then
      write (vtk_filename, "(A,I2,A)") 'mesh_with_gpmat_', i_IMC,'.vtk'
    elseif (i_IMC.lt.1000) then
      write (vtk_filename, "(A,I3,A)") 'mesh_with_gpmat_', i_IMC,'.vtk'
    endif    

    call VTKopen(trim(vtk_filename))
    call VTKPlotMesh('qu4',pq%enodgp,pq%coordgp)


    ! Write material_gp to mesh_with_gpmat.vtk
    ! call integ2nod(pq%nodplot_gp,real(lssys%material_gp, real(dp)),'qu4',pq%enodgp)
    ! call VTKPlotNodeVal(pq%nodplot_gp,'material_at_gp')

    ! Write phi_gp data to mesh_with_gpmat.vtk
    do g=1,lssys%ngrains
      call integ2nod(pq%nodplot_gp,lssys%phi_gp(:,:,g),'qu4',pq%enodgp)
      write(fieldname,'(A7,I1)') 'phi_gp_', g
      call VTKPlotNodeVal(pq%nodplot_gp, trim(fieldname))
    enddo

    ! Write E_save_gp to mesh_with_gpmat.vtk
    call integ2nod(pq%nodplot_gp,pq%E_save_gp,'qu4',pq%enodgp)
    call VTKPlotNodeVal(pq%nodplot_gp,'E_save_gp')

    ! call VTKPlotNodeVal(pq%gp_mat_label_nodplot,'gp_mat_label')
    ! call VTKPlotNodeVal(pq%gp_phi_label_nodplot,'gp_phi_active')
    ! call VTKPlotNodeVal(reshape(pq%plastic([1,2,4,3],:),[nnodgp]),'gp_plastic')
    call VTKclose()
    
    return
  end subroutine plot_vtk




  subroutine write_stress_iter_to_matlab(pq, mesh, lssys, input_location, itot)
    ! Routine for writing to matlab in Newton iteration
    implicit none

    ! Intent inout
    type(plot_system), intent(inout)          :: pq

    ! Intent in
    type(mesh_system), intent(in)             :: mesh
    type(ls_system), intent(in)               :: lssys
    character(len=:), allocatable, intent(in) :: input_location
    integer, intent(in)                       :: itot

    ! Subroutine variables
    character(len=40)                         :: matfilename

    ! Change to input location
    call chdir(input_location)
    call chdir('mat_files')

    ! Filename
    write (matfilename, "(A,I0,A)") 'stress_', itot,'.mat'
    print*, 'mat_iter_file: ', matfilename

    ! --- Write stress to file --- 

    ! Full matrices
    call matWrt2f(trim(matfilename), pq%vm           , 'vm'              ,'w')

    ! Volume average biaxial stress
    call matWrt2f(trim(matfilename), [pq%biax_Sn_avg], 'biax_Sn_avg'     ,'u')

    ! Vm    
    call integ2nod(pq%nodplot,pq%vm,'qu4',mesh%enod)
    call matWrt2f(trim(matfilename), pq%nodplot      , 'vm_nodplot'      ,'u')

    ! Vm Cu
    call integ2nod(pq%nodplot,pq%vm_Cu,'qu4',mesh%enod)
    call matWrt2f(trim(matfilename), pq%nodplot      , 'vm_Cu_nodplot'      ,'u')

    ! Vm IMC
    call integ2nod(pq%nodplot,pq%vm_IMC,'qu4',mesh%enod)
    call matWrt2f(trim(matfilename), pq%nodplot      , 'vm_IMC_nodplot'      ,'u')

    ! Vm Sn
    call integ2nod(pq%nodplot,pq%vm_Sn,'qu4',mesh%enod)
    call matWrt2f(trim(matfilename), pq%nodplot      , 'vm_Sn_nodplot'      ,'u')

    ! Biax
    call integ2nod(pq%nodplot,pq%biax,'qu4',mesh%enod)    
    call matWrt2f(trim(matfilename), pq%nodplot      , 'biax_nodplot'    ,'u')

    ! Biax Cu
    call integ2nod(pq%nodplot,pq%biax_Cu,'qu4',mesh%enod)
    call matWrt2f(trim(matfilename), pq%nodplot      , 'biax_Cu_nodplot'      ,'u')

    ! Biax IMC
    call integ2nod(pq%nodplot,pq%biax_IMC,'qu4',mesh%enod)    
    call matWrt2f(trim(matfilename), pq%nodplot      , 'biax_IMC_nodplot'      ,'u')

    ! Biax Sn
    call integ2nod(pq%nodplot,pq%biax_Sn,'qu4',mesh%enod)    
    call matWrt2f(trim(matfilename), pq%nodplot      , 'biax_Sn_nodplot'      ,'u')

    ! Biax Sn 2
    call integ2nod(pq%nodplot,pq%biax_Sn2,'qu4',mesh%enod)    
    call matWrt2f(trim(matfilename), pq%nodplot      , 'biax_Sn2_nodplot'      ,'u')

    ! Hydrostatic pressure
    call matWrt2f(trim(matfilename), pq%p_nod        , 'p_nod'      ,'u')

    call matWrt2f(trim(matfilename), mesh%newex      , 'newex'        , 'u')
    call matWrt2f(trim(matfilename), mesh%newey      , 'newey'        , 'u')
    call matWrt2f(trim(matfilename), mesh%newcoord   , 'newcoord'     , 'u')    


    ! call integ2nod(pq%nodplot,pq%vm_Cu,'qu4',mesh%enod)
    ! call matWrt2f(trim(matfilename), pq%nodplot          , 'vm_Cu_nodplot'   ,'u')
    ! call integ2nod(pq%nodplot,pq%vm_IMC,'qu4',mesh%enod)
    ! call matWrt2f(trim(matfilename), pq%nodplot          , 'vm_IMC_nodplot'  ,'u')
    ! call integ2nod(pq%nodplot,pq%vm_Sn,'qu4',mesh%enod)
    ! call matWrt2f(trim(matfilename), pq%nodplot          , 'vm_Sn_nodplot'   ,'u')    
    ! call integ2nod(pq%nodplot,pq%biax,'qu4',mesh%enod)
    ! call matWrt2f(trim(matfilename), pq%nodplot          , 'biax_nodplot'    ,'u')
    ! call integ2nod(pq%nodplot,pq%biax_Cu,'qu4',mesh%enod)
    ! call matWrt2f(trim(matfilename), pq%nodplot          , 'biax_Cu_nodplot' ,'u')
    ! call integ2nod(pq%nodplot,pq%biax_IMC,'qu4',mesh%enod)
    ! call matWrt2f(trim(matfilename), pq%nodplot          , 'biax_IMC_nodplot','u')
    ! call integ2nod(pq%nodplot,pq%biax_Sn,'qu4',mesh%enod)
    ! call matWrt2f(trim(matfilename), pq%nodplot          , 'biax_Sn_nodplot' ,'u')
    

    return
  end subroutine write_stress_iter_to_matlab


  subroutine write_level_set_init_to_matlab(mesh, IMC_steps, ngrains, input_location)
    ! Routine for writing initial level set data to matlab
    implicit none

    ! Intent in    
    type(mesh_system), intent(in)             :: mesh
    integer, intent(in)                       :: IMC_steps, ngrains
    character(len=:), allocatable, intent(in) :: input_location

    ! Subroutine variables
    character(len=40)                         :: matfilename

    call chdir(input_location)
    call chdir('mat_files')

    ! Filename
    write (matfilename, "(A,I0,A)") 'level_setinit.mat'
    write(*,'(A16,A)') 'Saving matfile: ', matfilename    

    ! Write to file
    call matWrt2f(trim(matfilename), mesh%coord  , 'coord'    , 'w')
    call matWrt2f(trim(matfilename), mesh%enod   , 'enod'     , 'u')
    call matWrt2f(trim(matfilename), [mesh%nodel], 'nodel'    , 'u')
    call matWrt2f(trim(matfilename), mesh%ex     , 'ex'       , 'u')
    call matWrt2f(trim(matfilename), mesh%ey     , 'ey'       , 'u')
    call matWrt2f(trim(matfilename), mesh%gpx    , 'gpx'      , 'u')
    call matWrt2f(trim(matfilename), mesh%gpy    , 'gpy'      , 'u')
    call matWrt2f(trim(matfilename), [IMC_steps] , 'IMC_steps', 'u')
    call matWrt2f(trim(matfilename), [ngrains]   , 'ngrains'  , 'u')

    return
  end subroutine write_level_set_init_to_matlab


  subroutine write_level_set_iter_to_matlab(lssys, mesh, i_IMC, input_location)
    ! Routine for writing level set to matlab
    implicit none

    ! Intent in
    type(ls_system), intent(in)               :: lssys
    type(mesh_system), intent(in)             :: mesh
    integer, intent(in)                       :: i_IMC
    character(len=:), allocatable, intent(in) :: input_location    

    ! Subroutine variables
    character(len=40)                         :: matfilename


    call chdir(input_location)
    call chdir('mat_files')

    ! Filename
    if (i_IMC.lt.10) then
      write (matfilename, "(A10,I1,A4)") 'level_set_', i_IMC,'.mat'
    elseif (i_IMC.lt.100) then
      write (matfilename, "(A10,I2,A4)") 'level_set_', i_IMC,'.mat'
    elseif (i_IMC.lt.1000) then
      write (matfilename, "(A10,I3,A4)") 'level_set_', i_IMC,'.mat'
    endif
    write(*,'(A16,A)') 'Saving matfile: ', matfilename

    ! Write to file    
    call matWrt2f(trim(matfilename), lssys%a                                , 'a'            , 'w')
    call matWrt2f(trim(matfilename), lssys%line_ex                          , 'line_ex'      , 'u')
    call matWrt2f(trim(matfilename), lssys%line_ey                          , 'line_ey'      , 'u')
    call matWrt2f(trim(matfilename), lssys%line_seg                         , 'line_seg'     , 'u')
    call matWrt2f(trim(matfilename), lssys%line_coord                       , 'line_coord'   , 'u')
    call matWrt2f(trim(matfilename), lssys%line_coordN                      , 'line_coordN'  , 'u')
    call matWrt2f(trim(matfilename), lssys%vp(1,:,:)                        , 'vpx'          , 'u')
    call matWrt2f(trim(matfilename), lssys%vp(2,:,:)                        , 'vpy'          , 'u')
    call matWrt2f(trim(matfilename), lssys%sep_lines                        , 'sep_lines'    , 'u')
    call matWrt2f(trim(matfilename), [lssys%time]                           , 'time'         , 'u')
    call matWrt2f(trim(matfilename), lssys%material                         , 'material'     , 'u')
    call matWrt2f(trim(matfilename), [lssys%IMC_area_init]                  , 'IMC_area_init', 'u')
    call matWrt2f(trim(matfilename), [lssys%IMC_area]                       , 'IMC_area'     , 'u')
    call matWrt2f(trim(matfilename), [lssys%time]                           , 'time'         , 'u')
    call matWrt2f(trim(matfilename), [lssys%time]                           , 'time'         , 'u')
    call matWrt2f(trim(matfilename), [mesh%elmsize_x]                       , 'elsize_x'     , 'u')
    call matWrt2f(trim(matfilename), [mesh%elmsize_y]                       , 'elsize_y'     , 'u')    
    call matWrt2f(trim(matfilename), mesh%newex                             , 'newex'        , 'u')
    call matWrt2f(trim(matfilename), mesh%newey                             , 'newey'        , 'u')
    call matWrt2f(trim(matfilename), mesh%newcoord                          , 'newcoord'     , 'u')
    call matWrt2f(trim(matfilename), mesh%enod                              , 'enod'         , 'u')
    call matWrt2f(trim(matfilename), lssys%tp_points([1:lssys%ntp_points],:), 'tppoints'     , 'u')
    call matWrt2f(trim(matfilename), lssys%hphi                             , 'hphi'     , 'u')
    

    return
  end subroutine write_level_set_iter_to_matlab

  subroutine write_level_set_unsorted_lines_iter_to_matlab(line_ex, line_ey, i_IMC, input_location)
    ! Routine for writing level set to matlab
    implicit none

    character(len=:), allocatable, intent(in) :: input_location

    integer, intent(in)                       :: i_IMC

    real(dp), intent(in)                      :: line_ex(:,:), line_ey(:,:)
    character(len=40)                         :: matfilename

    call chdir(input_location)
    call chdir('mat_files')

    ! Filename
    write (matfilename, "(A,I0,A)") 'level_setunsorted_lines_', i_IMC,'.mat'
    write(*,'(A16,A)') 'Saving matfile: ', matfilename

    ! Write to file
    call matWrt2f(trim(matfilename), line_ex       , 'line_ex'       , 'w')
    call matWrt2f(trim(matfilename), line_ey       , 'line_ey'       , 'u')

    return
  end subroutine write_level_set_unsorted_lines_iter_to_matlab


  subroutine write_diffusion_iter_to_matlab(a, r, ed, jint, j_flux, i_IMC, input_location, g)
    ! Routine for writing level set to matlab
    implicit none

    ! Intent in
    character(len=:), allocatable, intent(in) :: input_location
    integer, intent(in)                       :: i_IMC, g
    real(dp), intent(in)                      :: a(:), r(:), ed(:,:), jint(:), j_flux(:,:)
    character(len=25)                         :: matfilename

    call chdir(input_location)
    call chdir('mat_files')

    ! Filename
    if (i_IMC.lt.10) then
      write (matfilename, "(A10,I1,A1,I1,A4)") 'diffusion_', i_IMC, '_', g, '.mat'
    elseif (i_IMC.lt.100) then
      write (matfilename, "(A10,I2,A1,I1,A4)") 'diffusion_', i_IMC, '_', g, '.mat'
    elseif (i_IMC.lt.1000) then
      write (matfilename, "(A10,I3,A1,I1,A4)") 'diffusion_', i_IMC, '_', g, '.mat'
    endif
    
    ! Write to file
    call matWrt2f(trim(matfilename), a     , 'a'     , 'w')
    call matWrt2f(trim(matfilename), r     , 'r'     , 'u')
    call matWrt2f(trim(matfilename), ed    , 'ed'    , 'u')
    call matWrt2f(trim(matfilename), jint  , 'jint'  , 'u')
    call matWrt2f(trim(matfilename), j_flux, 'j_flux', 'u')

    write(*,'(A16,A)') 'Saving matfile: ', matfilename

    return
  end subroutine write_diffusion_iter_to_matlab

  subroutine write_diffusion_glob_iter_to_matlab(input_location, i_IMC, ls)
    ! Routine for writing level set to matlab
    implicit none

    ! Intent in
    character(len=:), allocatable, intent(in) :: input_location
    integer, intent(in)                       :: i_IMC
    real(dp), intent(in)                      :: ls(:,:)
    character(len=25)                         :: matfilename

    call chdir(input_location)
    call chdir('mat_files')

    ! Filename
    if (i_IMC.lt.10) then
      write (matfilename, "(A12,I1,A4)") 'triangle_ls_', i_IMC, '.mat'
    elseif (i_IMC.lt.100) then
      write (matfilename, "(A12,I2,A4)") 'triangle_ls_', i_IMC, '.mat'
    elseif (i_IMC.lt.1000) then
      write (matfilename, "(A12,I3,A4)") 'triangle_ls_', i_IMC, '.mat'
    endif

    write(*,'(A16,A)') 'Saving matfile: ', matfilename

    ! Write to file
    call matWrt2f(trim(matfilename), ls     , 'ls'     , 'w')

    return
  end subroutine write_diffusion_glob_iter_to_matlab

  ! subroutine append_diffusion_sn_vel_iter_to_matlab(vnod, i_IMC, input_location)
  !   ! Routine for writing level set to matlab
  !   implicit none

  !   ! Intent in
  !   character(len=:), allocatable, intent(in) :: input_location
  
  !   integer, intent(in)                       :: i_IMC
  !   real(dp), intent(in)                      :: vnod(:)
  !   character(len=40)                         :: matfilename

  !   call chdir(input_location)
  !   call chdir('mat_files')

  !   ! Filename
  !   write (matfilename, "(A,I0,A)") 'diffusion_sn', i_IMC,'.mat'
  !   print*, 'Saving matfile: ', matfilename

  !   ! Append to file
  !   call matWrt2f(trim(matfilename), vnod, 'vnod', 'u')

  !   return
  ! end subroutine append_diffusion_sn_vel_iter_to_matlab


  subroutine write_times_to_matlab(exec_time, init_time, stiffness_time, assem_stiffness_time, spaputval_time, solveq_time, &
                                   force_time, assem_force_time, mat_time, plotvtk_time, input_location)
    ! Routine for writing times to matlab
    implicit none

    ! Intent in
    real(dp), intent(in)    ::  exec_time, init_time, stiffness_time, assem_stiffness_time, spaputval_time, solveq_time
    real(dp), intent(in)    ::  force_time, assem_force_time, mat_time, plotvtk_time
    character(len=:), allocatable, intent(in) :: input_location

    ! Subroutine variables
    character(len=40)       :: matfilename

    call chdir(input_location)
    call chdir('mat_files')

    ! Filename
    write (matfilename, "(A,I0,A)") 'times.mat'
    print*, 'Saving matfile: ', matfilename

    ! Write to file
    call matWrt2f(trim(matfilename), [exec_time]           , 'exec_time'           , 'w')
    call matWrt2f(trim(matfilename), [init_time]           , 'init_time'           , 'u')
    call matWrt2f(trim(matfilename), [stiffness_time]      , 'stiffness_time'      , 'u')
    call matWrt2f(trim(matfilename), [assem_stiffness_time], 'assem_stiffness_time', 'u')
    call matWrt2f(trim(matfilename), [spaputval_time]      , 'spaputval_time'      , 'u')
    call matWrt2f(trim(matfilename), [solveq_time]         , 'solveq_time'         , 'u')
    call matWrt2f(trim(matfilename), [force_time]          , 'force_time'          , 'u')
    call matWrt2f(trim(matfilename), [assem_force_time]    , 'assem_force_time'    , 'u')
    call matWrt2f(trim(matfilename), [mat_time]            , 'mat_time'            , 'u')
    call matWrt2f(trim(matfilename), [plotvtk_time]        , 'plotvtk_time'        , 'u')

    return
  end subroutine write_times_to_matlab
  

end module output_data
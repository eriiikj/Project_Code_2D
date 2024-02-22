module diffusion
    ! Last modified
    ! E. Jacobsson 2023-04-11
      
    ! mf_datatypes
    use mf_datatypes

    ! Somelib
    use matrix_util

    ! whiskerlib
    use diffusion_types
    use output_data
    

    implicit none

    ! This is an attempt to predefine some vectors for the 4-node element, to obtain speed
    double precision, parameter  :: G1=0.577350269189626D0
    double precision, parameter  :: xsi4(4)=(/-1D0,  1D0, -1D0, 1D0/)*G1
    double precision, parameter  :: eta4(4)=(/-1D0, -1D0,  1D0, 1D0/)*G1
    private xsi4, eta4, G1

    ! 4 Gauss points and 4 nodes
    double precision, parameter :: NR4_col11(4) = (1D0-XSI4)*(1D0-ETA4)/4D0
    double precision, parameter :: NR4_col22(4) = (1D0+XSI4)*(1D0-ETA4)/4D0
    double precision, parameter :: NR4_col33(4) = (1D0+XSI4)*(1D0+ETA4)/4D0
    double precision, parameter :: NR4_col44(4) = (1D0-XSI4)*(1D0+ETA4)/4D0
    double precision, parameter :: NR4(4,4)    = reshape([NR4_col11, NR4_col22, NR4_col33, NR4_col44], [4,4])
    private NR4
    private NR4_col11, NR4_col22, NR4_col33, NR4_col44

    ! derivate of shape functions with respect to xsi
    double precision, parameter :: DNR4_182_col1(4) = -(1D0-ETA4)/4D0
    double precision, parameter :: DNR4_182_col2(4) =  (1D0-ETA4)/4D0
    double precision, parameter :: DNR4_182_col3(4) =  (1D0+ETA4)/4D0
    double precision, parameter :: DNR4_182_col4(4) = -(1D0+ETA4)/4D0

    ! derivate of shape functions with respect to eta
    double precision, parameter :: DNR4_282_col1(4) = -(1D0-XSI4)/4D0
    double precision, parameter :: DNR4_282_col2(4) = -(1D0+XSI4)/4D0
    double precision, parameter :: DNR4_282_col3(4) =  (1D0+XSI4)/4D0
    double precision, parameter :: DNR4_282_col4(4) =  (1D0-XSI4)/4D0

    ! Collect in DNR4
    double precision, parameter :: DNR4(8,4)=reshape([DNR4_182_col1(1), DNR4_282_col1(1), DNR4_182_col1(2), DNR4_282_col1(2), &
                                                    DNR4_182_col1(3), DNR4_282_col1(3), DNR4_182_col1(4), DNR4_282_col1(4), &
                                                    DNR4_182_col2(1), DNR4_282_col2(1), DNR4_182_col2(2), DNR4_282_col2(2), &
                                                    DNR4_182_col2(3), DNR4_282_col2(3), DNR4_182_col2(4), DNR4_282_col2(4), &
                                                    DNR4_182_col3(1), DNR4_282_col3(1), DNR4_182_col3(2), DNR4_282_col3(2), &
                                                    DNR4_182_col3(3), DNR4_282_col3(3), DNR4_182_col3(4), DNR4_282_col3(4), &
                                                    DNR4_182_col4(1), DNR4_282_col4(1), DNR4_182_col4(2), DNR4_282_col4(2), &
                                                    DNR4_182_col4(3), DNR4_282_col4(3), DNR4_182_col4(4), DNR4_282_col4(4)], [8,4])
    private DNR4
    private DNR4_182_col1, DNR4_182_col2, DNR4_182_col3, DNR4_182_col4
    private DNR4_282_col1, DNR4_282_col2, DNR4_282_col3, DNR4_282_col4


contains


subroutine init_diffusion_system(diffsys, ngrains)
  ! Subroutine for initiating diffusion system
  implicit none

  ! Intent inout
  type(diffusion_system), intent(inout) :: diffsys 

  ! Intent in
  integer, intent(in) :: ngrains
 
  ! Subroutine variables
  integer ::ierr

  ! Allocate grain mesh system
  allocate(diffsys%grain_meshes(ngrains),stat=ierr)  

  return
end subroutine init_diffusion_system


subroutine solve_diffusion_problem_global(diffsys,i_IMC,lssys,mesh,pq,input_location,omp_run)
  ! Solve diffusion problem and compute new velocity field

  ! Intent inout
  type(diffusion_system), intent(inout)      :: diffsys 
  type(ls_system), intent(inout)             :: lssys 

  ! Intent in
  type(mesh_system), intent(in)              :: mesh
  character(len=:), allocatable, intent(in)  :: input_location
  logical, intent(in)                        :: omp_run
  type(plot_system), intent(inout)           :: pq
  integer, intent(in)                        :: i_IMC

  ! Subroutine variables
  integer :: g
  
  ! Obtain global diffusion mesh
  call generate_global_diffusion_mesh(input_location, diffsys, mesh, lssys, i_IMC)

  ! Solve diffusion problem for all grains seperately
  do g=1,lssys%ngrains    
    call diffusion_grain(diffsys,i_IMC, mesh, diffsys%grain_meshes(g), g, input_location, lssys, omp_run, pq)
  enddo 

  ! ! Compute IMC area
  ! call compute_IMC_area(lssys,diffsys)
  ! if (i_IMC.eq.1) then
  !   lssys%IMC_area_init = lssys%IMC_area
  ! endif  

  return
end subroutine solve_diffusion_problem_global



subroutine generate_global_diffusion_mesh(input_location, diffsys, mesh, lssys, i_IMC)
  ! --- Routine for generating and partiotioning a global diffusion mesh ---
  implicit none

  ! Intent inout
  type(diffusion_system), intent(inout) :: diffsys 

  ! Intent in
  type(mesh_system), intent(in)         :: mesh
  type(ls_system), intent(in)           :: lssys
  integer, intent(in)                   :: i_IMC
  character(len=:), allocatable, intent(in)  :: input_location

  ! Subroutine variables
  character(len=40)                     :: command_line  
  character(len=80)                     :: diffusion_mesh_script_location 
  character(len=:),allocatable          :: diffusion_mesh_location  
  integer                               :: g, ierr, i, g_cols(2), nelm, ii(1), inod
  integer, allocatable :: elmidx(:), nodidx(:), nodidx2(:), uniqueArr(:)
  logical :: ch

  diffusion_mesh_script_location = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output/phase_mesh_script'  
  diffusion_mesh_location        = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output/single_study/phase_meshes'

  ! Generete mesh by running python script
  call chdir(diffusion_mesh_script_location)
  write(command_line,'(A28)'), 'python triangle_mesh_Main.py'
  call execute_command_line(trim(command_line))

  ! Read global diffusion mesh
  if (allocated(diffsys%ls)) then
    deallocate(diffsys%ls)
  endif
  if (allocated(diffsys%coord)) then
    deallocate(diffsys%coord)
  endif
  if (allocated(diffsys%enod)) then
    deallocate(diffsys%enod)
  endif  

  ! Generate new diffusion mesh
  call read_json_diffusion_mesh(diffusion_mesh_location, i_IMC, diffsys%coord, diffsys%enod)  

  diffsys%nnod  = size(diffsys%coord,2)
  diffsys%nelm  = size(diffsys%enod,2)
  diffsys%nodel = size(diffsys%enod,1)

  ! Allocated ls in new mesh
  if (allocated(diffsys%ls)) then
    deallocate(diffsys%ls)
  endif
  if (allocated(diffsys%ls_ed)) then
    deallocate(diffsys%ls_ed)
  endif
  if (allocated(diffsys%closest_line)) then
    deallocate(diffsys%closest_line)
  endif
  allocate(diffsys%ls(diffsys%nnod,lssys%ngrains),stat=ierr)
  allocate(diffsys%ls_ed(diffsys%nodel,diffsys%nelm,lssys%ngrains),stat=ierr)
  allocate(diffsys%closest_line(diffsys%nnod,lssys%ngrains),stat=ierr)
  diffsys%ls    = 0d0
  diffsys%ls_ed = 0d0

  ! Interpolate old level set values to new mesh
  do g=1,lssys%ngrains
    call interpolate_scalar_mesh1_to_mesh2(mesh%enod, mesh%coord, lssys%ed(:,:,g), diffsys%enod, diffsys%coord, diffsys%ls(:,g))
  enddo

  if (allocated(elmidx)) then
    deallocate(elmidx)
  endif
  allocate(elmidx(diffsys%nelm))
  do i=1,diffsys%nelm
    elmidx(i)=i
  enddo

  if (allocated(nodidx)) then
    deallocate(nodidx)
  endif
  allocate(nodidx(diffsys%nnod))
  do i=1,diffsys%nnod
    nodidx(i)=i
  enddo

  ! --- Reinitialize triangle mesh ---
  do g = 1,lssys%ngrains
    g_cols = [2*(g-1) + 1, 2*(g-1) + 2]
    call reinit_level_set_spatial2(diffsys%ls(:,g),lssys%line_ex(:,g_cols(1):g_cols(2)),&
    lssys%line_ey(:,g_cols(1):g_cols(2)), diffsys%closest_line(:,g),lssys%line_seg(g),diffsys%coord)
    where(abs(diffsys%ls(:,g)).lt.1d-8) diffsys%ls(:,g) = 0d0

    ! Post-processing sign of reinitialized mesh as this is not always true from the last state
    print *, 'changing signs...'
    ch = .true.
    do while (ch)
      ch = .false.  
      nodidx2 = pack(nodidx,diffsys%ls(:,g).lt.0d0)    
      do i=1,size(nodidx2)
        inod = nodidx2(i)
        ! Find elements containing node
        nelm = count(count(diffsys%enod.eq.inod,1).gt.0)
        call FindUniqueIntegers(reshape(diffsys%enod(:,pack(elmidx, count(diffsys%enod.eq.inod,1).gt.0)),[nelm*diffsys%nodel]), &
        uniqueArr)
        ii = inod
        call setdiff_int_1D(uniqueArr, ii, uniqueArr)
        ! if (abs(diffsys%coord(1,i)-0.334833d-3).lt.1d-7 .and. abs(diffsys%coord(2,i)-0.58539d-3).lt.1d-7) then
        !   print *, 'coord', diffsys%coord(:,i)
        !   print *, 'as', diffsys%ls(uniqueArr,g)
        ! endif
        if (count((diffsys%ls(uniqueArr,g)).lt.0d0).ge.count((diffsys%ls(uniqueArr,g)).gt.0d0) .and.diffsys%ls(inod,g).gt.0d0) then
          diffsys%ls(inod,g) = -1d0*abs(diffsys%ls(inod,g))
          ch = .true.
        endif
      enddo
    enddo

    ! Update ed
    call extract(diffsys%ls_ed(:,:,g),diffsys%ls(:,g),diffsys%enod,1)
  enddo

  ! Write interpolated ls data to file
  call write_diffusion_glob_iter_to_matlab(input_location, i_IMC, diffsys%ls)

  ! Generete partitioned grain meshes by running python script
  call chdir(diffusion_mesh_script_location)
  write(command_line,'(A33)'), 'python extract_grain_data_Main.py'
  call execute_command_line(trim(command_line))

  return
end subroutine generate_global_diffusion_mesh




subroutine diffusion_grain(diffsys, i_IMC, global_mesh, grain_mesh, g, input_location, lssys, omp_run, pq)
  ! --- Routine for solving diffusion equation in grain ---
    implicit none
    
    ! Intent inout
    type(diffusion_system), intent(inout)      :: diffsys 
    type(grain_mesh_system),intent(inout)      :: grain_mesh    
    
    ! Intent in
    type(mesh_system), intent(in)              :: global_mesh
    type(ls_system), intent(in)                :: lssys 
    integer, intent(in)                        :: i_IMC, g
    character(len=:), allocatable, intent(in)  :: input_location
    logical, intent(in)                        :: omp_run
    type(plot_system), intent(inout)           :: pq

    ! Subroutine variables
    integer :: g_cols(2)

    ! Locations
    character(len=255)                         :: main_location               ! Location of this Fortran program
    character(len=255)                         :: phase_mesh_script_location  ! Location of phase mesh script
    character(len=:),allocatable               :: phase_mesh_location         ! Location of phase mesh json files
    character(len=40)                          :: command_line         

    ! Times
    real(dp)                                   :: t1, t2, tmp(2) = 1d0
    

    write(*,*) 
    write(*,'(A30,I1,A4)') '--- Generating mesh for grain ', g, ' ---'

    ! Reset mesh quantities
    print *, 'Before mesh reset'
    call mesh_reset(grain_mesh)
    print *, 'After mesh reset'

    
    ! Location of phase mesh script
    call getcwd(main_location)
    phase_mesh_script_location = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output/phase_mesh_script'
    phase_mesh_location        = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output/single_study/phase_meshes'


    ! ! Generete mesh by running python script
    ! call chdir(phase_mesh_script_location)
    ! call clock_time(t1, omp_run)
    ! write(command_line,'(A37,I1)'), 'python grain_mesh_triangular_Main.py ', g
    ! call execute_command_line(trim(command_line))
    ! call clock_time(t2, omp_run)
    ! diffsys%generate_mesh_time = diffsys%generate_mesh_time + (t2 - t1)    
    ! call chdir(main_location)

    ! Read mesh from json file
    call clock_time(t1, omp_run)    
    call read_json_phase_mesh2(phase_mesh_location, i_IMC, grain_mesh%coord, grain_mesh%enod, grain_mesh%bcnod, grain_mesh%bcval, &
    grain_mesh%bcval_idx, g)
    call clock_time(t2, omp_run)
    diffsys%read_mesh_time = diffsys%read_mesh_time + (t2 - t1)

    ! Init mesh system of grain
    call init_grain_mesh_system(grain_mesh, lssys, diffsys, g)

    ! Print mesh info
    call print_diffusion_mesh_info(grain_mesh)

    ! ! If Sn phase, interpolate hydrostatic pressure field from global mesh (stored in pq) to grain_mesh
    ! if (lssys%material(g).eq.3) then
    !   call interpolate_glob_to_grain(global_mesh,grain_mesh,pq%p_ed)
    ! endif

    ! Solve diffusion equation
    print *, 'Before solve diffusion problem'
    call clock_time(t1, omp_run)
    call solve_diffusion_problem(grain_mesh, global_mesh, lssys%a(:,g), lssys%material(g))
    call clock_time(t2, omp_run)
    diffsys%solve_time = diffsys%solve_time + (t2-t1)
    print *, 'After solve diffusion problem'


    ! Obtain jint of grain
    ! g_cols = [2*(g-1) + 1, 2*(g-1) + 2]
    ! print *, 'Before compute jint'
    ! call compute_jint2(grain_mesh, I_IMC)
    ! print *, 'After compute jint'
  
    ! Save diffusion results    
    call write_diffusion_iter_to_matlab(grain_mesh%a, grain_mesh%r,grain_mesh%ed, grain_mesh%jint, grain_mesh%j_flux, &
    grain_mesh%p, i_IMC, input_location, g)

    print *, 'Leaving diffusion_grain'

    return
end subroutine diffusion_grain


subroutine mesh_reset(mesh)
  implicit none

  ! Intent inout
  type(grain_mesh_system), intent(inout) :: mesh 


  ! Deallocate variables
  if (allocated(mesh%coord)) then
    deallocate(mesh%coord)
    deallocate(mesh%enod)
    deallocate(mesh%a)
    deallocate(mesh%r)
    deallocate(mesh%bcnod)
  endif

  if (allocated(mesh%bcval)) then
    deallocate(mesh%bcval)
  endif

  if (allocated(mesh%bcval_idx)) then
    deallocate(mesh%bcval_idx)
  endif

  if (allocated(mesh%jint)) then
    deallocate(mesh%jint)
  endif


  ! Reset variables
  ! mesh%nelm   = 0
  ! mesh%nnod   = 0
  ! mesh%coord  = 0d0
  ! mesh%enod   = 0
  ! mesh%a      = 0d0
  ! mesh%r      = 0d0
  ! mesh%ed     = 0d0
  ! mesh%j_flux = 0d0
  ! mesh%jint   = 0d0
  ! mesh%bcnod  = 0
  ! mesh%bcval  = 0d0

  return
end subroutine mesh_reset


subroutine init_grain_mesh_system(mesh, lssys, diffsys, g)
  implicit none

  ! Intent inout
  type(diffusion_system), intent(inout)  :: diffsys
  type(grain_mesh_system), intent(inout) :: mesh

  ! Intent in
  type(ls_system), intent(in)            :: lssys  
  integer, intent(in)                    :: g

  ! Subroutine variables    
  integer :: ierr  

  ! Number of dofs in each nod and nodel
  mesh%nelm   = size(mesh%enod,2)
  mesh%nnod   = size(mesh%coord,2)
  mesh%dofnod = 1
  mesh%nodel  = 3

  ! Allocate variables in grain_mesh. OBS allocation made larger in order to avoid reallocation
  allocate(mesh%a(mesh%nnod),stat=ierr)
  allocate(mesh%r(mesh%nnod),stat=ierr)
  allocate(mesh%ed(mesh%nodel,mesh%nelm), stat=ierr)
  allocate(mesh%j_flux(2,mesh%nelm),stat=ierr)  

  ! Bulk diffusion mobility
  mesh%M = diffsys%D_diff_material(lssys%material(g))  / &
  (diffsys%molar_volumes(lssys%material(g))*diffsys%thermo_parameterA(lssys%material(g)))

  ! Init variables
  mesh%a      = 0d0
  mesh%r      = 0d0
  mesh%ed     = 0d0
  mesh%j_flux = 0d0
  
  return
end subroutine init_grain_mesh_system


subroutine solve_diffusion_problem(grain_mesh, global_mesh, lssys_ag, material)
  implicit none

  ! Intent inout
  type(grain_mesh_system), intent(inout) :: grain_mesh

  ! Intent in  
  type(mesh_system), intent(in)          :: global_mesh
  real(dp), intent(in)                   :: lssys_ag(:)
  integer, intent(in)                    :: material

  ! Subroutine variables
  real(dp)                               :: M(2,2)
  real(dp)                               :: ke3(3,3), scale_factor=1d30
  type(sparse)                           :: K
  integer                                :: ie, ierr
  real(dp), allocatable                  :: Ksave(:)

  real(dp)                               :: mean_distance, Mt(2,2), Mgb(2,2), alpha=9d3, mean_x, mean_y, model_edge_x_distance
  real(dp)                               :: Mgb_fac = 1d0
  integer :: minIdx(1)  

  ! Mobility [mol^2/(mm*s*J)]
  M(1,:) =  (/grain_mesh%M,0d0/)
  M(2,:) =  (/0d0,grain_mesh%M/)

  ! Scale M such that equation system is not below number precision
  M   = M*scale_factor  
  Mgb = Mgb_fac*M
  Mt  = M

  ! Stiffness K
  call spaTopDef(K,grain_mesh%enod,grain_mesh%dofnod)
  K%a = 0d0
  do ie=1,grain_mesh%nelm
    
    ! Change diffusion coefficient according to grain boundary diffusion
    ! Find point in global mesh closest to (mean_x,mean_y)
    ! Mean x and y of element
    ! mean_x = sum(grain_mesh%coord(1,grain_mesh%enod(:,ie)))/size(grain_mesh%enod(:,ie))
    ! mean_y = sum(grain_mesh%coord(2,grain_mesh%enod(:,ie)))/size(grain_mesh%enod(:,ie))
    ! minIdx = minloc((global_mesh%coord(1,:) - mean_x)**2 + (global_mesh%coord(2,:) - mean_y)**2)
    ! mean_distance         = abs(lssys_ag(minIdx(1)))
    ! model_edge_x_distance = minval([abs(global_mesh%model_width-mean_x), abs(0d0 - mean_x)])
    ! mean_distance         = minval([mean_distance,model_edge_x_distance]) 
    ! Mt = M + (Mgb - M)*exp(-alpha*model_edge_x_distance)
    ! if (material.eq.3) then
    !   mean_x                = sum(grain_mesh%coord(1,grain_mesh%enod(:,ie)))/size(grain_mesh%enod(:,ie))
    !   model_edge_x_distance = minval([abs(global_mesh%model_width-mean_x), abs(0d0 - mean_x)])
    !   Mt                    = M + (Mgb - M)*exp(-alpha*model_edge_x_distance)
    ! endif    
    
    call diffusion2D3kc(ke3, grain_mesh%coord(:,grain_mesh%enod(:,ie)),Mt)
    call assem(K,Ke3,grain_mesh%enod(:,ie),grain_mesh%dofnod)
  enddo

  ! Save K%a in Ksave
  allocate(Ksave(size(K%a)), stat=ierr)
  Ksave = K%a

  ! Solve Kc = f, OBS this alters K
  grain_mesh%a = 0d0 ! fb
  call solveq(K, grain_mesh%a, grain_mesh%bcnod, grain_mesh%bcval*scale_factor, grain_mesh%dofnod)

  ! Residual r = Kc-f = Kc (f=0)
  K%a = Ksave
  grain_mesh%r = matmul(K, grain_mesh%a)

  ! Scale back values and reset M
  grain_mesh%a = grain_mesh%a/scale_factor
  grain_mesh%r = grain_mesh%r/scale_factor**2
  M   = M/scale_factor  
  Mgb = Mgb_fac*M
  Mt  = M

  ! Compute element values
  call extract(grain_mesh%ed,grain_mesh%a,grain_mesh%enod,grain_mesh%dofnod)

  ! ! Compute flux J and jmean in each element
  ! do ie=1,grain_mesh%nelm

  !   ! Change diffusion coefficient according to grain boundary diffusion
  !   ! Find point in global mesh closest to (mean_x,mean_y)
  !   ! Mean x and y of element
  !   ! mean_x = sum(grain_mesh%coord(1,grain_mesh%enod(:,ie)))/size(grain_mesh%enod(:,ie))
  !   ! mean_y = sum(grain_mesh%coord(2,grain_mesh%enod(:,ie)))/size(grain_mesh%enod(:,ie))
  !   ! minIdx = minloc((global_mesh%coord(1,:) - mean_x)**2 + (global_mesh%coord(2,:) - mean_y)**2)
  !   ! mean_distance         = abs(lssys_ag(minIdx(1)))
  !   ! model_edge_x_distance = minval([abs(global_mesh%model_width-mean_x), abs(0d0 - mean_x)])
  !   ! mean_distance         = minval([mean_distance,model_edge_x_distance])
  !   ! Mt                    = M + (Mgb - M)*exp(-alpha*model_edge_x_distance)
  !   ! if (material.eq.3) then
  !   !   mean_x                = sum(grain_mesh%coord(1,grain_mesh%enod(:,ie)))/size(grain_mesh%enod(:,ie))
  !   !   model_edge_x_distance = minval([abs(global_mesh%model_width-mean_x), abs(0d0 - mean_x)])
  !   !   Mt                    = M + (Mgb - M)*exp(-alpha*model_edge_x_distance)
  !   ! endif

  !   call diffusion2D3j(grain_mesh%j_flux(:,ie),grain_mesh%coord(:,grain_mesh%enod(:,ie)),grain_mesh%ed(:,ie),Mt)
  ! enddo


  ! Deallocate
  call sparemove(K)
  deallocate(Ksave)
  

  return
end subroutine solve_diffusion_problem



! --- 3-node elmement routines ---

subroutine diffusion2D3kc(Ke,elcoord,M)
  ! Routine for computing the element stiffness for a 3-node element
  implicit none

  ! Intent inout
  real(dp), intent(out) :: Ke(3,3)

  ! Intent in
  real(dp), intent(in)  :: elcoord(:,:), M(2,2)

  ! Subroutine variables    
  real(dp)              :: B(2,3)
  real(dp)              :: x(3), y(3)
  real(dp)              :: detA, A

  ! Area
  x=elcoord(1,:)
  y=elcoord(2,:)
  detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detA*0.5D0

  ! B-matrix
  B(1,:) = [y(2)-y(3), y(3)-y(1), y(1)-y(2)]
  B(2,:) = [x(3)-x(2), x(1)-x(3), x(2)-x(1)]
  B      = B/(2*A)

  ! Stiffness matrix
  Ke=matmul(matmul(transpose(B),M),B)*A

  return
end subroutine diffusion2D3kc


subroutine diffusion2D3j(jvec,elcoord,ed,M)
  ! Routine for computing diffusion flux for a 3-node element
  implicit none

  ! Intent inout
  real(dp), intent(inout) :: jvec(:)

  ! Intent in
  real(dp), intent(in)  :: elcoord(:,:), M(2,2), ed(3)

  ! Subroutine variables    
  real(dp)              :: B(2,3)
  real(dp)              :: x(3), y(3)
  real(dp)              :: detA, A

  ! Area
  x=elcoord(1,:)
  y=elcoord(2,:)
  detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detA*0.5D0

  ! B-matrix
  B(1,:) = [y(2)-y(3), y(3)-y(1), y(1)-y(2)]
  B(2,:) = [x(3)-x(2), x(1)-x(3), x(2)-x(1)]
  B      = B/(2*A)

  ! Flux  
  jvec = -matmul(M,matmul(B,ed))

  return
end subroutine diffusion2D3j


subroutine lvlsetIMCarea2D3_grain(IMC_area_grain, coord, enod, nelm)
  ! Routine for computing IMC area of elm
  implicit none

  ! Intent inout
  real(dp), intent(inout)       :: IMC_area_grain

  ! Intent in 
  real(dp), intent(in)          :: coord(:,:)
  integer, intent(in)           :: enod(:,:), nelm

  ! Subroutine variables
  real(dp)              :: x(3), y(3)
  real(dp)              :: detA, A, elcoord(2,3)
  integer               :: ie
  

  IMC_area_grain = 0d0
  do ie=1,nelm

    ! Area
    elcoord = coord(:,enod(:,ie))
    x=elcoord(1,:)
    y=elcoord(2,:)
    detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
    A=detA*0.5D0

    IMC_area_grain = IMC_area_grain + A

  enddo

  return
end subroutine lvlsetIMCarea2D3_grain


! --- 4-node elmement routines ---

subroutine diffusion2D4kc(Kce,coord,M)
  ! Routine for computing the element stiffness due to curvature motion
  implicit none

  double precision, intent(in)  :: coord(:,:), M(:,:)
  double precision, intent(out) :: Kce(:,:)

  double precision              :: JT(8,2), JTinv(2,2), DNX(2,4), B(2,4), Dgp(2,2)
  double precision              :: DetJ
  integer                       :: GP_NR, indx(2)
  integer, parameter            :: NGP=4


  JT=MATMUL(DNR4,transpose(coord))
  Kce=0D0
  do GP_NR=1,NGP

    indx = (/ 2*gp_nr-1, 2*gp_nr /)
    detJ = det2(JT(indx,:))
    call inv2(JTinv,JT(indx,:))
    dNX  = matmul(JTinv,DNR4(indx,:))

    B   = dNX

    Kce = Kce+MATMUL(TRANSPOSE(B),MATMUL(M,B))*detJ

  enddo

  return
end subroutine diffusion2D4kc


subroutine diffusion2D4kc_enhanced_gb_diffusion(Kce,coord,M,model_width)
  ! Routine for computing the element stiffness due to curvature motion
  implicit none

  real(dp), intent(in)  :: coord(:,:), M(:,:), model_width
  real(dp), intent(out) :: Kce(:,:)

  real(dp)              :: JT(8,2), JTinv(2,2), DNX(2,4), B(2,4), N(4)
  real(dp)              :: DetJ
  integer               :: GP_NR, indx(2)
  integer, parameter    :: NGP=4
  real(dp)              :: DC_L, DC_gb, DC, ww, xlim, xb, gp_x, D(2,2), x0

  ! Mid of model
  x0 = model_width/2d0

  JT=MATMUL(DNR4,transpose(coord))
  Kce=0D0
  do GP_NR=1,NGP

    indx = (/ 2*gp_nr-1, 2*gp_nr /)
    detJ = det2(JT(indx,:))
    call inv2(JTinv,JT(indx,:))
    dNX  = matmul(JTinv,DNR4(indx,:))

    ! B
    B = dNX

    ! N
    N = NR4(GP_NR,:)

    ! x-coordinates of gauss point
    gp_x = dot_product(N,coord(1,:))

    ! Diffusion coefficients
    DC_L  = M(1,1)      ! mm/h
    DC_gb = 1d0*DC_L  ! mm/h

    ! Coefficient determing width of Gauss function
    ww   = 30d0
    xlim = 0.2d0*model_width ! Region of model where Gauss function is applied

    ! Distribution of diffusion coefficient
    if (abs(gp_x-x0).gt.xlim) then
      DC = DC_L
    else 
      ! Middle gb
      if (abs(gp_x-x0).lt.xlim) then
        xb = (gp_x-x0)/xlim
      endif
      DC  = (DC_L + (DC_gb - DC_L)*exp(-ww*xb**2d0))
    endif 

    ! Diffusion coefficient matrix
    D(1,:) =  (/DC,0d0/)
    D(2,:) =  (/0d0,DC/)

    Kce = Kce+MATMUL(TRANSPOSE(B),MATMUL(D,B))*detJ

  enddo

  return
end subroutine diffusion2D4kc_enhanced_gb_diffusion


subroutine diffusion2D4j(jvec,coord,ed,M)
  ! Routine for computing diffusion flux for a 4-node element
  implicit none

  real(dp), intent(inout)        :: jvec(:,:)

  double precision, intent(in)   :: coord(:,:), ed(:), M(2,2)

  double precision               :: JT(8,2), JTinv(2,2), DNX(2,4), B(2,4)
  double precision               :: DetJ
  integer                        :: GP_NR, indx(2)
  integer, parameter             :: NGP=4

  JT = MATMUL(DNR4,transpose(coord))
  do GP_NR=1,NGP

    indx = (/ 2*gp_nr-1, 2*gp_nr /)
    detJ = det2(JT(indx,:))
    call inv2(JTinv,JT(indx,:))
    dNX = matmul(JTinv,DNR4(indx,:))
    B = dNX

    ! Flux
    jvec(:,gp_nr) = -matmul(M,matmul(B,ed))

  enddo

  return
end subroutine diffusion2D4j


subroutine lvlsetIMCarea2D4(IMC_area_grain, coord, enod, nelm)
  ! Routine for computing diffusion flux for a 4-node element
  implicit none

  ! Intent inout
  real(dp), intent(inout) :: IMC_area_grain

  ! Intent in 
  real(dp), intent(in)    :: coord(:,:)
  integer, intent(in)     :: enod(:,:), nelm

  ! Subroutine variables
  double precision        :: JT(8,2)
  double precision        :: DetJ
  integer                 :: GP_NR, indx(2)
  integer, parameter      :: NGP=4
  real(dp)                :: elcoord(2,4)
  integer                 :: ie

  IMC_area_grain = 0d0
  do ie=1,nelm

    ! Area
    elcoord = coord(:,enod(:,ie))

    JT = MATMUL(DNR4,transpose(elcoord))
    do GP_NR=1,NGP

      indx = (/ 2*gp_nr-1, 2*gp_nr /)
      detJ = det2(JT(indx,:))

      ! Flux
      IMC_area_grain = IMC_area_grain + detJ
    enddo
  enddo

  return
end subroutine lvlsetIMCarea2D4


! --- Other routines ---


subroutine compute_jint(mesh, nseg, nsep_lines, sep_lines)
  ! --- Routine for computing jint in level set interface ---
  implicit none

  ! Intent inout
  type(grain_mesh_system), intent(inout) :: mesh

  ! Intent in
  integer, intent(in)                    :: nseg, nsep_lines, sep_lines(:,:)


  ! Subroutine variables
  integer                                :: sep_idx, nseg_sep, iseg_start, iseg_end, add_extra

  ! Mass matrix
  type(sparse)                           :: M
  integer, allocatable                   :: enod_m(:,:)
  integer                                :: nbnods, ierr, ie, n1, n2
  real(dp)                               :: x1, y1, x2, y2, Le, Me(2,2)


  ! Allocate jint
  allocate(mesh%jint(nseg+1),stat=ierr)

  ! Compute jint for each sepatate line
  add_extra = 0 ! need to take extra step between lines
  do sep_idx = 1,nsep_lines

    iseg_start = sep_lines(sep_idx,1)
    iseg_end   = sep_lines(sep_idx,2)
    nseg_sep   = iseg_end - iseg_start + 1

    ! Compute mass matrix
    nbnods = nseg_sep + 1
    allocate(enod_m(2,nseg_sep),stat=ierr)
    do ie=1,nseg_sep
      enod_m(:,ie) = (/ie,ie+1/)
    enddo
    call spaTopDef(M,enod_m,1)

    

    ! Assemble mass matrix
    M%a = 0d0
    do ie=1,nseg_sep

      ! Nods
      n1 = mesh%bcnod(ie,1)
      n2 = mesh%bcnod(ie+1,1)

      ! Coordinates and length of element
      x1 = mesh%coord(1,n1)
      y1 = mesh%coord(2,n1)
      x2 = mesh%coord(1,n2)
      y2 = mesh%coord(2,n2)
      Le = norm2([(x2-x1), (y2-y1)]) ! Length of elm

      ! Mass matrix
      Me(1,:) = [2d0, 1d0]
      Me(2,:) = [1d0, 2d0]
      Me      = Le/6d0*Me
      call assem(M,Me,enod_m(:,ie),1) 
    enddo

    ! Solve M*vnod = -r
    mesh%jint(iseg_start+add_extra:iseg_end+1+add_extra) = -mesh%r(mesh%bcnod(iseg_start+add_extra:iseg_end+1+add_extra,1))
    call solveq(M, mesh%jint(iseg_start+add_extra:iseg_end+1+add_extra))
    add_extra = add_extra + 1

    ! Deallocate
    call sparemove(M)    

  enddo

  return
end subroutine compute_jint


subroutine compute_jint2(grain_mesh, I_IMC)
  ! --- Routine for computing jint in level set interface ---
  implicit none

  ! Intent inout
  type(grain_mesh_system), intent(inout) :: grain_mesh

  ! Subroutine variables
  integer                                :: nseg_sep, iseg_start, iseg_end, nnods_bcseg, I_IMC

  ! Mass matrix
  type(sparse)                           :: M
  integer, allocatable                   :: enod_m(:,:)
  integer                                :: nbnods, ierr, ie, n1, n2, counter, bcval_c, bcval_cc
  real(dp)                               :: x1, y1, x2, y2, Le, Me(2,2)

  ! Number of bc
  nbnods = size(grain_mesh%bcval)

  ! Allocate jint
  print*, 'nbnods: ', nbnods
  if (allocated(grain_mesh%jint)) then
    deallocate(grain_mesh%jint)
  endif
  allocate(grain_mesh%jint(nbnods),stat=ierr)
  grain_mesh%jint = 0d0  

  ! Compute jint for each sepatate line
  counter = 1
  do while (counter.lt.nbnods)
    
    iseg_start  = counter
    bcval_c     = grain_mesh%bcval_idx(counter)
    nnods_bcseg = 0
    bcval_cc    = 1
    do while (abs(bcval_cc-bcval_c).lt.1d-15 .or. nnods_bcseg.eq.0)
        nnods_bcseg = nnods_bcseg + 1
        bcval_cc    = grain_mesh%bcval_idx(counter + nnods_bcseg-1)
    enddo
    nnods_bcseg = nnods_bcseg - 1
    iseg_end    = counter + nnods_bcseg - 1
    nseg_sep    = nnods_bcseg - 1

    print *, 'iseg_start', iseg_start
    print *, 'iseg_end', iseg_end

    if (nnods_bcseg.gt.1) then
      ! Compute mass matrix
      allocate(enod_m(2,nseg_sep),stat=ierr)
      do ie=1,nseg_sep
        enod_m(:,ie) = (/ie,ie+1/)
      enddo

      print *, 'before allocating M'
      call spaTopDef(M,enod_m,1)

      ! Assemble mass matrix      
      M%a = 0d0
      do ie=1,nseg_sep

        ! Nods
        n1 = grain_mesh%bcnod(counter + ie - 1,1)
        n2 = grain_mesh%bcnod(counter + ie,1)

        ! Coordinates and length of element
        x1 = grain_mesh%coord(1,n1)
        y1 = grain_mesh%coord(2,n1)
        x2 = grain_mesh%coord(1,n2)
        y2 = grain_mesh%coord(2,n2)
        Le = norm2([(x2-x1), (y2-y1)]) ! Length of elm

        ! Mass matrix
        Me(1,:) = [2d0, 1d0]
        Me(2,:) = [1d0, 2d0]
        Me      = Le/6d0*Me
        call assem(M,Me,enod_m(:,ie),1)
      enddo

      
      ! Solve M*vnod = -r
      ! if (I_IMC.eq.37) then
      !   print *, 'as'
      ! endif
      ! print *, 'iseg_start: ', iseg_start
      ! print *, 'iseg_end: ', iseg_end
      ! print *, 'size enodm: ', size(enod_m)
      ! print *, 'max enodm: ', maxval(enod_m)
      print *, 'size grain_mesh%r', size(grain_mesh%r)
      grain_mesh%jint(iseg_start:iseg_end) = -grain_mesh%r(grain_mesh%bcnod(iseg_start:iseg_end,1))
      call solveq(M, grain_mesh%jint(iseg_start:iseg_end))

      ! Deallocate
      ! call sparemove(M)
      if (allocated(enod_m)) then
        deallocate(enod_m)
      endif

    endif

    ! Update counter
    counter = counter + nnods_bcseg

  enddo


  ! If jint much bigger than other jint (such as at triple junction - reduce)
  ! jint_mean = sum(abs(mesh%jint))/size(mesh%jint)
  ! do ie=1,size(mesh%jint)
  !   if (abs(mesh%jint(ie)).gt.1d-10) then
  !   ! if (abs(mesh%jint(ie)).gt.jint_mean*20d0) then
  !     mesh%jint(ie) = sign(1d-11, mesh%jint(ie))
  !   endif
  ! enddo


  return
end subroutine compute_jint2
  

subroutine interpolate_glob_to_grain(global_mesh,grain_mesh,p_ed)
  ! Interpolate pq%p from global mesh to grain_mesh
  implicit none

  ! Intent inout
  type(grain_mesh_system),intent(inout) :: grain_mesh

  ! Intent in
  type(mesh_system), intent(in)         :: global_mesh
  real(dp), intent(in)                  :: p_ed(:,:)         ! Pressure field (per nod in global mesh)

  ! Subroutine variables
  integer  :: inod, ie, ierr
  real(dp) :: x, y, xecoords(4), yecoords(4), x1,x2,x3,x4,y1,y2,y3,y4, detT, c1,c2,c3, pevals(4)
  real(dp) :: xsi,eta, Nvec(4)
  logical  :: pwt


  ! Allocate the hydrostatic pressure in grain_mesh, grain_mesh%p
  if (allocated(grain_mesh%p)) then
    deallocate(grain_mesh%p)
  endif
  allocate(grain_mesh%p(grain_mesh%nnod),stat=ierr)

  ! Interpolate pq%p from global mesh to grain_mesh
  do inod=1,grain_mesh%nnod

    ! Global coordinates of point
    x = grain_mesh%coord(1,inod)
    y = grain_mesh%coord(2,inod)
    
    ! Find which element in the global mesh that the point (x,y) belong to
    pwt = .false.
    do ie=1,global_mesh%nelm
      xecoords = global_mesh%coord(1,global_mesh%enod(:,ie))
      yecoords = global_mesh%coord(2,global_mesh%enod(:,ie))

      ! Divide element in two triangles and use point in triangle function

      ! 4 --- 3      4 --- 3         3
      ! |     |      |    /       /  |
      ! |     |  ->  |(1)/       /(2)|
      ! |     |      |  /       /    |
      ! 1 --- 2      1         1 --- 2

      ! First triangle
      x1=xecoords(1); x2=xecoords(3); x3=xecoords(4)
      y1=yecoords(1); y2=yecoords(3); y3=yecoords(4)
      call point_within_triangle(pwt,x,y,x1,x2,x3,y1,y2,y3)
      if (pwt) then
        exit
      endif  
      
      ! Second triangle
      x1=xecoords(1); x2=xecoords(2); x3=xecoords(3)
      y1=yecoords(1); y2=yecoords(2); y3=yecoords(3)
      call point_within_triangle(pwt,x,y,x1,x2,x3,y1,y2,y3)
      if (pwt) then
        exit
      endif 

    enddo

    ! Element in global mesh found. Global element number is ie and node coordinates are given in xecoords and yecoords 
    ! Now the task is, given (x,y), interpolate the node values of the element to (x,y)

    ! Global element coordinates
    x1 = xecoords(1); x2 = xecoords(2); x3 = xecoords(3); x4 = xecoords(4);
    y1 = yecoords(1); y2 = yecoords(2); y3 = yecoords(3); y4 = yecoords(4);

    ! Global node values
    pevals = p_ed(:,ie)

    ! Natural coordinates (isoparametric coordinates) of (x,y). Assuming fixed parallel grid!
    xsi = 2d0*(x-x1)/(x2-x1) - 1d0
    eta = 2d0*(y-y1)/(y3-y1) - 1d0

    ! Shape functions evaluated at (xsi,eta). Nvec = [N1,N2,N3,N4]
    Nvec(1) = (1d0-xsi)*(1d0-eta)/4d0
    Nvec(2) = (1d0+xsi)*(1d0-eta)/4d0
    Nvec(3) = (1d0+xsi)*(1d0+eta)/4d0
    Nvec(4) = (1d0-xsi)*(1d0+eta)/4d0

    ! Interpolated pressure at (x,y)
    grain_mesh%p(inod) = dot_product(Nvec, pevals)

  enddo

  return
end subroutine interpolate_glob_to_grain


subroutine interpolate_scalar_mesh1_to_mesh2(mesh1_enod, mesh1_coord, mesh1_ed, mesh2_enod,mesh2_coord, mesh2_a)
  ! Interpolate scalar field mesh1_a from mesh1 to mesh2
  implicit none

  ! Intent inout
  real(dp), intent(inout)              :: mesh2_a(:)

  ! Intent in
  real(dp), intent(in)                 :: mesh1_coord(:,:), mesh2_coord(:,:)
  real(dp), intent(in)                 :: mesh1_ed(:,:)
  integer, intent(in)                  :: mesh1_enod(:,:), mesh2_enod(:,:)

  ! Subroutine variables
  integer  :: inod, ie
  real(dp) :: x, y, xecoords(4), yecoords(4), x1, x2, x3, x4, y1, y2, y3, y4, pevals(4)
  real(dp) :: xsi,eta, Nvec(4)
  logical  :: pwt
  integer  :: nelm1, nelm2, nnod1, nnod2

  !
  nelm1 = size(mesh1_enod,2)
  nelm2 = size(mesh2_enod,2)
  nnod1 = size(mesh1_coord,2)
  nnod2 = size(mesh2_coord,2)

  ! Interpolate pq%p from global mesh to grain_mesh
  do inod=1,nnod2

    ! Coordinates of point in new mesh
    x = mesh2_coord(1,inod)
    y = mesh2_coord(2,inod)
    
    ! Find which element in the global mesh that the point (x,y) belong to
    pwt = .false.
    do ie=1,nelm1
      xecoords = mesh1_coord(1,mesh1_enod(:,ie))
      yecoords = mesh1_coord(2,mesh1_enod(:,ie))

      ! Divide element in two triangles and use point in triangle function

      ! 4 --- 3      4 --- 3         3
      ! |     |      |    /       /  |
      ! |     |  ->  |(1)/       /(2)|
      ! |     |      |  /       /    |
      ! 1 --- 2      1         1 --- 2

      ! First triangle
      x1=xecoords(1); x2=xecoords(3); x3=xecoords(4)
      y1=yecoords(1); y2=yecoords(3); y3=yecoords(4)
      call point_within_triangle(pwt,x,y,x1,x2,x3,y1,y2,y3)
      if (pwt) then
        exit
      endif  
      
      ! Second triangle
      x1=xecoords(1); x2=xecoords(2); x3=xecoords(3)
      y1=yecoords(1); y2=yecoords(2); y3=yecoords(3)
      call point_within_triangle(pwt,x,y,x1,x2,x3,y1,y2,y3)
      if (pwt) then
        exit
      endif 

    enddo

    if (pwt.eqv..false.) then
      print *, 'Error in interpolating. No element found'
      call exit(0)
    endif

    ! Element in global mesh found. Global element number is ie and node coordinates are given in xecoords and yecoords 
    ! Now the task is, given (x,y), interpolate the node values of the element to (x,y)

    ! Global element coordinates
    x1 = xecoords(1); x2 = xecoords(2); x3 = xecoords(3); x4 = xecoords(4);
    y1 = yecoords(1); y2 = yecoords(2); y3 = yecoords(3); y4 = yecoords(4);

    ! Global node values
    pevals = mesh1_ed(:,ie)

    ! Natural coordinates (isoparametric coordinates) of (x,y). Assuming fixed parallel grid!
    xsi = 2d0*(x-x1)/(x2-x1) - 1d0
    eta = 2d0*(y-y1)/(y3-y1) - 1d0

    ! Shape functions evaluated at (xsi,eta). Nvec = [N1,N2,N3,N4]
    Nvec(1) = (1d0-xsi)*(1d0-eta)/4d0
    Nvec(2) = (1d0+xsi)*(1d0-eta)/4d0
    Nvec(3) = (1d0+xsi)*(1d0+eta)/4d0
    Nvec(4) = (1d0-xsi)*(1d0+eta)/4d0

    ! Interpolated pressure at (x,y)
    mesh2_a(inod) = dot_product(Nvec, pevals)

  enddo

  return
end subroutine interpolate_scalar_mesh1_to_mesh2


subroutine point_within_triangle(condition_satisfied,x,y,x1,x2,x3,y1,y2,y3)
  ! Subroutine for determining if point (x,y) is within triangle with vertices xte and yte
  implicit none

  ! Intent inout
  logical, intent(inout) :: condition_satisfied  

  ! Intent in
  real(dp), intent(in)   :: x, y, x1, x2, x3, y1, y2, y3

  ! Subroutine variables
  real(dp)               :: detT, c1, c2, c3, tol=1d-13

  ! Barycentric coordinates of (x,y) in the triangle
  detT = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3)
  c1   = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3))/detT
  c2   = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3))/detT
  c3   = 1d0 - c1 - c2

  if (((0d0-tol).le.c1 .and. c1.le.(1d0+tol)) .and. ((0d0-tol).le.c2 .and. c2.le.(1d0+tol)) .and. &
  ((0d0-tol).le.c3 .and. c3.le.(1d0+tol)) .and. ((c1+c2+c3).le.1d0)) then
    condition_satisfied = .true.
  endif

  return
end subroutine point_within_triangle


subroutine compute_IMC_area(lssys,diffsys)
  ! --- Routine for computing the total IMC area ---

  ! Intent inout
  type(diffusion_system), intent(in) :: diffsys
  type(ls_system), intent(inout) :: lssys

  ! Intent in
  ! type(mesh_system), intent(in)  :: mesh

  ! Subroutine variables
  integer                        :: g
  real(dp)                       :: IMC_area_grain


  ! Compute area of IMC grains from grain meshes
  lssys%IMC_area = 0d0
  do g=1,lssys%ngrains
      if (lssys%material(g).eq.2) then
      ! IMC grain
      call lvlsetIMCarea2D3_grain(IMC_area_grain, diffsys%grain_meshes(g)%coord, diffsys%grain_meshes(g)%enod, &
      diffsys%grain_meshes(g)%nelm)
      lssys%IMC_area = lssys%IMC_area + IMC_area_grain
      endif
  enddo  

  return
end subroutine compute_IMC_area


subroutine print_diffusion_mesh_info(mesh)
  implicit none

  ! Intent inout
  type(grain_mesh_system), intent(inout) :: mesh

  ! ! Check that enough memory is allocated
  ! if ((mesh%nalloc_nelm.lt.mesh%nelm) .or. (mesh%nalloc_nnod.lt.mesh%nnod) .or. (mesh%nalloc_nbc.lt.mesh%nbc)) then      
  !    print *, 'Not allocating enough size in diffusion problem'
  !    call exit(0)
  ! endif

  write(*,'(A,A18)') ' grain mesh data: '
  write(*,'(A40,I5)')  'Number of nods in grain                :', mesh%nnod
  write(*,'(A40,I5)')  'Number of elements in grain            :', mesh%nelm  


  return
end subroutine print_diffusion_mesh_info

subroutine reinit_level_set_spatial2(a_g,line_ex,line_ey,closest_line,nseg,coord)
  ! --- Routine for reinitializing level set function in deformed configuration ---
  implicit none

  ! Intent inout
  real(dp),intent(inout)         :: a_g(:)
  integer, intent(inout)         :: closest_line(:)

  ! Intent in   
  integer, intent(in)            :: nseg
  real(dp), intent(in)           :: coord(:,:), line_ex(:,:), line_ey(:,:)

  ! Subroutine variables
  integer                        :: inod, iIntElm, ierr
  real(dp)                       :: P(2), A(2), B(2), E(2),v(2), u(2), t, seg_dists(nseg)
  ! real(dp), allocatable          :: seg_dists(nseg)


  ! Loop through all nodes. For each node -> compute distance to all line segments
  do inod = 1,size(coord,2)

      ! Coordinates of node
      P = coord(:,inod)

      ! Compute distance to all line segments for given node P
      seg_dists  = 0d0
      do iIntElm = 1,nseg

          A = [line_ex(iIntElm,1), line_ey(iIntElm,1)]
          B = [line_ex(iIntElm,2), line_ey(iIntElm,2)]
          
          v = B - A
          u = A - P
          t = -dot_product(v, u)/dot_product(v, v)
                  
          ! Determine point D on line segment closest to P
          if (t.le.0d0) then
              E = A
          elseif (t.ge.1d0) then
              E = B
          else
              E = (1d0-t)*A + t*B
          endif
          
          ! Compute closest distance from P to line segment (to point D)
          seg_dists(iIntElm) = norm2(E - P)

          if (isnan(norm2(E - P))) then
              print *, 'distance is nan'
              ! call exit(0)
          endif
      enddo


      ! Signed distance function
      a_g(inod) = sign(1d0,a_g(inod))*abs(minval(seg_dists))

      ! if (abs(a_g(inod)).lt.1d-9) then
      !   iIntElm = 1
      !   A = [line_ex(iIntElm,1), line_ey(iIntElm,1)]
      !   B = [line_ex(iIntElm,2), line_ey(iIntElm,2)]
      !   v = B - A
      !   u = A - P
      !   t = -dot_product(v, u)/dot_product(v, v)
      !   ! print *, 'as'
      ! endif

      ! Store closest line
      closest_line(inod) = minloc(seg_dists, 1)
  enddo

  return
end subroutine reinit_level_set_spatial2


SUBROUTINE FindUniqueIntegers(arr, uniqueArr)
  INTEGER, INTENT(IN)               :: arr(:)      ! Input integer array
  INTEGER, ALLOCATABLE, INTENT(OUT) :: uniqueArr(:)  ! Output array for unique integers
  INTEGER, ALLOCATABLE              :: uniqueArr2(:)  ! Output array for unique integers
  INTEGER :: i, j, counter


  ALLOCATE(uniqueArr2(size(arr)))
  uniqueArr2 = 0

  ! Initialize count to zero
  uniqueArr2(1) = arr(1)
  counter       = 1
  
  ! Check each element in the input array
  DO i = 2, size(arr)
    ! Check if the current element is unique
    IF (all((uniqueArr2(1:counter) - arr(i)) .ne. 0)) THEN
      ! Increment count and add the unique element to the output array
      counter = counter + 1
      uniqueArr2(counter) = arr(i)
    END IF
  END DO

  ! Allocate the output array with the correct size
  ALLOCATE(uniqueArr(counter))
  uniqueArr = uniqueArr2(1:counter)
  deallocate(uniqueArr2)

  return

END SUBROUTINE FindUniqueIntegers


end module diffusion
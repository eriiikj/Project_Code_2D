module json_read

  use mf_datatypes
  use json_module
  use ls_types
  use diffusion_types

  implicit none


  type(json_value), pointer :: output_root
  type(json_core)           :: jsonoutc

  character(len=40)         :: output_filename

  private jsonoutc, output_root, output_filename


  interface json_write
    module procedure write_json_real, write_json_int, write_json_realvec, write_json_intvec, write_json_realmat, write_json_intmat
  end interface json_write

contains

subroutine read_json_input_location(input_location)
  implicit none


  ! Path
  character(len=255)                          :: filename, cwd

  ! json
  type(json_file)                             :: my_json_file

  ! Other
  logical                                     :: found, file_e

  ! Out
  character(len=:),allocatable, intent(inout) :: input_location

  filename = 'python_input_location.json'
  

  call getcwd(cwd)           ! Save path

  ! ----- Locate file -----
  inquire(file=filename, exist=file_e)
  if ( .not. file_e) then
    print*, 'Input location file does not exist'
    call exit(0)
  endif

  ! ----- Initiate and load file -----
  call my_json_file%initialize()
  call my_json_file%load_file(filename); if (my_json_file%failed()) stop

  ! ----- Input location ----
  call my_json_file%get('input_location', input_location, found); if (.not. found) call exit(0)
  write(*,'(A20,A)') 'input location    : ', input_location

  ! ----- Destroy -----
  call my_json_file%destroy()

  return
end subroutine read_json_input_location


! subroutine read_json_input(input_location, imc_bound, IMC_vol_transf, IMC_steps, load_steps)
!   implicit none

!   ! Input location
!   character(len=:),allocatable, intent(in) :: input_location

!   ! Path
!   character(len=255)                       :: filename, cwd

!   ! json
!   type(json_file)                          :: my_json_file

!   ! Other
!   logical                                  :: found, file_e

!   real(dp)                                 :: Cu_height, IMC_side_height_start, IMC_mid_height_start, IMC_radius_inc, interface_w
!   real(dp)                                 :: IMC_final_vol_per_area

!   ! Out
!   type(IMC_boundary_system), intent(inout) :: imc_bound
!   real(dp), intent(out)                    :: IMC_vol_transf

!   integer , intent(inout)                  :: IMC_steps, load_steps

!   filename = 'python_input.json'
  
!   ! ----- Locate file -----
!   call getcwd(cwd)           ! Save current location
!   call chdir(input_location)

!   inquire(file=filename, exist=file_e)
!   if (file_e) then
!     print*, 'Input data file does exist'
!   else
!     print*, 'Input data file does not exist'
!     call exit(0)
!   endif

!   ! ----- Initiate and load file -----
!   call my_json_file%initialize()
!   call my_json_file%load_file(filename); if (my_json_file%failed()) stop

!   ! --- Cu_height ---
!   call my_json_file%get('y3', Cu_height, found); if (.not. found) call exit(0)

!   ! --- IMC_side_height_start ---
!   call my_json_file%get('imc_side_height_start', IMC_side_height_start, found); if (.not. found) call exit(0)

!   ! --- IMC_mid_height_start ---
!   call my_json_file%get('imc_mid_height_start', IMC_mid_height_start, found); if (.not. found) call exit(0)

!   ! --- IMC_radius_inc ---
!   call my_json_file%get('imc_radius_inc', IMC_radius_inc, found); if (.not. found) call exit(0)

!   ! --- Interface width ---
!   call my_json_file%get('interface_w', interface_w, found); if (.not. found) call exit(0)

!   ! --- IMC volume per area ---
!   call my_json_file%get('imc_final_vol_per_area', IMC_final_vol_per_area, found); if (.not. found) call exit(0)

!   ! --- Convert to mm ---
!   imc_bound%Cu_height              = Cu_height*1d-3
!   imc_bound%IMC_side_height_start  = IMC_side_height_start*1d-3
!   imc_bound%IMC_mid_height_start   = IMC_mid_height_start*1d-3
!   imc_bound%IMC_radius_inc         = IMC_radius_inc*1d-3
!   imc_bound%interface_w            = interface_w*1d-3
!   imc_bound%IMC_final_vol_per_area = IMC_final_vol_per_area*1d-3


!   ! ----- IMC_vol_transf -----
!   call my_json_file%get('imc_vol_transf', IMC_vol_transf, found); if (.not. found) call exit(0)

!   ! ----- IMC_steps -----
!   call my_json_file%get('imc_steps', IMC_steps, found); if (.not. found) call exit(0)

!   ! ----- load_steps -----
!   call my_json_file%get('load_steps', load_steps, found); if (.not. found) call exit(0)

!   ! ----- Destroy -----
!   call my_json_file%destroy()

!   ! ----- Go back to program -----
!   call chdir(cwd)

!   return
! end subroutine read_json_input

subroutine read_json_input_2(input_location, IMC_vol_transf, IMC_steps, lssys, diffsys)
  implicit none

  ! Input location
  character(len=:),allocatable, intent(in) :: input_location

  ! Path
  character(len=255)                       :: filename, cwd

  ! json
  type(json_file)                          :: my_json_file

  ! Other
  logical                                  :: found, file_e

  real(dp), intent(out)                    :: IMC_vol_transf
  integer , intent(inout)                  :: IMC_steps
  type(ls_system), intent(inout)           :: lssys
  type(diffusion_system), intent(inout)    :: diffsys  

  filename = 'python_input.json'  
  
  ! ----- Locate file -----
  call getcwd(cwd)           ! Save current location
  call chdir(input_location)

  inquire(file=filename, exist=file_e)
  if (.not. file_e) then
    print*, 'Input data file does not exist'
    call exit(0)
  endif

  ! ----- Initiate and load file -----
  call my_json_file%initialize()
  call my_json_file%load_file(filename); if (my_json_file%failed()) stop


  ! ----- IMC_vol_transf -----
  call my_json_file%get('imc_vol_transf', IMC_vol_transf, found); if (.not. found) call exit(0)

  ! ----- IMC_steps -----
  call my_json_file%get('imc_steps', IMC_steps, found); if (.not. found) call exit(0)

  ! ----- ls gamma -----
  call my_json_file%get('ls_gamma', lssys%gamma, found); if (.not. found) call exit(0)

  ! ----- ls m -----
  call my_json_file%get('ls_m', lssys%m, found); if (.not. found) call exit(0)

  ! ----- ls zeta -----
  call my_json_file%get('ls_zeta', lssys%zeta, found); if (.not. found) call exit(0)

  ! ----- ls mzetagb -----
  call my_json_file%get('ls_mzetagb', lssys%mzetagb, found); if (.not. found) call exit(0)

  ! ----- ls phiw -----
  call my_json_file%get('ls_phiw', lssys%w, found); if (.not. found) call exit(0)

  ! ----- ls malpha -----
  call my_json_file%get('ls_malpha', lssys%malpha, found); if (.not. found) call exit(0)

  ! ----- diff DIMC -----
  call my_json_file%get('ls_DIMC', diffsys%D_diff_material(2), found); if (.not. found) call exit(0)

  ! ----- Destroy -----
  call my_json_file%destroy()

  ! ----- Go back to program -----
  call chdir(cwd)

  return
end subroutine read_json_input_2

! subroutine read_default_input(imc_bound, IMC_vol_transf, IMC_steps, load_steps)
!   implicit none

!   ! Out
!   type(IMC_boundary_system), intent(inout) :: imc_bound
!   real(dp), intent(out)                    :: IMC_vol_transf

!   integer , intent(inout)                  :: IMC_steps, load_steps

!    ! --- Input data - standard values ---

!    ! Define default values for IMC boundary increase
!   imc_bound%Cu_height             = 1d-3     ! Height of Cu-layer
!   imc_bound%IMC_side_height_start = 0.03d-3  ! IMC start height at ends, must be smaller than grain_width/2 (all grains same width)
!   imc_bound%IMC_mid_height_start  = 0.02d-3  ! IMC mid height (center of grain) at start
!   imc_bound%IMC_radius_inc        = 0.002d-3 ! IMC radius reduction by increment
!   imc_bound%interface_w           = 0.4d-3   ! IMC radius reduction by increment

!   ! Define default value of IMC volume increase
!   IMC_vol_transf = 1.1d0

!   ! Define standard values of IMC load loop
!   IMC_steps   = 2
!   load_steps  = 5
  
!   return
! end subroutine read_default_input



subroutine read_json_mesh(input_location, coord, newcoord, enod, bcnod, bcval, bcnod_all)
  implicit none

  ! Input location
  character(len=:),allocatable, intent(in) :: input_location

  ! Path
  character(len=255)                       :: filename, cwd

  ! json
  type(json_file)                          :: my_json_file
  type(json_core)                          :: core
  type(json_value),pointer                 :: readpointer,child

  ! Read variables
  real(dp), allocatable                    :: coordrow(:)
  integer,  allocatable                    :: enodrow(:), bcnodrow(:)

  ! Other
  logical                                  :: found, file_e
  integer                                  :: rows, cols, i, var_type

  ! Out
  real(dp), allocatable, intent(inout)     :: coord(:,:), newcoord(:,:), bcval(:)
  integer , allocatable, intent(inout)     :: enod(:,:), bcnod(:,:), bcnod_all(:)

  filename = 'python_mesh.json'
  
  ! ----- Locate file -----
  call getcwd(cwd)           ! Save path
  call chdir(input_location)

  inquire(file=filename, exist=file_e)
  if (.not. file_e) then
    print*, 'Input mesh file does not exist'
    call exit(0)
  endif

  ! ----- Initiate and load file -----
  call my_json_file%initialize()
  call my_json_file%load_file(filename); if (my_json_file%failed()) stop

  ! ----- Coord -----
  call my_json_file%info('coord',found,var_type,rows)
  call my_json_file%info('coord(1)',found,var_type,cols)

  allocate(coord(rows,cols))
  allocate(newcoord(rows,cols))

  ! Pointer to the coord matrix:
  call my_json_file%get('coord',readpointer)

  ! Read coord and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,coordrow) 
    if (size(coordrow)/=cols) error stop 'error: column is wrong size'
    coord(i,:) = coordrow
  end do  

  ! ----- enod -----
  call my_json_file%info('enod',found,var_type,rows)
  call my_json_file%info('enod(1)',found,var_type,cols)

  allocate(enod(rows,cols))

  ! Pointer to the enod matrix:
  call my_json_file%get('enod',readpointer)

  ! Read enod and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,enodrow) 
    if (size(enodrow)/=cols) error stop 'error: column is wrong size'
    enod(i,:) = enodrow
  end do

  ! ----- bcnod -----
  call my_json_file%info('bcnod',found,var_type,rows)
  call my_json_file%info('bcnod(1)',found,var_type,cols)

  allocate(bcnod(rows,cols))

  ! Pointer to the bcnod matrix:
  call my_json_file%get('bcnod',readpointer)

  ! Read bcnod and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,bcnodrow) 
    if (size(bcnodrow)/=cols) error stop 'error: column is wrong size'
    bcnod(i,:) = bcnodrow
  end do


  ! ----- bcval -----
  call my_json_file%get('bcval', bcval, found); if (.not. found) call exit(0)  
  
  ! ----- bcnod_all -----
  call my_json_file%get('bcnod_all', bcnod_all, found); if (.not. found) call exit(0)

  ! ----- Destroy -----
  call my_json_file%destroy()

  ! --- Convert to mm ---
  coord    = coord*1d-3
  newcoord = coord

  ! ----- Go back to program -----
  call chdir(cwd)

  return
end subroutine read_json_mesh


subroutine read_json_phase_mesh(input_location, i_IMC, nelm, nnod, nbc, coord, enod, bcnod, bcval, grain)
  implicit none

  ! Input location
  character(len=:),allocatable, intent(in) :: input_location
  integer, intent(in)                      :: i_IMC, grain

  ! Path
  character(len=255)                       :: filename, cwd

  ! json
  type(json_file)                          :: my_json_file
  type(json_core)                          :: core
  type(json_value),pointer                 :: readpointer,child

  ! Read variables
  real(dp), allocatable                    :: coordrow(:)
  integer,  allocatable                    :: enodrow(:), bcnodrow(:)

  ! Other
  logical                                  :: found, file_e
  integer                                  :: rows, cols, i, var_type

  real(dp),allocatable :: bcval2(:)
  ! Out
  integer, intent(inout)                   :: nelm, nnod, nbc
  real(dp), intent(inout)                  :: coord(:,:), bcval(:)
  integer, intent(inout)                   :: enod(:,:), bcnod(:,:)  
  

  ! Filename phase mesh file  
  if (i_IMC.lt.10) then
    write(filename, '(A,I1,A,I1,A)' ), 'phase_mesh_', i_IMC, '_', grain, '.json'
  elseif (i_IMC.ge.10 .and. i_IMC.lt.100) then
    write(filename, '(A,I2,A,I1,A)' ), 'phase_mesh_', i_IMC, '_', grain, '.json'
  elseif (i_IMC.ge.100 .and. i_IMC.lt.1000) then
    write(filename, '(A,I3,A,I1,A)' ), 'phase_mesh_', i_IMC, '_', grain, '.json'
  endif



  ! ----- Locate file -----
  call getcwd(cwd)           ! Save path
  call chdir(input_location)

  inquire(file=trim(filename), exist=file_e)
  if (.not. file_e) then
    print*, 'Phase mesh file does not exist'
    call exit(0)
  endif

  ! ----- Initiate and load file -----
  call my_json_file%initialize()
  call my_json_file%load_file(trim(filename)); if (my_json_file%failed()) stop

  ! ----- Coord -----
  call my_json_file%info('coord',found,var_type,rows)
  call my_json_file%info('coord(1)',found,var_type,cols)

  ! Pointer to the coord matrix:
  call my_json_file%get('coord',readpointer)

  ! Read coord and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,coordrow) 
    if (size(coordrow)/=cols) error stop 'error: column is wrong size'
    coord(i,1:size(coordrow)) = coordrow
  end do  

  nnod = size(coordrow)

  ! ----- enod -----
  call my_json_file%info('enod',found,var_type,rows)
  call my_json_file%info('enod(1)',found,var_type,cols)

  ! allocate(enod(rows,cols))

  ! Pointer to the enod matrix:
  call my_json_file%get('enod',readpointer)

  ! Read enod and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,enodrow) 
    if (size(enodrow)/=cols) error stop 'error: column is wrong size'
    enod(i,1:size(enodrow)) = enodrow
  end do

  nelm = size(enodrow)

  ! ----- bcnod -----
  call my_json_file%info('bcnod',found,var_type,rows)
  call my_json_file%info('bcnod(1)',found,var_type,cols)

  ! Pointer to the bcnod matrix:
  call my_json_file%get('bcnod',readpointer)

  ! Read bcnod and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,bcnodrow) 
    if (size(bcnodrow)/=cols) error stop 'error: column is wrong size'
    bcnod(i,1:size(bcnodrow)) = bcnodrow
  end do

  nbc = rows

  ! ----- bcval -----
  call my_json_file%get('bcval', bcval2, found); if (.not. found) call exit(0)
  bcval(1:nbc) = bcval2

  ! ----- Destroy -----
  call my_json_file%destroy()

  ! --- Convert to mm ---
  coord = coord*1d-3

  ! ----- Go back to program -----
  call chdir(cwd)

  return
end subroutine read_json_phase_mesh



subroutine read_json_diffusion_mesh(input_location, i_IMC, coord, enod)
  implicit none

  ! Input location
  character(len=:),allocatable, intent(in) :: input_location
  integer, intent(in)                      :: i_IMC

  ! Path
  character(len=255)                       :: filename, cwd

  ! json
  type(json_file)                          :: my_json_file
  type(json_core)                          :: core
  type(json_value),pointer                 :: readpointer,child

  ! Read variables
  real(dp), allocatable                    :: coordrow(:)
  integer,  allocatable                    :: enodrow(:)

  ! Other
  logical                                  :: found, file_e
  integer                                  :: rows, cols, i, var_type  

  ! Out
  real(dp), allocatable, intent(inout)     :: coord(:,:)
  integer, allocatable, intent(inout)      :: enod(:,:)
  

  ! Filename phase mesh file  
  if (i_IMC.lt.10) then
    write(filename, '(A,I1,A,I1,A)' ), 'triangle_mesh_', i_IMC, '.json'
  elseif (i_IMC.ge.10 .and. i_IMC.lt.100) then
    write(filename, '(A,I2,A,I1,A)' ), 'triangle_mesh_', i_IMC, '.json'
  elseif (i_IMC.ge.100 .and. i_IMC.lt.1000) then
    write(filename, '(A,I3,A,I1,A)' ), 'triangle_mesh_', i_IMC, '.json'
  endif

  ! ----- Locate file -----
  call getcwd(cwd)           ! Save path
  call chdir(input_location)

  inquire(file=trim(filename), exist=file_e)
  if (.not. file_e) then
    print*, 'Phase mesh file does not exist'
    call exit(0)
  endif

  ! ----- Initiate and load file -----
  call my_json_file%initialize()
  call my_json_file%load_file(trim(filename)); if (my_json_file%failed()) stop

  ! ----- Coord -----
  call my_json_file%info('coord',found,var_type,rows)
  call my_json_file%info('coord(1)',found,var_type,cols)  

  ! Pointer to the coord matrix:
  allocate(coord(rows,cols))
  call my_json_file%get('coord',readpointer)

  ! Read coord and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,coordrow) 
    if (size(coordrow)/=cols) error stop 'error: column is wrong size'
    coord(i,1:size(coordrow)) = coordrow
  end do

  ! ----- enod -----
  call my_json_file%info('enod',found,var_type,rows)
  call my_json_file%info('enod(1)',found,var_type,cols)  

  ! Pointer to the enod matrix:  
  allocate(enod(rows,cols))
  call my_json_file%get('enod',readpointer)

  ! Read enod and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,enodrow) 
    if (size(enodrow)/=cols) error stop 'error: column is wrong size'
    enod(i,1:size(enodrow)) = enodrow
  end do

  ! ----- Destroy -----
  call my_json_file%destroy()

  ! --- Convert to mm ---
  coord = coord*1d-3

  ! ----- Go back to program -----
  call chdir(cwd)

  return
end subroutine read_json_diffusion_mesh

subroutine read_json_phase_mesh2(input_location, i_IMC, coord, enod, bcnod, bcval, bcval_idx, grain)
  implicit none

  ! Input location
  character(len=:),allocatable, intent(in) :: input_location
  integer, intent(in)                      :: i_IMC, grain

  ! Path
  character(len=255)                       :: filename, cwd

  ! json
  type(json_file)                          :: my_json_file
  type(json_core)                          :: core
  type(json_value),pointer                 :: readpointer,child

  ! Read variables
  real(dp), allocatable                    :: coordrow(:)
  integer,  allocatable                    :: enodrow(:), bcnodrow(:)

  ! Other
  logical                                  :: found, file_e
  integer                                  :: rows, cols, i, var_type  

  ! Out
  real(dp), allocatable, intent(inout)     :: coord(:,:), bcval(:), bcval_idx(:)
  integer, allocatable, intent(inout)      :: enod(:,:), bcnod(:,:)  
  

  ! Filename phase mesh file  
  if (i_IMC.lt.10) then
    write(filename, '(A,I1,A,I1,A)' ), 'phase_mesh_', i_IMC, '_', grain, '.json'
  elseif (i_IMC.ge.10 .and. i_IMC.lt.100) then
    write(filename, '(A,I2,A,I1,A)' ), 'phase_mesh_', i_IMC, '_', grain, '.json'
  elseif (i_IMC.ge.100 .and. i_IMC.lt.1000) then
    write(filename, '(A,I3,A,I1,A)' ), 'phase_mesh_', i_IMC, '_', grain, '.json'
  endif

  ! ----- Locate file -----
  call getcwd(cwd)           ! Save path
  call chdir(input_location)

  inquire(file=trim(filename), exist=file_e)
  if (.not. file_e) then
    print*, 'Phase mesh file does not exist'
    call exit(0)
  endif

  ! ----- Initiate and load file -----
  call my_json_file%initialize()
  call my_json_file%load_file(trim(filename)); if (my_json_file%failed()) stop

  ! ----- Coord -----
  call my_json_file%info('coord',found,var_type,rows)
  call my_json_file%info('coord(1)',found,var_type,cols)  

  ! Pointer to the coord matrix:
  allocate(coord(rows,cols))
  call my_json_file%get('coord',readpointer)

  ! Read coord and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,coordrow) 
    if (size(coordrow)/=cols) error stop 'error: column is wrong size'
    coord(i,1:size(coordrow)) = coordrow
  end do

  ! ----- enod -----
  call my_json_file%info('enod',found,var_type,rows)
  call my_json_file%info('enod(1)',found,var_type,cols)  

  ! Pointer to the enod matrix:  
  allocate(enod(rows,cols))
  call my_json_file%get('enod',readpointer)

  ! Read enod and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,enodrow) 
    if (size(enodrow)/=cols) error stop 'error: column is wrong size'
    enod(i,1:size(enodrow)) = enodrow
  end do

  ! ----- bcnod -----
  call my_json_file%info('bcnod',found,var_type,rows)
  call my_json_file%info('bcnod(1)',found,var_type,cols)

  ! Pointer to the bcnod matrix:
  allocate(bcnod(rows,cols))
  call my_json_file%get('bcnod',readpointer)

  ! Read bcnod and store in matrix
  do i=1,rows
    call core%get_child(readpointer,i,child)
    if (.not. associated(child)) error stop 'error: column not found'
    call core%get(child,bcnodrow) 
    if (size(bcnodrow)/=cols) error stop 'error: column is wrong size'
    bcnod(i,1:size(bcnodrow)) = bcnodrow
  end do


  ! ----- bcval -----
  call my_json_file%get('bcval', bcval, found); if (.not. found) call exit(0)

  ! ----- bcval_idx -----
  call my_json_file%get('bcval_idx', bcval_idx, found); if (.not. found) call exit(0)

  ! ----- Destroy -----
  call my_json_file%destroy()

  ! --- Convert to mm ---
  coord = coord*1d-3

  ! ----- Go back to program -----
  call chdir(cwd)

  return
end subroutine read_json_phase_mesh2


subroutine read_json_init_level_set(input_location, a_ls)
  implicit none

  ! Input location
  character(len=:),allocatable, intent(in) :: input_location

  ! Path
  character(len=255)                       :: filename, cwd

  ! json
  type(json_file)                          :: my_json_file
  logical                                  :: found, file_e

  ! Out
  real(dp), allocatable, intent(inout)     :: a_ls(:)
  

  ! Filename phase mesh file
  write(filename, '(A,I,A)' ), 'level_set_init.json'

  ! ----- Locate file -----
  call getcwd(cwd)           ! Save path
  call chdir(input_location)

  inquire(file=filename, exist=file_e)
  if (.not. file_e) then
    print*, 'Phase mesh file does not exist'
    call exit(0)
  endif

  ! ----- Initiate and load file -----
  call my_json_file%initialize()
  call my_json_file%load_file(filename); if (my_json_file%failed()) stop

  ! ----- Get a_ls -----
  call my_json_file%get('a_ls', a_ls, found); if (.not. found) call exit(0)

  ! ----- Destroy -----
  call my_json_file%destroy()

  ! ----- Go back to program -----
  call chdir(cwd)

  return
end subroutine read_json_init_level_set


subroutine json_write_open(filename)
  implicit none

  character(len=*), intent(in) :: filename

  output_filename = filename

  call jsonoutc%create_object(output_root,'')  ! create root object

  return
end subroutine json_write_open



subroutine json_write_close()
  implicit none

  call jsonoutc%destroy(output_root)

  return
end subroutine json_write_close



subroutine write_json_real(real, real_name)
  implicit none

  character(len=*), intent(in) :: real_name
  real(dp), intent(in)         :: real

  ! --- Add real ---
  call jsonoutc%add(output_root,real_name,real)
  
  ! --- Print to file ---
  call jsonoutc%print(output_root,output_filename)

  return
end subroutine write_json_real



subroutine write_json_int(int, int_name)
  implicit none

  character(len=*), intent(in) :: int_name
  integer, intent(in)          :: int

  ! --- Add int ---
  call jsonoutc%add(output_root,int_name,int)
  
  ! --- Print to file ---
  call jsonoutc%print(output_root,output_filename)

  return
end subroutine write_json_int



subroutine write_json_realvec(realvec, realvec_name)
  implicit none

  character(len=*), intent(in) :: realvec_name
  real(dp), intent(in)         :: realvec(:)

  ! --- Add realvec ---
  call jsonoutc%add(output_root,realvec_name,realvec)
  
  ! --- Print to file ---
  call jsonoutc%print(output_root,output_filename)

  return
end subroutine write_json_realvec



subroutine write_json_intvec(intvec, intvec_name)
  implicit none

  character(len=*), intent(in) :: intvec_name
  integer, intent(in)          :: intvec(:)

  ! --- Add intvec ---
  call jsonoutc%add(output_root,intvec_name,intvec)
  
  ! --- Print to file ---
  call jsonoutc%print(output_root,output_filename)

  return
end subroutine write_json_intvec



subroutine write_json_realmat(realmat, realmat_name)
  implicit none

  character(len=*), intent(in) :: realmat_name
  real(dp), intent(in)         :: realmat(:,:)

  type(json_value),pointer     :: mat_pointer

  integer :: rows, i

  ! --- Add mat ---
  call jsonoutc%create_array(mat_pointer,realmat_name)
  rows = size(realmat,1) 
  do i=1,rows
    call jsonoutc%add(mat_pointer,'',realmat(i,:))
  enddo

  call jsonoutc%add(output_root,mat_pointer) !add mat_pointer array to the root structure
  
  ! --- Print to file ---
  call jsonoutc%print(output_root,output_filename)

  return
end subroutine write_json_realmat


subroutine write_json_intmat(intmat, intmat_name)
  implicit none

  character(len=*), intent(in) :: intmat_name
  integer, intent(in)          :: intmat(:,:)

  type(json_value),pointer     :: mat_pointer

  integer :: rows, i

  ! --- Add mat ---
  call jsonoutc%create_array(mat_pointer,intmat_name)
  rows = size(intmat,1) 
  do i=1,rows
    call jsonoutc%add(mat_pointer,'',intmat(i,:))
  enddo

  call jsonoutc%add(output_root,mat_pointer) !add mat_pointer array to the root structure
  
  ! --- Print to file ---
  call jsonoutc%print(output_root,output_filename)

  return
end subroutine write_json_intmat
  

end module json_read

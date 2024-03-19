module mesh_module
  
  ! OpenMP
  use omp_lib

  ! Mflib
  use mf_datatypes

  ! Somelib
  use fem_system

  ! Whiskerlib
  use gp_util  

  implicit none

  type mesh_system
    integer,  allocatable :: enod(:,:),  bcnod1(:,:), bcdof1(:), bcnod(:,:), bcdof(:), bcnod_all(:)
    integer, allocatable  :: bcnods_left_side(:), bcnods_right_side(:)
    real(dp), allocatable :: coord(:,:), newcoord(:,:), bcval1(:), bcval(:)
    real(dp), allocatable :: ex(:,:), ey(:,:), newex(:,:), newey(:,:)
    integer               :: dofnod, nelm, enlmx, nelmx, nelmy, nrgp, nodel, nnod, ndof, dofel, nnodgp
    real(dp)              :: model_width, model_height, elmsize_x, elmsize_y, min_elm_width, axisbc(4)
    integer               :: nset
    integer, allocatable  :: element_set(:), set_groups(:,:), elmneighb(:,:), counter_set(:)
    integer, allocatable  :: elmidx(:), material_nod(:)
    real(dp), allocatable :: gpx(:,:), gpy(:,:)

    ! Level set subdomain, grain 1
    integer, allocatable  :: enod_ls(:,:), nods_ls(:), elms_ls(:)
    real(dp), allocatable :: coord_ls(:,:), newcoord_ls(:,:), ex_ls(:,:), ey_ls(:,:), newex_ls(:,:), newey_ls(:,:)
    real(dp)              :: P0_coord(2)
    integer               :: dofnod_ls, nnod_ls, nelm_ls, P0_nod
    integer, allocatable  :: element_set_ls(:), set_groups_ls(:,:), counter_set_ls(:)
  end type mesh_system
  
contains


  subroutine init_mesh_system(mesh)
    implicit none

    ! Intent inout
    type(mesh_system), intent(inout) :: mesh

    ! Intent in

    ! Subroutine variables    
    integer                          :: ierr, k
    
    ! Extract nod topology
    mesh%dofnod_ls = 1
    mesh%dofnod    = 2
    mesh%nelm      = size(mesh%enod,2)
    mesh%nodel     = size(mesh%enod,1)
    mesh%nnod      = maxval(mesh%enod)
    mesh%ndof      = mesh%nnod*mesh%dofnod
    mesh%dofel     = mesh%nodel*mesh%dofnod    
    
    write(*,'(A17)') 'Global mesh data:'
    write(*,'(A31,I5)') 'Number of nods                :',mesh%nnod
    write(*,'(A31,I5)') 'Number of dofs                :',mesh%ndof
    write(*,'(A31,I5)') 'Number of elements            :',mesh%nelm
    write(*,'(A31,I5)') 'Number of nods in a element   :',mesh%nodel
    write(*,'(A31,I5)') 'Number of dofs in a element   :',mesh%dofel
    write(*,'(A31,I5)') 'Number of dofs in a node      :',mesh%dofnod
    write(*,'(A31,I5)') 'Number of boundary conditions :',size(mesh%bcval)    

    ! Define nbr of Gauss points
    mesh%nrgp = 4 

    ! Define nbr of nodes in gp mesh
    mesh%nnodgp = mesh%nelm*mesh%nrgp

    ! Allocate ex and ey
    allocate(mesh%ex(4,mesh%nelm), stat=ierr)
    allocate(mesh%ey(4,mesh%nelm), stat=ierr)
    allocate(mesh%newex(4,mesh%nelm), stat=ierr)
    allocate(mesh%newey(4,mesh%nelm), stat=ierr)
    
    ! Extract ex and ey
    call coordxtr(mesh%ex,mesh%ey,mesh%coord,mesh%enod)
    mesh%newex= mesh%ex
    mesh%newey= mesh%ey

    ! Allocate gpx and gpy
    allocate(mesh%gpx(mesh%nrgp,mesh%nelm), stat=ierr)
    allocate(mesh%gpy(mesh%nrgp,mesh%nelm), stat=ierr)

    ! Extract gpx and gpy
    call gpxtr(mesh%gpx, mesh%gpy, mesh%coord, mesh%enod)    

    
    ! Nelm in x and y dir
    mesh%model_width   = maxval([mesh%ex])
    mesh%model_height  = maxval([mesh%ey])
    mesh%elmsize_x     = abs(mesh%ex(2,1)-mesh%ex(1,1))
    mesh%elmsize_y     = abs(mesh%ey(3,1)-mesh%ey(1,1))
    mesh%min_elm_width = minval(maxval(mesh%ex,1) - minval(mesh%ex,1))
    mesh%nelmx         = nint(mesh%model_width/mesh%elmsize_x)
    mesh%nelmy         = nint(mesh%model_height/mesh%elmsize_y)

    ! Get element sets
    call get_elm_sets(mesh)

    ! axisbc
    mesh%axisbc = [minval(mesh%ex), maxval(mesh%ex), minval(mesh%ey), maxval(mesh%ey)]

    ! Elmidx
    allocate(mesh%elmidx(mesh%nelm),stat=ierr)
    do k=1,mesh%nelm
      mesh%elmidx(k) = k
    enddo
    
    return
  end subroutine init_mesh_system


  subroutine get_elm_sets(mesh)
    implicit none
    
    type(mesh_system), intent(inout)    :: mesh
    integer                             :: ie, inod, currset, ierr, nod
    integer, allocatable                :: nodelmat(:,:), counter_nodelmat(:), elmneighbours(:)
    logical                             :: flag=.true.
 
 
    ! Set nummber of sets
    mesh%nset = 6
 
    ! Allocate set variables
    allocate(nodelmat(mesh%nodel,mesh%nnod), stat=ierr)
    allocate(counter_nodelmat(mesh%nnod), stat=ierr)
    allocate(elmneighbours(mesh%nodel*mesh%nodel), stat=ierr)
 
    ! Allocate mesh variables
    allocate(mesh%element_set(mesh%nelm), stat=ierr)
 
    ! Create nodelmat
    nodelmat         = 0
    counter_nodelmat = 1
    do ie=1,mesh%nelm
       do inod=1,mesh%nodel
          nod = mesh%enod(inod,ie)
          nodelmat(counter_nodelmat(nod),nod) = ie
          counter_nodelmat(nod) = counter_nodelmat(nod) + 1
       enddo
    enddo
    
    do while (flag)
       allocate(mesh%set_groups(mesh%nelm,mesh%nset), stat=ierr)
       allocate(mesh%counter_set(mesh%nset), stat=ierr)
       mesh%set_groups  = 0
       mesh%counter_set = 1
       mesh%element_set = 0
       do ie=1,mesh%nelm
          elmneighbours = reshape(nodelmat(:,mesh%enod(:,ie)),[size(elmneighbours)])
          do currset=1,mesh%nset
             if (all(mesh%element_set(pack(elmneighbours,elmneighbours.ne.0)).ne.currset)) then
                mesh%element_set(ie)                               = currset
                mesh%set_groups(mesh%counter_set(currset),currset) = ie
                mesh%counter_set(currset)                          = mesh%counter_set(currset) + 1
                exit
             end if            
          enddo
          if (mesh%element_set(ie).eq.0) then
             print *, 'error at elm ', ie
             mesh%nset = mesh%nset + 1
             deallocate(mesh%set_groups)
             deallocate(mesh%counter_set)
             exit
          endif
       enddo
 
       if (all(mesh%element_set.ne.0)) then
         flag = .false.
       endif
    enddo
 
    ! Change counter_set such that it holds info on number of elms in each set
    mesh%counter_set = mesh%counter_set - 1
    
    deallocate(nodelmat)
    deallocate(counter_nodelmat)
    deallocate(elmneighbours)
 
    write(*,'(A6,I2)') 'nset: ', mesh%nset


    return
  end subroutine get_elm_sets


  subroutine init_mesh_lvlset_subdomain(lvlsety0, lvlsety1, grain_width, mesh)
    implicit none
  
    ! Intent inout
    type(mesh_system), intent(inout) :: mesh

    ! Intent in
    real(dp)                         :: lvlsety0, lvlsety1, grain_width
  
    ! Subroutine variables
    integer, allocatable             :: elm_indecies(:), nod_indecies(:)
    real(dp), allocatable            :: nod_indecies_ed(:,:), nodlbl(:)
    real(dp),parameter               :: tol = 1d-9
    integer                          :: i, ierr
  
  
    ! 1) Find nodes in subdomain
    mesh%nnod_ls  = count((mesh%coord(2,:).ge.(lvlsety0-tol)) .and. (mesh%coord(2,:).lt.(lvlsety1+tol)) & 
                          .and. (mesh%coord(1,:).lt.(grain_width+tol)))
    allocate(mesh%nods_ls(mesh%nnod_ls), stat=ierr)

    allocate(nod_indecies(mesh%nnod), stat=ierr)
    nod_indecies = [(i, i=1,mesh%nnod)]
    mesh%nods_ls = pack(nod_indecies,((mesh%coord(2,:).ge.(lvlsety0-tol)) .and. (mesh%coord(2,:).lt.(lvlsety1+tol)) & 
                                       .and. (mesh%coord(1,:).lt.(grain_width+tol))))

    ! 2) Find elms in subdomain
    allocate(nodlbl(mesh%nnod), stat=ierr)
    nodlbl = 0d0
    where ((mesh%coord(2,:).ge.(lvlsety0-tol)) .and. (mesh%coord(2,:).lt.(lvlsety1+tol)) & 
           .and. (mesh%coord(1,:).lt.(grain_width+tol))) nodlbl = 1d0
    allocate(nod_indecies_ed(mesh%nodel,mesh%nelm), stat=ierr)
    call extract(nod_indecies_ed, nodlbl, mesh%enod,mesh%dofnod_ls)
    mesh%nelm_ls  = count(sum(nod_indecies_ed,1).eq.mesh%nodel)

    allocate(elm_indecies(mesh%nelm), stat=ierr)
    elm_indecies  = [(i, i=1,mesh%nelm)]
    allocate(mesh%elms_ls(mesh%nelm_ls), stat=ierr)
    mesh%elms_ls = pack(elm_indecies,(sum(nod_indecies_ed,1).eq.mesh%nodel))

    ! 3) Find coord in subdomain
    allocate(mesh%coord_ls(2,mesh%nnod_ls), stat=ierr)
    allocate(mesh%newcoord_ls(2,mesh%nnod_ls), stat=ierr)
    mesh%coord_ls = mesh%coord(:,mesh%nods_ls)
    mesh%newcoord_ls = mesh%coord_ls    
    call coordxtr(mesh%ex,mesh%ey,mesh%coord,mesh%enod)


    ! 4) Find enod in subdomain
    call get_enod_ls(mesh)

    ! 5) Find ex and ey in subdomain
    allocate(mesh%ex_ls(mesh%nrgp,mesh%nelm_ls), stat=ierr)
    allocate(mesh%ey_ls(mesh%nrgp,mesh%nelm_ls), stat=ierr)
    allocate(mesh%newex_ls(mesh%nrgp,mesh%nelm_ls), stat=ierr)
    allocate(mesh%newey_ls(mesh%nrgp,mesh%nelm_ls), stat=ierr)
    mesh%ex_ls = mesh%ex(:,mesh%elms_ls)
    mesh%ey_ls = mesh%ey(:,mesh%elms_ls)
    mesh%newex_ls = mesh%ex_ls
    mesh%newey_ls = mesh%ey_ls

    ! Level set sets
    call get_elm_sets_ls(mesh)

    ! Deallocate
    deallocate(elm_indecies)
    deallocate(nod_indecies)
    deallocate(nod_indecies_ed)
    deallocate(nodlbl)

    return
  end subroutine init_mesh_lvlset_subdomain



  subroutine get_enod_ls(mesh)
    ! Routine for extractin enod labeled 1:size max nod in sliced subdomain with same topology connection as original enod
    implicit none

    ! Intent inout
    type(mesh_system), intent(inout) :: mesh

    ! Subroutine variables
    integer, allocatable             :: unique_nums(:)
    integer                          :: i


    ! Find unique nods in sliced enod
    call find_unique_numbers_2D_integer(mesh%enod(:,mesh%elms_ls), unique_nums)


    ! Sort unique nods from low to high
    call bubble_sort_1D_integer(unique_nums)


    ! Substitute nods with nods from 1:size(sliced enod)
    call substitute_number_2D_integer(mesh%enod_ls, mesh%enod(:,mesh%elms_ls), unique_nums, [(i, i=1,size(unique_nums))])

    ! Deallocate
    deallocate(unique_nums)


    return
  end subroutine get_enod_ls


  subroutine find_unique_numbers_2D_integer(arr, unique_nums)
    ! Routine for finding unique numbers in a 2D integer array
    implicit none

    ! Intent inout
    integer, intent(inout), allocatable :: unique_nums(:)
    
    ! Intent in
    integer, intent(in)                 :: arr(:,:)

    ! Subroutine variables
    integer, allocatable                :: unique_nums_tmp(:)
    integer                             :: i, j, k, counter, num, n_rows, n_cols, ierr, nunique

    ! Rows and cols in arr intent in
    n_rows = size(arr,1)
    n_cols = size(arr,2)

    ! Allocate unique_nums_tmp
    allocate(unique_nums_tmp(n_rows*n_cols),stat=ierr)
    unique_nums_tmp = 0

    ! Find unique numbers in arr and store in unique_nums_tmp
    counter = 0
    do j=1,n_cols
      do i=1,n_rows
        num = arr(i,j)
        do k=1,counter
          if (num .eq. unique_nums_tmp(k)) then
            exit
          endif
        enddo
        if (k .gt. counter) then
          counter = counter + 1
          unique_nums_tmp(counter) = num
        endif
      enddo
    enddo

    ! Extract nnz numbers in unique_nums_tmp and send out
    nunique = count(unique_nums_tmp.gt.0)

    ! Allocate unique_nums and put in nnz unique_nums_tmp
    allocate(unique_nums(nunique),stat=ierr)
    unique_nums = unique_nums_tmp(1:counter)

    ! Deallocate
    deallocate(unique_nums_tmp)

    return
  end subroutine find_unique_numbers_2D_integer



  subroutine bubble_sort_1D_integer(a)
    ! Routine for sorting a 1D integer array from low to high
    implicit none

    ! Intent inout
    integer, intent(inout) :: a(:)

    ! Subroutine variables
    integer                :: n, i, j, tmp

    ! Size
    n = size(a)

    ! Sort array from lowest to highest number
    do i=1,n-1
      do j=1,n-i
        if (a(j).gt.a(j+1)) then
          tmp    = a(j)
          a(j)   = a(j+1)
          a(j+1) = tmp
        endif
      enddo
    enddo

    return
  end subroutine bubble_sort_1D_integer


  subroutine substitute_number_2D_integer(new_arr, old_arr, old_nums, new_nums)
    ! Routine for creating a new array from an old array and substituting values in list old_nums with
    ! values in list new_nums. OBS old_nums and new_nums must be of equal size. OBS new_arr and old_arr
    ! have the same format
    implicit none

    ! Intent inout
    integer, intent(inout), allocatable :: new_arr(:,:)    

    ! Intent in
    integer, intent(in)                 :: old_arr(:,:), old_nums(:), new_nums(:)

    ! Subroutine variables
    integer                             :: n_rows, n_cols, nsubs, old_number, new_number, i, j, k, ierr

    ! Rows and cols in arr intent in
    n_rows = size(old_arr,1)
    n_cols = size(old_arr,2)

    ! Number of values to substitute
    nsubs = size(new_nums)

    ! Obs  old_nums and new_num must have same length
    if (nsubs.ne.size(old_nums)) then
      print *, 'OBS old_nums and new_nums is not the same size!'
      return
    endif

    ! Allocate new_arr and set equal to old_array
    allocate(new_arr(n_rows,n_cols), stat=ierr)
    new_arr = old_arr

    do k=1,nsubs
      old_number = old_nums(k)
      new_number = new_nums(k)
      do j=1,n_cols
        do i=1,n_rows
          if (new_arr(i,j) .eq. old_number) then
            new_arr(i,j) = new_number
          endif
        enddo
      enddo
    enddo

    return
  end subroutine substitute_number_2D_integer

  subroutine get_elm_sets_ls(mesh)
    implicit none
    
    type(mesh_system), intent(inout)    :: mesh
    integer                             :: ie, inod, currset, ierr, nod
    integer, allocatable                :: nodelmat(:,:), counter_nodelmat(:), elmneighbours(:)
    logical                             :: flag=.true.
 
 
    ! Set nummber of sets
    mesh%nset = 6
 
    ! Allocate set variables
    allocate(nodelmat(mesh%nodel,mesh%nnod_ls), stat=ierr)
    allocate(counter_nodelmat(mesh%nnod_ls), stat=ierr)
    allocate(elmneighbours(mesh%nodel*mesh%nodel), stat=ierr)
 
    ! Allocate mesh variables
    allocate(mesh%element_set_ls(mesh%nelm_ls), stat=ierr)
 
    ! Create nodelmat
    nodelmat         = 0
    counter_nodelmat = 1
    do ie=1,mesh%nelm_ls
       do inod=1,mesh%nodel
          nod = mesh%enod_ls(inod,ie)
          nodelmat(counter_nodelmat(nod),nod) = ie
          counter_nodelmat(nod) = counter_nodelmat(nod) + 1
       enddo
    enddo
    
    do while (flag)
       allocate(mesh%set_groups_ls(mesh%nelm,mesh%nset), stat=ierr)
       allocate(mesh%counter_set_ls(mesh%nset), stat=ierr)
       mesh%set_groups_ls  = 0
       mesh%counter_set_ls = 1
       mesh%element_set_ls = 0
       do ie=1,mesh%nelm_ls
          elmneighbours = reshape(nodelmat(:,mesh%enod_ls(:,ie)),[size(elmneighbours)])
          do currset=1,mesh%nset
             if (all(mesh%element_set_ls(pack(elmneighbours,elmneighbours.ne.0)).ne.currset)) then
                mesh%element_set_ls(ie)                                  = currset
                mesh%set_groups_ls(mesh%counter_set_ls(currset),currset) = ie
                mesh%counter_set_ls(currset)                             = mesh%counter_set_ls(currset) + 1
                exit
             end if            
          enddo
          if (mesh%element_set_ls(ie).eq.0) then
             print *, 'error at elm ', ie
             mesh%nset = mesh%nset + 1
             deallocate(mesh%set_groups_ls)
             deallocate(mesh%counter_set_ls)
             exit
          endif
       enddo
 
       if (all(mesh%element_set_ls.ne.0)) then
         flag = .false.
       endif
    enddo
 
    ! Change counter_set such that it holds info on number of elms in each set
    mesh%counter_set_ls = mesh%counter_set_ls - 1
    
    ! Deallocate
    deallocate(nodelmat)
    deallocate(counter_nodelmat)
    deallocate(elmneighbours)
 
    return
   end subroutine get_elm_sets_ls


   subroutine clock_time(time, omp_run)
    ! Start clock with respect to omp run or not
    implicit none
 
    real(dp), intent(inout) :: time
    logical, intent(in)     :: omp_run
    
    if (omp_run) then 
      !$ time = OMP_GET_WTIME() 
    else 
      call cpu_time(time)
    endif 
 
   return
  end subroutine clock_time


end module mesh_module

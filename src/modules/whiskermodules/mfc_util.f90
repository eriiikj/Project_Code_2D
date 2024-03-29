module mfc_util
  ! Last modified
  ! E. Jacobsson 2022-04-28
  !    - Routines for handeling mfc lagrange multiplier

  ! mf_datatypes
  use mf_datatypes
  
  ! somelib
  use sparse_util

  ! Whiskerlib
  use mesh_module  

  implicit none

contains


  subroutine init_MFC(K_hat, K, Kcell, T_map_sparse, T_map_sparse_transpose, g, nMFC, ndof_hat, mesh)
    ! Subroutine for initializing MFC
    implicit none

    ! Intent inout    
    type(sparse), intent(inout)          :: K_hat, K, T_map_sparse, T_map_sparse_transpose
    real(dp), allocatable, intent(inout) :: g(:)
    integer, intent(inout)               :: nMFC, ndof_hat
    integer, allocatable, intent(inout)  :: Kcell(:,:)

    ! Intent in
    type(mesh_system), intent(in)        :: mesh

    ! Subroutine variables
    ! Sparse matrices
    type(sparse)                         :: K_exp, T_exp

    ! BC nods
    integer, allocatable                 :: right_boundary_nodes(:),  left_boundary_nodes(:)

    ! MFC parameters
    integer                              :: nMFC_nods, nMFC_dofs_per_constraint

    ! T
    integer , allocatable                :: Tcell(:,:), Tcell_transpose(:,:)
    real(dp), allocatable                :: Tval(:)

    ! T expanded
    integer , allocatable                :: Tcell_exp(:,:)
    real(dp), allocatable                :: Tval_exp(:)

    ! MFC dofs
    nMFC_dofs_per_constraint = 2

    ! Define Tcell, holding indecies of T. Also define g
    ! call init_MFC_T_ux_uy_equal(Tcell, Tcell_transpose, Tval, nMFC_nods, nMFC, nMFC_dofs_per_constraint, &
    !                             right_boundary_nodes, left_boundary_nodes, mesh)
    call init_MFC_T_ux_right_equal(Tcell, Tcell_transpose, Tval, nMFC, mesh)

    call init_MFC_g(g, nMFC)    

    ! Number of dofs in expanded format
    ndof_hat = mesh%ndof+nMFC

    ! T_map_sparse
    call spacellDef(T_map_sparse, Tcell, nMFC, mesh%ndof)
    call spaputval(T_map_sparse, Tcell, Tval)

    ! T_map_sparse transpose
    call spacellDef(T_map_sparse_transpose, Tcell_transpose, mesh%ndof, nMFC)
    call spaputval(T_map_sparse_transpose, Tcell_transpose, Tval)
    
    ! Expanded T
    call expand_MFC_T(Tcell_exp, Tval_exp, Tcell, Tval, mesh%ndof)
    call spacellDef(T_exp,Tcell_exp,ndof_hat,ndof_hat)
    call spaputval(T_exp, Tcell_exp, Tval_exp)   

    ! ----- Define sparse K -----
    call spaTopDef(K,mesh%enod,mesh%dofnod)

    ! K indecies collected in Kcell
    call get_Kcell(Kcell,K)

    ! Define K_exp
    call spaCellDef(K_exp, Kcell, ndof_hat, ndof_hat)

    ! Define K_hat where both K and T indecies are allocated
    call spaadd(K_hat,K_exp,T_exp)

    ! Deallocate
    ! deallocate(right_boundary_nodes)
    ! deallocate(left_boundary_nodes)
    deallocate(Tcell)
    deallocate(Tcell_transpose)
    deallocate(Tval)
    deallocate(Tcell_exp)
    deallocate(Tval_exp)

    ! Remove sparse matrices
    call sparemove(K_exp)
    call sparemove(T_exp)
    
    return
  end subroutine init_MFC  
    



  subroutine init_MFC_T_ux_uy_equal(Tcell, Tcell_transpose, Tval, nMFC_nods, nMFC, nMFC_dofs_per_constraint, &
                                    right_boundary_nodes, left_boundary_nodes, mesh)
    ! Routine for defining mfc such that ux and uy displacement equal on both sides
    implicit none

    ! Intent inout
    integer , allocatable, intent(inout) :: Tcell(:,:), Tcell_transpose(:,:), right_boundary_nodes(:),  left_boundary_nodes(:)
    real(dp), allocatable, intent(inout) :: Tval(:)
    integer, intent(inout)               :: nMFC_nods, nMFC

    ! Intent in
    integer , intent(in)                 :: nMFC_dofs_per_constraint
    type(mesh_system), intent(in)        :: mesh

    integer                              :: nTcell
    integer                              :: counter, col, row, ierr, nnod, i, mloc

    integer , allocatable                :: node_indecies(:)
    integer , allocatable                :: left_boundary_dofs_x(:) , left_boundary_dofs_y(:) , left_boundary_dofs(:), & 
                                            right_boundary_dofs_x(:), right_boundary_dofs_y(:), right_boundary_dofs(:)

    integer, allocatable                 :: Tdof(:,:)
    

    real(dp), allocatable                :: node_coord(:,:)
    logical, allocatable                 :: mask_vec(:)


    ! Topology
    nnod      = size(mesh%coord,2)

    ! Allocate node indecies
    allocate(node_indecies(nnod))
    node_indecies = [(i, i=1,nnod)]

    print *, 'model_width: ', mesh%model_width
    ! number of MFC
    ! nMFC_nods   = size(pack(node_indecies,mesh%coord(1,:).gt.(mesh%elmsize_x*(mesh%nelmx-0.01))))
    nMFC_nods   = size(pack(node_indecies,mesh%coord(1,:).gt.(mesh%model_width-1d-8)))
    nMFC        = nMFC_nods*mesh%dofnod

    ! print *, 'nMFC_nods',nMFC_nods
    write(*,'(A6,I3)') 'nMFC: ',nMFC
  
    ! Allocate left and right side nodes
    allocate(left_boundary_nodes(nMFC_nods) , stat=ierr)
    allocate(right_boundary_nodes(nMFC_nods), stat=ierr)

    ! Allocate right side dofs
    allocate(right_boundary_dofs_x(nMFC_nods), stat=ierr)
    allocate(right_boundary_dofs_y(nMFC_nods), stat=ierr)
    allocate(right_boundary_dofs(nMFC)       , stat=ierr)

    ! Allocate left side dofs
    allocate(left_boundary_dofs_x(nMFC_nods), stat=ierr)
    allocate(left_boundary_dofs_y(nMFC_nods), stat=ierr)
    allocate(left_boundary_dofs(nMFC)       , stat=ierr)

    ! --- Unsorted boundary nodes ---
    right_boundary_nodes = pack(node_indecies,mesh%coord(1,:).gt.(mesh%elmsize_x*(mesh%nelmx-0.5d0)))
    left_boundary_nodes  = pack(node_indecies,mesh%coord(1,:).lt.mesh%elmsize_x*0.5d0)
    
    ! --- Sort nodes from bottom to top ---
    allocate(mask_vec(nMFC_nods), stat=ierr)
    allocate(node_coord(nMFC_nods,2))

    ! Sort right boundary nodes
    mask_vec = .true.
    node_coord(:,1) = right_boundary_nodes          ! nodes
    node_coord(:,2) = mesh%coord(2,right_boundary_nodes) ! ycoords
    do i=1,size(node_coord,1)
        mloc = minloc(node_coord(:,2),1,mask_vec)
        right_boundary_nodes(i) = nint(node_coord(mloc,1))
        mask_vec(mloc) = .FALSE.
    enddo

    ! Sort left boundary nodes
    mask_vec = .true.
    node_coord(:,1) = left_boundary_nodes          ! nodes
    node_coord(:,2) = mesh%coord(2,left_boundary_nodes) ! ycoords
    do i=1,size(node_coord,1)
        mloc = minloc(node_coord(:,2),1,mask_vec)
        left_boundary_nodes(i) = nint(node_coord(mloc,1))
        mask_vec(mloc) = .FALSE.
    enddo

    ! --- Boundary dofs ---
    ! right side
    right_boundary_dofs_x = (right_boundary_nodes-1)*2+1
    right_boundary_dofs_y = (right_boundary_nodes-1)*2+2
    right_boundary_dofs(1:nMFC:2) = right_boundary_dofs_x
    right_boundary_dofs(2:nMFC:2) = right_boundary_dofs_y

    ! left side
    left_boundary_dofs_x = (left_boundary_nodes-1)*2+1
    left_boundary_dofs_y = (left_boundary_nodes-1)*2+2
    left_boundary_dofs(1:nMFC:2) = left_boundary_dofs_x
    left_boundary_dofs(2:nMFC:2) = left_boundary_dofs_y

    ! Tdof - topology of MFC collected. Each MFC dofs stored columnwise.
    allocate(Tdof(nMFC_dofs_per_constraint,nMFC))
    Tdof(1,:) = right_boundary_dofs
    Tdof(2,:) = left_boundary_dofs

    ! Tcell
    nTcell = nMFC*nMFC_dofs_per_constraint
    allocate(Tcell(nTcell,2))
    counter = 1
    do col=1,nMFC
        do row=1,nMFC_dofs_per_constraint
          Tcell(counter,:) = [col,Tdof(row,col)]
          counter = counter + 1
        enddo
    enddo 

    allocate(Tcell_transpose(nTcell,2))
    Tcell_transpose(:,1) = Tcell(:,2)
    Tcell_transpose(:,2) = Tcell(:,1)

    ! Tval - constraint values. Defined for each element in Tcell
    allocate(Tval(nMFC*2))
    Tval(1:nMFC*2:2) =  1d0
    Tval(2:nMFC*2:2) = -1d0

    ! Deallocate
    deallocate(Tdof)
    deallocate(node_indecies)
    deallocate(left_boundary_dofs_x)
    deallocate(left_boundary_dofs_y)
    deallocate(left_boundary_dofs)
    deallocate(right_boundary_dofs_x)
    deallocate(right_boundary_dofs_y)
    deallocate(right_boundary_dofs)
    deallocate(node_coord)
    deallocate(mask_vec)

    if (size(right_boundary_nodes).ne.nMFC_nods) then
      print *, 'Size right_boundary_nodes not correct'
      call exit()
    endif

    return
  end subroutine init_MFC_T_ux_uy_equal



  subroutine init_MFC_T_ux_right_equal(Tcell, Tcell_transpose, Tval, nMFC, mesh)
    ! Routine for defining mfc such that ux and uy displacement equal on both sides
    implicit none

    ! Intent inout
    integer , allocatable, intent(inout) :: Tcell(:,:), Tcell_transpose(:,:)
    real(dp), allocatable, intent(inout) :: Tval(:)
    integer, intent(inout)               :: nMFC

    ! Intent in
    type(mesh_system), intent(in)        :: mesh

    integer                              :: nTcell
    integer                              :: counter, col, row, ierr, nnod, nMFC_dofs_per_constraint
    integer , allocatable                :: right_boundary_dofs_x(:)

    integer, allocatable                 :: Tdof(:,:)


    ! Topology
    nMFC_dofs_per_constraint = 2
    nnod                     = size(mesh%coord,2)

    ! Number of MFC
    nMFC = size(mesh%bcnods_right_side) - 1
    write(*,'(A6,I3)') 'nMFC: ',nMFC

    ! Allocate right side dofs
    allocate(right_boundary_dofs_x(size(mesh%bcnods_right_side)), stat=ierr)

    ! --- Boundary dofs ---
    right_boundary_dofs_x = (mesh%bcnods_right_side-1)*2+1

    ! Tdof - topology of MFC collected. Each MFC dofs stored columnwise
    allocate(Tdof(2,nMFC))
    Tdof(1,:) = right_boundary_dofs_x(2:)
    Tdof(2,:) = right_boundary_dofs_x(1)

    ! Tcell
    nTcell = nMFC*nMFC_dofs_per_constraint
    allocate(Tcell(nTcell,2))
    counter = 1
    do col=1,nMFC
      do row=1,nMFC_dofs_per_constraint
        Tcell(counter,:) = [col,Tdof(row,col)]
        counter = counter + 1
      enddo
    enddo

    allocate(Tcell_transpose(nTcell,2))
    Tcell_transpose(:,1) = Tcell(:,2)
    Tcell_transpose(:,2) = Tcell(:,1)

    ! Tval - constraint values. Defined for each element in Tcell
    allocate(Tval(nMFC*2))
    Tval(1:nMFC*2:2) =  1d0
    Tval(2:nMFC*2:2) = -1d0

    ! Deallocate
    deallocate(Tdof)
    deallocate(right_boundary_dofs_x)

    return
  end subroutine init_MFC_T_ux_right_equal

  subroutine init_MFC_g(g, nMFC)
    ! Compute vector g used in MFC
    implicit none

    ! Intent inout
    real(dp), allocatable, intent(inout) :: g(:)

    ! Intent in
    integer, intent(in)                  :: nMFC

    ! Subroutine variables    
    integer                              :: ierr    
  
    ! Define MFC right hand side g, compare Eq 9.22 Felippa
    allocate(g(nMFc), stat=ierr)
    g = 0d0

    return
  end subroutine init_MFC_g


  subroutine expand_MFC_T(Tcell_exp, Tval_exp, Tcell, Tval, ndof)
    implicit none

    integer , intent(in)                 :: Tcell(:,:), ndof
    real(dp), intent(in)                 :: Tval(:)

    integer                              :: nTcell, nTcell_exp

    integer , allocatable, intent(inout) :: Tcell_exp(:,:)
    real(dp), allocatable, intent(inout) :: Tval_exp(:)

    ! Expanded Tcell
    nTcell = size(Tcell,1)
    nTcell_exp = nTcell*2
    allocate(Tcell_exp(nTcell_exp,2))
    Tcell_exp = reshape([[Tcell(:,1)+ndof,Tcell(:,2)],[Tcell(:,2),Tcell(:,1)+ndof]],[nTcell_exp,2])

    ! Expaned Tval
    allocate(Tval_exp(nTcell_exp))
    Tval_exp = [Tval,Tval]

    return
  end subroutine expand_MFC_T

  subroutine get_Kcell(Kcell, K)
    implicit none
    type(sparse), intent(in)            :: K

    integer                             :: n_K_rows, n_K_elms, colstaken, counter, row, col, elmsrowi

    integer, allocatable, intent(inout) :: Kcell(:,:)

    ! K indecies collected in Kcell
    n_K_rows  = size(K%ia)-1
    n_K_elms  = size(K%a)
    allocate(Kcell(n_K_elms,2))

    colstaken = 0
    counter = 1
    do row=1,n_K_rows
      elmsrowi = K%ia(row+1)-K%ia(row)
      do col=colstaken+1,colstaken + elmsrowi       
        Kcell(counter,:) = [row,K%ja(col)]
        counter = counter+1
      enddo
      colstaken = colstaken + elmsrowi
    enddo

    return
  end subroutine get_Kcell


end module mfc_util

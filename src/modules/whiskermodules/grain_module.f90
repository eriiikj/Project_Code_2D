module grain_module
    ! Last modified E.Jacobsson
    !  2021-11-30
    !    - Implemented print_nr_elms, breaking if not matching
    !
    !    2021-12-15
    !    - Grain initiation in a seperate subroutine 'init_grains'
    !
    !   2022-03-01
    !   - IMC boundary defined per Gauss point
    !  2022-03-28
    !  - replacing real operators (+,-,eq,ne) for logicals with logical operators (.or., .neqv.)
    
    ! General info:
    ! Model size: tot width = 12 microns, tot height = 4 microns, Cu_height = 1 micron
    ! 4 grains with equal width (3 microns)
    
    ! Intrinsic modules
    use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
 
    ! mf_datatypes
    use mf_datatypes
 
    ! somelib modules
    use fem_util
    use mater_large
    use matrix_util
    
    ! Whiskerlib
    use gp_util
    use mesh_module
        
    implicit none
 
    type grain_system
      integer               :: ngrains
      real(dp)              :: G1_xmid, G2_xmid
      logical,  allocatable :: G1_gps(:,:), G2_gps(:,:)
      real(dp)              :: grain_width
      logical, allocatable  :: Cu_gps(:,:), IMC_gps(:,:), IMCSn_gps(:,:), Sn_gps(:,:)
      logical, allocatable  :: Cu_elms(:), IMC_elms(:), IMCSn_elms(:), Sn_elms(:), mater_vec(:)
      integer, allocatable  :: ieloc_Cu(:), gploc_Cu(:), ieloc_IMC(:), gploc_IMC(:), ieloc_IMCSn(:), gploc_IMCSn(:)
      integer, allocatable  :: ieloc_Sn(:), gploc_Sn(:)
      integer               :: nloop_Cu, nloop_IMC, nloop_IMCSn, nloop_Sn
      integer               :: nIMCGPS
    end type grain_system
 
contains

   subroutine init_grains(grains, coord, enod, nrgp, gpx, gpy, Cu_height)
      implicit none
 
      ! Intent inout
      type(grain_system), intent(inout) :: grains

      ! Intent in
      real(dp), intent(in)              :: coord(:,:), gpx(:,:), gpy(:,:), Cu_height
      integer,  intent(in)              :: enod(:,:), nrgp

      ! Subroutine variables
      real(dp)                          :: tot_width
      integer                           :: ierr, nelm
   
      ! Get nelm
      nelm = size(enod,2)
      
      ! Allocate      
   
      allocate(grains%G1_gps(nrgp,nelm)   , stat=ierr)
      allocate(grains%G2_gps(nrgp,nelm)   , stat=ierr)
   
      allocate(grains%Cu_gps(nrgp,nelm)   , stat=ierr)
      allocate(grains%IMC_gps(nrgp,nelm)  , stat=ierr)
      allocate(grains%IMCSn_gps(nrgp,nelm), stat=ierr)
      allocate(grains%Sn_gps(nrgp,nelm)   , stat=ierr)
      allocate(grains%Cu_elms(nelm)       , stat=ierr)
      allocate(grains%IMC_elms(nelm)      , stat=ierr)
      allocate(grains%IMCSn_elms(nelm)    , stat=ierr)
      allocate(grains%Sn_elms(nelm)       , stat=ierr)
   
      allocate(grains%ieloc_Cu(nelm*nrgp)   , stat=ierr)
      allocate(grains%gploc_Cu(nelm*nrgp)   , stat=ierr)
      allocate(grains%ieloc_IMC(nelm*nrgp)  , stat=ierr)
      allocate(grains%gploc_IMC(nelm*nrgp)  , stat=ierr)
      allocate(grains%ieloc_IMCSn(nelm*nrgp), stat=ierr)
      allocate(grains%gploc_IMCSn(nelm*nrgp), stat=ierr)
      allocate(grains%ieloc_Sn(nelm*nrgp)   , stat=ierr)
      allocate(grains%gploc_Sn(nelm*nrgp)   , stat=ierr)
      allocate(grains%mater_vec(nelm*nrgp)  , stat=ierr)
      

      ! Specify number of grains
      grains%ngrains = 1

      ! Obtain geometry input
      tot_width          = maxval(coord(1,:)) - minval(coord(1,:))
      grains%grain_width = tot_width/grains%ngrains

      ! Determine Cu gps and Cu elms
      grains%Cu_gps  = .false.
      grains%Cu_elms = .false.
      grains%Cu_gps  = abs(gpy).lt.abs(Cu_height)
      grains%Cu_elms = any(grains%Cu_gps,1)
      
      ! Find gps beloning to individual grains
      grains%G1_gps = .false.
      ! grains%G2_gps = .false.      
      where (gpx.lt.grains%grain_width)                                     grains%G1_gps = .true.
      ! where (gpx.gt.grains%grain_width .and. gpx.lt.grains%grain_width*2d0) grains%G2_gps = .true.

      ! x positions of rotation centers of grains; center of grains
      grains%G1_xmid = grains%grain_width*0.5d0
      ! grains%G2_xmid = grains%grain_width*1.5d0

      ! Init nIMCGPs
      grains%nIMCGPS = 0
      
      return
   end subroutine init_grains
 
 
 subroutine update_IMC_gps(grains, phi_gp)
    implicit none

    ! Intent inout
    type(grain_system), intent(inout) :: grains
    
    ! Intent in
    real(dp), intent(inout)           :: phi_gp(:,:)

    ! Subroutine variables
    integer                           :: i

    ! --- IMC and Sn gps ---
    grains%IMC_gps   = .false.
    grains%IMCSn_gps = .false.
    grains%Sn_gps    = .false.
    where (phi_gp.gt.0.99d0    .and. .not. grains%Cu_gps)                               grains%IMC_gps      = .true.
    where (phi_gp.lt.0.99d0    .and. phi_gp.gt.0.01d0   .and. .not. grains%Cu_gps)      grains%IMCSn_gps    = .true.
    where (.not. grains%Cu_gps .and. .not. grains%IMC_gps .and. .not. grains%IMCSn_gps) grains%Sn_gps       = .true.
 
    ! IMC and Sn elms
    grains%Cu_elms    = any(grains%Cu_gps,1)
    grains%IMC_elms   = any(grains%IMC_gps,1)
    grains%IMCSn_elms = any(grains%IMCSn_gps,1)
    grains%Sn_elms    = any(grains%Sn_gps,1)

    ! --- Update indecies ---


    ! Update Cu indecies
    grains%mater_vec                                     = reshape(grains%Cu_gps, (/size(grains%mater_vec)/))
    grains%nloop_Cu                                      = count(grains%mater_vec)
    grains%gploc_Cu(1:grains%nloop_Cu)                   = pack([(i,i=1,size(grains%mater_vec))],grains%mater_vec.eqv..true.)
    where (grains%gploc_Cu .ne. 0) grains%ieloc_Cu       = (grains%gploc_Cu-1)/4 + 1
    where (grains%gploc_Cu .ne. 0) grains%gploc_Cu       = grains%gploc_Cu - (grains%ieloc_Cu-1)*4
 
    ! Update IMC indecies
    grains%mater_vec                                     = reshape(grains%IMC_gps, (/size(grains%mater_vec)/))
    grains%nloop_IMC                                     = count(grains%mater_vec)
    grains%gploc_IMC(1:grains%nloop_IMC)                 = pack([(i,i=1,size(grains%mater_vec))],grains%mater_vec.eqv..true.)
    where (grains%gploc_IMC .ne. 0) grains%ieloc_IMC     = (grains%gploc_IMC-1)/4 + 1
    where (grains%gploc_IMC .ne. 0) grains%gploc_IMC     = grains%gploc_IMC - (grains%ieloc_IMC-1)*4

    ! Update IMCSn indecies
    grains%mater_vec                                     = reshape(grains%IMCSn_gps, (/size(grains%mater_vec)/))
    grains%nloop_IMCSn                                   = count(grains%mater_vec)
    grains%gploc_IMCSn(1:grains%nloop_IMCSn)             = pack([(i,i=1,size(grains%mater_vec))],grains%mater_vec.eqv..true.)
    where (grains%gploc_IMCSn .ne. 0) grains%ieloc_IMCSn = (grains%gploc_IMCSn-1)/4 + 1
    where (grains%gploc_IMCSn .ne. 0) grains%gploc_IMCSn = grains%gploc_IMCSn - (grains%ieloc_IMCSn-1)*4
 
    ! Update Sn indecies
    grains%mater_vec                                     = reshape(grains%Sn_gps, (/size(grains%mater_vec)/))
    grains%nloop_Sn                                      = count(grains%mater_vec)
    grains%gploc_Sn(1:grains%nloop_Sn)                   = pack([(i,i=1,size(grains%mater_vec))],grains%mater_vec.eqv..true.)
    where (grains%gploc_Sn .ne. 0) grains%ieloc_Sn       = (grains%gploc_Sn-1)/4 + 1
    where (grains%gploc_Sn .ne. 0) grains%gploc_Sn       = grains%gploc_Sn - (grains%ieloc_Sn-1)*4
 
     return
   end subroutine update_IMC_gps
 
 
   subroutine print_nr_elms(grains, nelm)
     implicit none
     type(grain_system), intent(in) :: grains
     integer                        :: nelm
 
 
     print *, 'Elements: '
     write(*,*)'Number of Cu elms        :', count(grains%Cu_elms)
     write(*,*)'Number of IMC elms       :', count(grains%IMC_elms)
     write(*,*)'Number of IMC/Sn elms    :', count(grains%IMCSn_elms)
     write(*,*)'Number of Sn elms        :', count(grains%Sn_elms)
     write(*,*)'Number of elms           :', nelm
     
     return
   end subroutine print_nr_elms
 
 
   subroutine print_nr_gps(grains)
     implicit none

     type(grain_system), intent(inout) :: grains
     integer                           :: nCu_gps, nIMC_gps, nIMCSn_gps, nSn_gps, ngps, nIMC_gps_diff

     nCu_gps    = count(grains%Cu_gps)
     nIMC_gps   = count(grains%IMC_gps)
     nIMCSn_gps = count(grains%IMCSn_gps)
     nSn_gps    = count(grains%Sn_gps)
     ngps       = nCu_gps + nIMC_gps + nIMCSn_gps + nSn_gps
     
     nIMC_gps_diff = nIMC_gps + nIMCSn_gps - grains%nIMCGPS

     grains%nIMCGPS = nIMC_gps + nIMCSn_gps
     
     write(*,"(A13)") 'Gauss points:'
     write(*,"(A24,I8)") 'Number of Cu gps     : ', nCu_gps
     write(*,"(A24,I8)") 'Number of IMC gps    : ', nIMC_gps
     write(*,"(A24,I8)") 'Number of IMC/Sn gps : ', nIMCSn_gps
     write(*,"(A24,I8)") 'Number of Sn gps     : ', nSn_gps
     write(*,"(A24,I8)") 'Number of gps        : ', ngps
     write(*,"(A24,I8)") 'Number of new IMC gps: ', nIMC_gps_diff
 
     ! Check so that each gp is assigned precisely one material
    !  if (any((grains%Cu_gps.neqv.grains%old_IMC_gps.neqv.grains%new_IMC_gps.neqv.grains%Sn_gps).neqv..true.)) then
    !     !pause 'More than one material at one or more gps'
    !  endif
    !  if (ngps.ne.(nelm*nrgp)) then
    !     !pause 'All gps not assigned'
    !  endif
     
     return
   end subroutine print_nr_gps
 
 
   subroutine global_gpel_check(global_check, grains, nelm)
    implicit none
 
    ! Intent in
    real(dp), intent(in)                  :: global_check(:,:)
    type(grain_system), intent(in)        :: grains
    
    ! Subroutine variables
    integer                               :: ie, igp, nelm
 

    do ie=1,nelm
       do igp=1,4
          if (global_check(igp,ie).eq.1d0) then
             !print *, 'ie: ', ie
             !print *, 'global check accepted: ', global_check(:,ie)
          else
             print *, 'global error at: ', igp,ie
             print *, 'global_check(:,ie): ', global_check(:,ie)
             if (grains%Cu_elms(ie)) then
                !pause 'error in Cu'
             endif
             if (grains%IMC_elms(ie)) then
                !pause 'error in old_IMC'
             endif
             if (grains%Sn_elms(ie)) then
                !pause 'error in Sn'
             endif
          end if
       enddo
    enddo
 
    return
   end subroutine global_gpel_check
     
 
 end module grain_module
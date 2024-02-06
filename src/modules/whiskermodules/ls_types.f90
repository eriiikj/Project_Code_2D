module ls_types
    ! Last modified 
    ! E. Jacobsson 2023-09-04
      
    ! mf_datatypes
    use mf_datatypes

    ! Somelib
    use fem_system

    implicit none

    type ls_system
      integer               :: ngrains, nseg_alloc
      real(dp), allocatable :: a(:,:), ed(:,:,:), a_gp(:,:,:),jnod(:,:), bcval(:,:), a_prev(:,:)
      real(dp),allocatable  :: line_ex(:,:), line_ey(:,:), jed(:,:,:), tp_points(:,:)
      integer, allocatable  :: line_seg(:), line_elms(:,:), closest_line(:,:)
      logical, allocatable  :: int_elms(:,:), line_seg_cc(:), tplines(:,:)
      real(dp), allocatable :: xvec(:), yvec(:), rvec(:), EEvec(:)
      real(dp), allocatable :: vp(:,:,:), vpm(:,:)
      integer, allocatable  :: rm_lines(:,:), rm_lines_tmp(:,:), sep_lines(:,:), nsep_lines(:)
      real(dp)              :: alpha=1d4, IMC_area, IMC_area_init

      
      ! Update level set function
      real(dp)              :: theta = 1d0, gamma, m, zeta, mzeta, mzetagb, malpha
      real(dp)              :: vk, D(2,2)
      real(dp)              :: h                                                 ! Time step
      real(dp)              :: w                                                 ! Interpolation width
      real(dp)              :: time = 0d0
      type(sparse)          :: C_hat, K_hat

      ! Stress loop
      real(dp),allocatable  :: phi(:,:), phi_ed(:,:,:), phi_gp(:,:,:)
      real(dp),allocatable  :: hphi(:,:), hphi_ed(:,:,:), hphi_gp(:,:,:), hphi_sum(:)      
      integer, allocatable  :: material(:)
      logical, allocatable  :: old_IMC_gp_flag(:,:)

      ! Sn plot grain
      real(dp), allocatable :: sn_a(:), sn_ed(:,:), sn_a_gp(:,:)
      real(dp),allocatable  :: sn_hphi_plot(:), sn_hphi_ed_plot(:,:), sn_hphi_gp_plot(:,:), sn_hphi_gp_plot2(:,:) 
      
    end type ls_system


end module ls_types
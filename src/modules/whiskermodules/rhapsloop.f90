module rhapsloop
   
   ! OpenMP
   use omp_lib

   ! Mflib
   use mf_datatypes

   ! Somelib
   use sparse_util
   use fem_util
   use fem_system
   use matrix_util
   use mater_large   
   use mater_J2iso_Cu
   use mater_hyperel_IMC
   use mater_J2iso_Sn
   use elem_large_cont_2d

   ! Whiskerlib
   use mesh_module
   use mfc_util   
   use imc_transf
   use output_data
   use ls_types   
   

   implicit none

   ! Define private variables

   ! K_hat
   type(sparse)          :: K_hat
   type(sparse)          :: K
   integer, allocatable  :: Kcell(:,:)

   ! T
   type(sparse)          :: T_map_sparse, T_map_sparse_transpose
   real(dp), allocatable :: g(:), lambda(:), lambdaold(:)
   real(dp), allocatable :: res_hat(:), a_hat(:), a_hatold(:)
   integer               :: nMFC, ndof_hat  

   ! Normal a and res (and others)
   real(dp), allocatable :: ed(:,:), es(:,:,:), es_g(:,:,:,:), dg(:,:,:), dge(:,:,:), fint(:)
   real(dp), allocatable :: f_ext(:), a(:), res(:), cau(:,:,:), elm_label(:), global_check(:,:)

   ! Old state
   real(dp), allocatable :: edold(:,:), dgold(:,:,:), dgeold(:,:,:), aold(:)   

   ! Material parameters
   real(dp)              :: ep, E_Cu, v_Cu, mp_Cu(6), E_IMC, v_IMC, mp_IMC(6), E_Sn, v_Sn, mp_Sn(6)
   real(dp), allocatable :: Dgp(:,:,:,:), Dgp_g(:,:,:,:,:), E_save_gp(:,:), yo_save(:,:)

   ! Phi_gp_old and phi_gp_avg
   real(dp), allocatable :: phi_gp_old(:,:), phi_gp_avg(:,:,:)
   
   ! Time
   real(dp)              :: stiffness_time, assem_stiffness_time, assem_force_time, spaputval_time
   real(dp)              :: solveq_time, force_time
   real(dp)              :: computeplot_time, mat_time, plotvtk_time

   ! Private
   private K_hat, K, Kcell, T_map_sparse, T_map_sparse_transpose, g, lambda, res_hat, a_hat, nMFC, ndof_hat
   private ed, es, dg, dge, fint, f_ext, a, res, cau, elm_label, global_check
   private ep, E_Cu, v_Cu, mp_Cu, E_IMC, v_IMC, mp_IMC, E_Sn, v_Sn, mp_Sn, Dgp
   private phi_gp_old, phi_gp_avg
   private stiffness_time, assem_stiffness_time, assem_force_time, spaputval_time, solveq_time, force_time
   private computeplot_time, mat_time, plotvtk_time

contains

   subroutine init_NewtonEqIter(mesh, lssys)
      implicit none
      
      ! Intent inout
      type(mesh_system), intent(inout)         :: mesh

      ! Intent in
      type(ls_system), intent(in)              :: lssys

      ! Subroutine variables
      integer                                  :: ierr      
      
      ! Allocate private variables for Newton iteration
      allocate(ed(8,mesh%nelm), stat=ierr)      
      allocate(es(4,mesh%nrgp,mesh%nelm)        , stat=ierr)
      allocate(es_g(4,mesh%nrgp,mesh%nelm,lssys%ngrains), stat=ierr)
      allocate(dg(4,mesh%nrgp,mesh%nelm)        , stat=ierr)
      allocate(dge(9,mesh%nrgp,mesh%nelm)       , stat=ierr)
      allocate(fint(mesh%ndof)                  , stat=ierr)
      allocate(f_ext(mesh%ndof)                 , stat=ierr)
      allocate(a(mesh%ndof)                     , stat=ierr)
      allocate(res(mesh%ndof)                   , stat=ierr)
      allocate(mesh%bcdof(size(mesh%bcval))     , stat=ierr)
      allocate(Dgp(3,3,mesh%nrgp,mesh%nelm)     , stat=ierr)
      allocate(Dgp_g(3,3,mesh%nrgp,mesh%nelm,lssys%ngrains), stat=ierr)
      allocate(yo_save(mesh%nrgp,mesh%nelm)     , stat=ierr)      
      allocate(elm_label(mesh%nelm)             , stat=ierr)
      allocate(global_check(mesh%nrgp,mesh%nelm), stat=ierr)
      allocate(phi_gp_old(mesh%nrgp,mesh%nelm)  , stat=ierr)
      allocate(phi_gp_avg(mesh%nrgp,mesh%nelm,4)  , stat=ierr)

      ! Old state
      allocate(edold(8,mesh%nelm), stat=ierr)
      allocate(dgold(4,mesh%nrgp,mesh%nelm)        , stat=ierr)
      allocate(dgeold(9,mesh%nrgp,mesh%nelm)       , stat=ierr)
      allocate(aold(mesh%ndof)                     , stat=ierr)

      ! Init phi_gp_old
      phi_gp_old = 0d0

      ! Init MFC
      call init_MFC(K_hat, K, Kcell, T_map_sparse, T_map_sparse_transpose, g, nMFC, ndof_hat, mesh)
   
      ! Allocate expanded form
      allocate(lambda(nMFC)     , stat=ierr)
      allocate(res_hat(ndof_hat), stat=ierr)
      allocate(a_hat(ndof_hat)  , stat=ierr)

      allocate(lambdaold(nMFC)     , stat=ierr)
      allocate(a_hatold(ndof_hat)  , stat=ierr)


      ! -- Define material properties in mm, taken from 'Fundamentals of materials science and engineering' --

      ! Thickness
      ep = 1d0

      ! --- JH ---
      
      ! ! Cu, von Mises
      ! E_Cu     = 150d3
      ! v_Cu     = 0.35d0
      ! mp_Cu(1) = E_Cu/3d0/(1d0-2d0*v_Cu)     ! K_Cu
      ! mp_Cu(2) = E_Cu/2d0/(1d0+v_Cu)         ! G_Cu
      ! mp_Cu(3) = 170d0                       ! Initial yield stress
      ! mp_Cu(4) = 0.69d0                      ! Hardening modulus
      ! mp_Cu(5) = 0d0                         ! Hardening exponent
      ! mp_Cu(6) = 0d0                         ! Saturation stress
      ! call J2iso_Cu_init(mp_Cu, mesh%nelm, mesh%nrgp)

      ! ! IMC, von Mises
      ! E_IMC     = 112.3d3
      ! v_IMC     = 0.31d0
      ! mp_IMC(1) = E_IMC/3d0/(1d0-2d0*v_IMC)  ! K_Cu
      ! mp_IMC(2) = E_IMC/2d0/(1d0+v_IMC)      ! G_Cu

      ! ! Init Neo-Hooke for IMC material
      ! call neohooke_init(mp_IMC)

      ! ! Sn, von Mises
      ! E_Sn     = 19d3
      ! v_Sn     = 0.36d0
      ! mp_Sn(1) = E_Sn/3d0/(1d0-2d0*v_Sn)     ! K_Sn
      ! mp_Sn(2) = E_Sn/2d0/(1d0+v_Sn)         ! G_Sn
      ! mp_Sn(3) = 12d0                        ! Initial yield stress
      ! mp_Sn(4) = 0.75d0                      ! Hardening modulus
      ! mp_Sn(5) = 0d0                         ! Hardening exponent
      ! mp_Sn(6) = 0d0                         ! Saturation stress
      ! call J2iso_Sn_init(mp_Sn, mesh%nelm, mesh%nrgp)
      
      

      ! --- Chason ---

      ! Cu, von Mises
      E_Cu     = 117d3
      v_Cu     = 0.34d0
      mp_Cu(1) = E_Cu/3d0/(1d0-2d0*v_Cu)     ! K_Cu
      mp_Cu(2) = E_Cu/2d0/(1d0+v_Cu)         ! G_Cu
      mp_Cu(3) = 200d0                       ! Initial yield stress
      mp_Cu(4) = 0d0                         ! Hardening modulus
      mp_Cu(5) = 0d0                         ! Hardening exponent
      mp_Cu(6) = 0d0                         ! Saturation stress
      call J2iso_Cu_init(mp_Cu, mesh%nelm, mesh%nrgp)

      ! IMC, Neo-Hookean
      E_IMC     = 86d3
      v_IMC     = 0.30d0
      mp_IMC(1) = E_IMC/3d0/(1d0-2d0*v_IMC)  ! K_Cu
      mp_IMC(2) = E_IMC/2d0/(1d0+v_IMC)      ! G_Cu
      call neohooke_init(mp_IMC)

      ! Sn, von Mises
      E_Sn     = 50d3 !19d3
      v_Sn     = 0.36d0
      mp_Sn(1) = E_Sn/3d0/(1d0-2d0*v_Sn)     ! K_Sn
      mp_Sn(2) = E_Sn/2d0/(1d0+v_Sn)         ! G_Sn
      mp_Sn(3) = 14.5d0                      ! Initial yield stress
      mp_Sn(4) = 0d0                         ! Hardening modulus
      mp_Sn(5) = 0d0                         ! Hardening exponent
      mp_Sn(6) = 0d0                         ! Saturation stress
      call J2iso_Sn_init(mp_Sn, mesh%nelm, mesh%nrgp)


      ! E save
      E_save_gp = 0d0

      ! Initiate data
      es                   = 0d0
      ed                   = 0d0
      dg                   = 0d0
      dg(1,:,:)            = 1d0      
      dg(4,:,:)            = 1d0
      dge                  = 0d0     
      fint                 = 0d0
      f_ext                = 0d0
      a                    = 0d0
      a_hat                = 0d0
      lambda               = 0d0
      res                  = 0d0
      res_hat              = 0d0

      ! Initiate times
      stiffness_time       = 0d0
      assem_stiffness_time = 0d0
      assem_force_time     = 0d0
      spaputval_time       = 0d0
      solveq_time          = 0d0
      force_time           = 0d0
      computeplot_time     = 0d0
      mat_time             = 0d0
      plotvtk_time         = 0d0      

      return
   end subroutine init_NewtonEqIter

   subroutine update_coord(mesh)
      implicit none
      type(mesh_system), intent(inout)          :: mesh

      ! Update coordinates
      call updcoord(mesh%newcoord,mesh%coord,a)
      call coordxtr(mesh%newex,mesh%newey,mesh%newcoord,mesh%enod)

      return
   end subroutine update_coord

   subroutine loadloop(mesh, lssys, IMC_eps_star, pq, i_IMC, omp_run, input_location, newton_loop_conv)
      implicit none

      ! --- Intent inout ---
      type(plot_system), intent(inout)          :: pq
      type(mesh_system), intent(inout)          :: mesh

      ! Intent in      
      type(ls_system), intent(inout)            :: lssys
      integer, intent(in)                       :: i_IMC
      real(dp), intent(in)                      :: IMC_eps_star
      logical, intent(in)                       :: omp_run
      character(len=:), allocatable, intent(in) :: input_location
      logical, intent(inout)                    :: newton_loop_conv

      ! Subroutine variables      
      real(dp)                                  :: beta, t11, t22, dg_ifstar(9)
      integer                                   :: iload, ls=1

      if (i_IMC.le.1) then

         if (i_IMC.eq.1) then
            do iload = 1,ls
               write(*,*)
               write(*,*)'-------------------'
               write(*,"(A,I3,A,I3)") 'loadstep ', iload, ' of ', ls

               ! Newton iteration
               beta = real(iload)/ls
               write(*,'(A,F10.3)') 'beta: ', beta
               call get_dg_ifstar(dg_ifstar, IMC_eps_star, beta)
               call NewtonEqIter4(mesh, lssys, omp_run, pq, dg_ifstar, newton_loop_conv)

               if (newton_loop_conv.eqv..false.) then
                  return
               else
                  ! Accept new state
                  call J2iso_Cu_accept(dg)
                  call J2iso_Sn_accept(dg)
               endif

            enddo  
         endif
      else         
         beta = 1d0
         call get_dg_ifstar(dg_ifstar, IMC_eps_star, beta)
         call NewtonEqIter4(mesh, lssys, omp_run, pq, dg_ifstar, newton_loop_conv)

         if (newton_loop_conv.eqv..false.) then
            return
         else
            ! Accept new state
            call J2iso_Cu_accept(dg)
            call J2iso_Sn_accept(dg)
         endif
      endif

      ! -- Iteration loop finished --
      write(*,*)'End iteration'
      write(*,*)'-------------------'             

      ! Calculate cauchy stresses and save in plot_quantities
      call pushforward(pq%cau,es,dg,'j')

      ! --- Load loop finished ---    
      call update_coord(mesh)

      ! --- Compute plot quantities ---
      call compute_plot_quantities(pq, mesh, lssys)

      ! Send to stress state to matlab
      call clock_time(t11, omp_run)
      call write_stress_iter_to_matlab(pq, mesh, lssys, input_location, i_IMC)
      call clock_time(t22, omp_run)
      mat_time = mat_time + (t22-t11)
      
      return 
   end subroutine loadloop



   subroutine NewtonEqIter4(mesh, lssys, omp_run, pq, dg_ifstar, newton_loop_conv)
      ! Equilibrium iteration
      implicit none
   
      ! Intent inout
      type(plot_system), intent(inout) :: pq
      logical, intent(inout)           :: newton_loop_conv

      ! --- Intent in ---
      type(mesh_system), intent(in)    :: mesh
      type(ls_system), intent(in)      :: lssys      
      logical, intent(in)              :: omp_run
      real(dp), intent(in)             :: dg_ifstar(9)
      
       
      ! --- Subroutine variables ---
      real(dp)                         :: residual, residual_prev, ke(8,8), fe(8)
      integer                          :: ie, i, iter, igp, j, gg, grain_idx(lssys%ngrains), material
      real(dp)                         :: t11, t22
      ! integer                          :: inside_grain(1), closest_grain(1), other_grains(lssys%ngrains-1)
      ! integer                          :: inside_material, closest_material
      ! real(dp)                         :: spiola1(4), spiola2(4), esloc(4), phi1, phi2, D1(3,3), D2(3,3), Dloc(3,3)
      real(dp)                         :: dg_ifstar_old(9)
   

      ! Set newton_loop_conv and old quantities
      newton_loop_conv = .true.
      edold            = ed
      dgold            = dg
      dgeold           = dge
      aold             = a

      a_hatold  = a_hat
      lambdaold = lambda
      

      ! Define dg_ifstar for old imc
      dg_ifstar_old    = 0d0
      dg_ifstar_old(1) = 1d0
      dg_ifstar_old(5) = 1d0
      dg_ifstar_old(9) = 1d0


      ! Define grain idx
      do i=1,lssys%ngrains
         grain_idx(i) = i
      enddo

      ! Residual
      res_hat = 0d0
      ! res = 0d0
   
      ! Define residual (to get into iteration loop)
      residual      = 1d0
      residual_prev = 1d20
   
      write(*,*)'-------------------'
      write(*,*)'Begin iteration'
      iter = 1
      do while (residual.gt.1d-10)
         global_check = 0d0
   
         ! -- Global stiffness --
         if (iter.gt.1) then
            K%a = 0d0

            ! Time stiffness
            call clock_time(t11, omp_run)
            
            ! --- Compute material tangents  ---
            ! For each grain, compute material tangent at all gauss points in the mesh
            ! (Assume each grains can have different material properties)
            Dgp_g  = 0d0            
            !$omp parallel private(gg, ie, igp, material) shared(Dgp_g, dg, dge, dg_ifstar)
            !$omp do schedule(guided)
            do gg=1,lssys%ngrains

               ! Material of grain (Cu, IMC or Sn)
               material = lssys%material(gg)

               if (material.eq.1) then
                  ! -- Cu --
                  do ie=1,mesh%nelm
                     call dJ2iso_Cu('tl', Dgp_g(:,:,:,ie,gg), dg(:,:,ie), ie)
                  enddo
               elseif (material.eq.2) then
                  ! -- IMC --
                  do ie=1,mesh%nelm
                     do igp=1,mesh%nrgp
                        if (lssys%old_IMC_gp_flag(igp,ie)) then
                           call dneohooke('tl', Dgp_g(:,:,igp,ie,gg), dge(:,igp,ie), dg_ifstar_old)
                        else
                           call dneohooke('tl', Dgp_g(:,:,igp,ie,gg), dge(:,igp,ie), dg_ifstar)
                        endif
                     enddo
                  enddo 
               elseif (material.eq.3) then
                  ! -- Sn --                  
                  do ie=1,mesh%nelm
                     call dJ2iso_Sn('tl', Dgp_g(:,:,:,ie,gg), dg(:,:,ie), ie)
                  enddo                  
               endif
            enddo
            !$omp end do    
            !$omp end parallel                               

            ! Interpolate material tangents between the different grains
            Dgp = 0d0     
            do ie=1,mesh%nelm
               do igp=1,mesh%nrgp
                  do gg=1,lssys%ngrains
                     Dgp(:,:,igp,ie) = Dgp(:,:,igp,ie) + lssys%hphi_gp(igp,ie,gg)*Dgp_g(:,:,igp,ie,gg)
                  enddo
                  pq%E_save_gp(igp,ie) = norm2(Dgp(:,:,igp,ie))
               enddo
            enddo
            
        
            ! --- Assemble all stiffness ---

            ! Stiffness time
            call clock_time(t22, omp_run)
            stiffness_time = stiffness_time + (t22-t11)
   
            ! Global check
            ! call global_gpel_check2(global_check, lssys, mesh%nelm)
   
            ! Assemble stiffness              
            call clock_time(t11, omp_run)
            do i=1,mesh%nset
               !$omp parallel do private(ie,ke,j) shared(mesh,ep,Dgp,ed,es,K,i) schedule(guided)
               do j=1,mesh%counter_set(i)
                  ie = mesh%set_groups(j,i)
                  call c2dtl4_e(ke,mesh%coord(:,mesh%enod(:,ie)),ep,Dgp(:,:,:,ie),ed(:,ie),es([1,2,4],:,ie))
                  call assem(K,ke,mesh%enod(:,ie),mesh%dofnod)
               enddo
               !$omp end parallel do
            enddo
            call clock_time(t22, omp_run)
            assem_stiffness_time = assem_stiffness_time + (t22-t11)
   
            ! Expand K, put new K%a in K_hat, T_map is already in K_hat
            ! call clock_time(t11, omp_run)        
            call spaputval(K_hat, Kcell, K%a)
            ! call clock_time(t22, omp_run)          
            ! spaputval_time = spaputval_time + (t22-t11)
   
            ! Solve K_hat*a_hat=res_hat
            ! call clock_time(t11, omp_run)
            call solveq(K_hat,res_hat,mesh%bcnod,mesh%bcval,mesh%dofnod)
            ! call solveq(K,res,mesh%bcnod,mesh%bcval,mesh%dofnod)
            call clock_time(t22, omp_run)
            solveq_time = solveq_time + (t22-t11)
         else
            ! res = 0d0
            res_hat = 0d0

         endif
   
         ! Update a_hat and extract a and lambda
         ! a = a + res
         a_hat  = a_hat + res_hat
         a      = a_hat(1:mesh%ndof)
         lambda = a_hat(mesh%ndof+1:ndof_hat)
   
         ! Extract element displacements
         call extract(ed,a,mesh%enod,mesh%dofnod)
   
         ! Extract deformation gradient for all elms
         !$omp parallel do private(ie) shared(dg)
         do ie=1,mesh%nelm
            call c2dtl4_d(dg(:,:,ie), mesh%coord(:,mesh%enod(:,ie)), ed(:,ie))
         enddo
         !$omp end parallel do

         ! Compute dge
         do ie=1,mesh%nelm
            do igp=1,mesh%nrgp
               if (lssys%old_IMC_gp_flag(igp,ie)) then                  
                  call elastd(dge(:,igp,ie), dg(:,igp,ie), dg_ifstar_old)
               else
                  call elastd(dge(:,igp,ie), dg(:,igp,ie), dg_ifstar)
               endif
            enddo
         enddo
         
         ! -- Internal forces --
         fint         = 0d0
         global_check = 0d0
   
         ! Time force
         call clock_time(t11, omp_run)

         ! --- Compute stresses ---
         ! For each grain, compute stress tensor at all gauss points in the mesh
         ! (Assume each grains can have different material properties)
         es_g  = 0d0         
         !$omp parallel private(gg, ie, igp, material) shared(es_g, dg, dge, dg_ifstar)
         !$omp do schedule(guided)
         do gg=1,lssys%ngrains

            ! Material of grain gg (Cu, IMC, Sn)
            material = lssys%material(gg)

            if (material.eq.1) then
               ! -- Cu --               
               do ie=1,mesh%nelm
                  call J2iso_Cu('2ndPiola', es_g(:,:,ie,gg), dg(:,:,ie), ie)
               enddo  
            elseif (material.eq.2) then
               ! -- IMC --
               do ie=1,mesh%nelm       
                  do igp=1,mesh%nrgp  
                     if (lssys%old_IMC_gp_flag(igp,ie)) then
                        call neohooke('2ndPiola',es_g(:,igp,ie,gg),dge(:,igp,ie),dg_ifstar_old)
                     else
                        call neohooke('2ndPiola',es_g(:,igp,ie,gg),dge(:,igp,ie),dg_ifstar)
                     endif
                  enddo
               enddo    
            elseif (material.eq.3) then
               ! -- Sn --               
               do ie=1,mesh%nelm
                  call J2iso_Sn('2ndPiola', es_g(:,:,ie,gg), dg(:,:,ie), ie)
               enddo               
            endif
         enddo
         !$omp end do   
         !$omp end parallel
              

         ! Interpolate stresses between the different grains         
         es = 0d0
         do ie=1,mesh%nelm
            do igp=1,mesh%nrgp
               do gg=1,lssys%ngrains
                  es(:,igp,ie) = es(:,igp,ie) + lssys%hphi_gp(igp,ie,gg)*es_g(:,igp,ie,gg)
               enddo               
            enddo
         enddo               
                            
         ! Compute force time
         call clock_time(t22, omp_run)
         force_time = force_time + (t22-t11)
   
         ! Invoke global check assembled in el loops
         ! call global_gpel_check2(global_check, lssys, mesh%nelm)
   
         ! Assemble global internal forces
         call clock_time(t11, omp_run)
         do i=1,mesh%nset
            !$omp parallel do private(ie,fe,j) shared(mesh,ep,ed,es,fint,i) schedule(guided)
            do j=1,mesh%counter_set(i)
               ie = mesh%set_groups(j,i)
               call c2dtl4_f(fe,mesh%coord(:,mesh%enod(:,ie)),ep,ed(:,ie),es([1,2,4],:,ie))
               call insert(fint,fe,mesh%enod(:,ie),mesh%dofnod)
            enddo
            !$omp end parallel do
         enddo
         call clock_time(t22, omp_run)
         assem_force_time = assem_force_time + (t22-t11)
         
         ! Residual      
         ! res             = -fint
         ! res(mesh%bcdof) = 0d0
   
         ! Define res_hat
         res_hat                       = 0d0
         res_hat(1:mesh%ndof)          = -(matmul(T_map_sparse_transpose, lambda) + fint)
         res_hat(mesh%ndof+1:ndof_hat) = -(matmul(T_map_sparse,a) -  g)
         res_hat(mesh%bcdof)           = 0d0
   
         ! Norm of residual
         residual=dot_product(res_hat,res_hat)
         ! residual=dot_product(res,res)
         write(*,'(A20, E11.4)') 'Residual          : ', residual
   
         ! Check if residual is increasing
         if (residual.gt.residual_prev) then
            print *,'Residual increasing.'            
            ed     = edold
            dg     = dgold
            dge    = dgeold
            a      = aold

            ! call exit
            a_hat  = a_hatold
            lambda = lambdaold

            newton_loop_conv = .false.
            return
         endif

         ! Set previous residual
         residual_prev = residual
   
         ! Update iteration counter
         iter = iter + 1      
   
      enddo
   
      return
     end subroutine NewtonEqIter4
   
  
  
     subroutine global_gpel_check2(global_check, lssys, nelm)
      implicit none
   
      ! Intent in
      type(ls_system), intent(in) :: lssys
      real(dp), intent(in)        :: global_check(:,:)
      
      ! Subroutine variables
      integer                     :: ie, igp, nelm
      logical                     :: check_passed=.true.
   
  
      do ie=1,nelm
         do igp=1,4
            if (global_check(igp,ie).eq.1d0) then
               !print *, 'ie: ', ie
               ! print *, 'global check accepted: ', global_check(:,ie)
            else
               check_passed= .false.
               print *, 'global error at: ', igp,ie
               print *, 'global_check(:,ie): ', global_check(:,ie)
               call exit(0)
            end if
         enddo
      enddo

      ! if (check_passed) then
      !    print *, 'Global check passed'
      ! endif
   
      return
     end subroutine global_gpel_check2


  subroutine printTimes(exec_time, init_time, diffp1_generate_mesh_time, diffp1_read_mesh_time, diffp1_solve_time, input_location)
   ! Routine for printing times in Newton loop
   implicit none

   ! Intent in
   real(dp), intent(in) :: exec_time, init_time, diffp1_generate_mesh_time, diffp1_read_mesh_time, diffp1_solve_time
   character(len=:), allocatable, intent(in) :: input_location
   
   write(*,*)
   write(*,'(A79)') '--------------------------- Program execution times ---------------------------'
   print *, 'Execution time                      = ', exec_time, 'seconds'
   print *, 'Initiation time                     = ', init_time, 'seconds'
   print *, 'Initiation time (%)                 = ', 100*init_time/exec_time
   print *, 'diffusion generate mesh time        = ', diffp1_generate_mesh_time, 'seconds'
   print *, 'diffusion generate mesh time (%)    = ', 100*diffp1_generate_mesh_time/exec_time
   print *, 'diffusion read mesh time            = ', diffp1_read_mesh_time, 'seconds'
   print *, 'diffusion read mesh time (%)        = ', 100*diffp1_read_mesh_time/exec_time
   print *, 'diffusion solve time                = ', diffp1_solve_time, 'seconds'
   print *, 'diffusion solve time (%)            = ', 100*diffp1_solve_time/exec_time
   print *, 'Stiffness time                      = ', stiffness_time, 'seconds'
   print *, 'Stiffness time (%)                  = ', 100*stiffness_time/exec_time
   print *, 'Assem stiffness time                = ', assem_stiffness_time, 'seconds'
   print *, 'Assem stiffness time (%)            = ', 100*assem_stiffness_time/exec_time
   print *, 'Spaputval time                      = ', spaputval_time, 'seconds'
   print *, 'Spaputval time (%)                  = ', 100*spaputval_time/exec_time
   print *, 'Solveq time                         = ', solveq_time, 'seconds'
   print *, 'Solveq time (%)                     = ', 100*solveq_time/exec_time
   print *, 'Force time                          = ', force_time, 'seconds'
   print *, 'Force time (%)                      = ', 100*force_time/exec_time
   print *, 'Assem force time                    = ', assem_force_time, 'seconds'
   print *, 'Assem force time (%)                = ', 100*assem_force_time/exec_time
   print *, 'Send to matlab time                 = ', mat_time, 'seconds'
   print *, 'Send to matlab time (%)             = ', 100*mat_time/exec_time
   print *, 'Send to vtk time                    = ', plotvtk_time, 'seconds'
   print *, 'Send to vtk time (%)                = ', 100*plotvtk_time/exec_time
   
   !print *, 'Press Enter to continue.'; read(stdin,*)


   ! Save times
   call write_times_to_matlab(exec_time, init_time, stiffness_time, assem_stiffness_time, spaputval_time, solveq_time, force_time,&
   assem_force_time, mat_time, plotvtk_time, input_location)

   return
  end subroutine printTimes

end module rhapsloop

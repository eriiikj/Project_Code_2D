module ls_vp_routines
    ! Last modified 
    ! E. Jacobsson 2023-09-13
      
    ! mf_datatypes
    use mf_datatypes

    ! Whiskerlib
    use ls_types
    use diffusion_types
    use matrix_util

    implicit none

     ! This is an attempt to predefine some vectors for the 4-node element, to obtain speed
    double precision, parameter  :: G1=0.577350269189626D0
    double precision, parameter  :: xsi4(4)=(/-1D0,  1D0, -1D0, 1D0/)*G1
    double precision, parameter  :: eta4(4)=(/-1D0, -1D0,  1D0, 1D0/)*G1
    private xsi4, eta4, G1

    ! 4 Gauss points and 4 nodes
    real(dp), parameter :: NR4_col11(4) = (1D0-XSI4)*(1D0-ETA4)/4D0
    real(dp), parameter :: NR4_col22(4) = (1D0+XSI4)*(1D0-ETA4)/4D0
    real(dp), parameter :: NR4_col33(4) = (1D0+XSI4)*(1D0+ETA4)/4D0
    real(dp), parameter :: NR4_col44(4) = (1D0-XSI4)*(1D0+ETA4)/4D0
    real(dp), parameter :: NR4(4,4)    = reshape([NR4_col11, NR4_col22, NR4_col33, NR4_col44], [4,4])
    private NR4


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

subroutine lvlset2D4_global_vp(vp_e, coord, ed_e_gr, jgp_e_gr, agp_e_gr, bcvalgp_e_gr, Igrain_gp, lssys, diffsys, ie)
    ! --- ---
    implicit none

    ! Intent out
    real(dp), intent(inout)     :: vp_e(:,:)

    ! Intent in
    type(ls_system), intent(in) :: lssys
    type(diffusion_system), intent(in) :: diffsys
    real(dp), intent(in)        :: coord(:,:), ed_e_gr(:,:), jgp_e_gr(:,:), agp_e_gr(:,:), bcvalgp_e_gr(:,:)
    integer, intent(in)         :: Igrain_gp(:), ie

    ! Subroutine variables    
    real(dp)                    :: JT(8,2), JTinv(2,2), DNX(2,4), B(2,4), nj(2), normn, jdiff, cdiff, vel_cont, f_j, ni(2)
    real(dp)                    :: DetJ, inside_material, other_material, Vmi, Vmj
    integer                     :: GP_NR, indx(2), inside_grain, grain_idx(lssys%ngrains), k, other_grain
    integer, allocatable        :: other_grains(:)
    integer, parameter          :: NGP=4
    real(dp)                    :: xi, xj, ci, cj
    

    ! Define grain idx
    do k=1,lssys%ngrains
        grain_idx(k) = k
    enddo


    JT=MATMUL(DNR4,transpose(coord))
    vp_e = 0d0
    do GP_NR=1,NGP
        indx = (/ 2*gp_nr-1, 2*gp_nr /)
        detJ = det2(JT(indx,:))
        call inv2(JTinv,JT(indx,:))
        dNX  = matmul(JTinv,DNR4(indx,:))
        B    = dNX

        ! Inside grain at current gp
        inside_grain = Igrain_gp(GP_NR)

        ! Material of inside grain
        inside_material = lssys%material(inside_grain)

        ! Inside molar volume
        Vmi = diffsys%molar_volumes(inside_material)

        ! Inside molar fraction
        xi =  (bcvalgp_e_gr(GP_NR,inside_grain) - diffsys%thermo_parameterB(inside_material)) &
        /diffsys%thermo_parameterA(inside_material) + diffsys%thermo_parameterxb(inside_material)

        ! Inside concentration
        ci = xi/diffsys%molar_volumes(inside_material)

        ! Normal direction at gp for grain
        ni    = matmul(B,ed_e_gr(:,inside_grain))
        normn = norm2(ni)
        ni    = ni/normn

        ! Grains to sum over at current gp
        other_grains = pack(grain_idx, (grain_idx-inside_grain).ne.0)


        do k = 1,size(other_grains)

            ! Grain number
            other_grain = other_grains(k)
            
            ! Normal direction at gp for other grain
            nj    = matmul(B,ed_e_gr(:,other_grain))
            normn = norm2(nj)
            
            if (normn.gt.1d-13) then
                nj = nj/normn
            else
                print *, 'Level set normal is zero'
            endif

            ! Material other grain            
            other_material  = lssys%material(other_grain)   
            
            ! Other molar volume
            Vmj = diffsys%molar_volumes(other_material)

            if (inside_material.ne.other_material) then        
                
                ! Other molar fraction
                xj =  (bcvalgp_e_gr(GP_NR,other_grain) - diffsys%thermo_parameterB(other_material)) &
                /diffsys%thermo_parameterA(other_material) + diffsys%thermo_parameterxb(other_material)

                ! Other concentration
                cj = xj/diffsys%molar_volumes(other_material)

                ! Jdiff
                jdiff = (jgp_e_gr(GP_NR,other_grain) + jgp_e_gr(GP_NR,inside_grain))

                ! Cdiff
                ! cdiff = (bcvalgp_e_gr(GP_NR,other_grain)-bcvalgp_e_gr(GP_NR,inside_grain))
                cdiff = cj - ci
                
                ! Velocity contribution ij
                vel_cont = (1/cdiff)*jdiff
                
                
    
                if (abs(bcvalgp_e_gr(GP_NR,other_grain)-bcvalgp_e_gr(GP_NR,inside_grain)).lt.1d-12) then
                    ! ! print*, 'Dividing by zero in vp construction routine'
                    ! jdiff = (jgp_e_gr(GP_NR,other_grain) + jgp_e_gr(GP_NR,inside_grain))/10000000d0
                    ! ! print *, 'as'
                    ! ! call exit(0)
                    ! if (jdiff.ne.0d0 .and. abs(f_j).gt.0.5d0) then
                    !     print *, 'as'
                    ! endif
                    print *, 'jdiff inf'
                endif
            else
                vel_cont = 0d0                
            endif

            ! s = sign(1d-1,bcvalgp_e_gr(GP_NR,other_grain)-bcvalgp_e_gr(GP_NR,inside_grain))
            ! jdiff = s*(jgp_e_gr(GP_NR,inside_grain) + jgp_e_gr(GP_NR,other_grain))


            ! if (inside_grain.eq.2) then
            !     print *, 'a'
            ! endif

            ! if (abs((coord(1,1)-0.000475d0)).lt.1d-5 .and. abs((coord(2,1)-0.0005205d0)).lt.1d-5) then
            !     print *, 'as'
            ! endif

            
            ! Decaying function       
            f_j = exp(-lssys%alpha*abs(agp_e_gr(GP_NR,other_grain)))

            ! Define vector vpvec  
            vp_e(:,GP_NR) = vp_e(:,GP_NR) + f_j*vel_cont*nj

            if (isnan(vel_cont)) then
                print *, 'vel_cont is NAN'
            endif

        enddo


    enddo

    ! Deallocate
    deallocate(other_grains)


    return
end subroutine lvlset2D4_global_vp


subroutine lvlset2D4_global_vp2(vp_e, coord, ed_e_gr, jgp_e_gr, agp_e_gr, Igrain_gp, lssys, diffsys, ie)
    ! --- ---
    implicit none

    ! Intent out
    real(dp), intent(inout)     :: vp_e(:,:)

    ! Intent in
    type(ls_system), intent(in) :: lssys
    type(diffusion_system), intent(in) :: diffsys
    real(dp), intent(in)        :: coord(:,:), ed_e_gr(:,:), jgp_e_gr(:,:), agp_e_gr(:,:)
    integer, intent(in)         :: Igrain_gp(:), ie

    ! Subroutine variables    
    real(dp)                    :: JT(8,2), JTinv(2,2), DNX(2,4), B(2,4), nj(2), normn, jdiff, cdiff, vel_cont, f_j, ni(2)
    real(dp)                    :: DetJ, inside_material, other_material, Vmi, Vmj
    integer                     :: GP_NR, indx(2), inside_grain, grain_idx(lssys%ngrains), k, other_grain, Ograin(1)
    integer, allocatable        :: other_grains(:)
    integer, parameter          :: NGP=4
    real(dp)                    :: xi, xj, ci, cj
    

    ! Define grain idx
    do k=1,lssys%ngrains
        grain_idx(k) = k
    enddo


    JT=MATMUL(DNR4,transpose(coord))
    vp_e = 0d0
    do GP_NR=1,NGP
        indx = (/ 2*gp_nr-1, 2*gp_nr /)
        detJ = det2(JT(indx,:))
        call inv2(JTinv,JT(indx,:))
        dNX  = matmul(JTinv,DNR4(indx,:))
        B    = dNX

        ! Inside grain at current gp
        inside_grain = Igrain_gp(GP_NR)

        ! Material of inside grain
        inside_material = lssys%material(inside_grain)

        ! Normal direction at gp for inside grain
        ni    = matmul(B,ed_e_gr(:,inside_grain))
        normn = norm2(ni)
        ni    = ni/normn

        ! Grains to sum over at current gp
        other_grains = pack(grain_idx, (grain_idx-inside_grain).ne.0)

        ! Closest other grain
        Ograin = minloc(agp_e_gr(GP_NR,other_grains))
        other_grain = other_grains(Ograin(1))

        ! do k = 1,size(other_grains)

        !     ! Grain number
        !     other_grain = other_grains(k)

            ! Material other grain            
            other_material  = lssys%material(other_grain) 

            if (inside_material.ne.other_material) then
            
                ! Normal direction at gp for other grain
                nj    = matmul(B,ed_e_gr(:,other_grain))
                normn = norm2(nj)
                
                if (normn.gt.1d-13) then
                    nj = nj/normn
                else
                    print *, 'Level set normal is zero'
                endif
                
                ! Equilibrium concentrations for inside grain and other grain for the inside/other phase interface
                xi = diffsys%eq_x(inside_material, other_material)
                ci = xi/diffsys%molar_volumes(inside_material)
                xj = diffsys%eq_x(other_material, inside_material)
                cj = xj/diffsys%molar_volumes(other_material)

                ! Jdiff in Gauss point
                jdiff = (jgp_e_gr(GP_NR,other_grain) + jgp_e_gr(GP_NR,inside_grain))

                ! Cdiff in Gauss point            
                cdiff = cj - ci
                if (cdiff.eq.0d0) then
                    print *, 'cdiff is 0'
                    print *, 'coord: '          , coord
                    print *, 'inside_grain: '   , inside_grain
                    print *, 'inside material: ', inside_material
                    print *, 'other_grain: '    , other_grain
                    print *, 'other material: ' , other_material
                    print *, 'cdiff: '          , cdiff
                    print *, 'jdiff: '          , jdiff
                    call exit(0)
                endif
                
                ! Velocity contribution ij
                vel_cont = (1/cdiff)*jdiff

                if (isnan(vel_cont)) then
                    print *, 'vel_cont is NAN'                
                    print *, 'coord: '          , coord
                    print *, 'inside_grain: '   , inside_grain
                    print *, 'inside material: ', inside_material
                    print *, 'other_grain: '    , other_grain
                    print *, 'other material: ' , other_material
                    print *, 'cdiff: '          , cdiff
                    print *, 'jdiff: '          , jdiff
                    call exit(0)
                endif

                ! if (abs((coord(1,1)-0.000524966d0)).lt.1d-6 .and. abs((coord(2,1)-0.000624925d0)).lt.1d-6) then
                !     print *, 'GP: ', GP_NR
                !     print *, 'Xcoord: ', coord(1,:)                    
                !     print *, 'Ycoord: ', coord(2,:)                    
                !     print *, 'other grain: ', other_grain
                !     print *, 'cdiff: '          , cdiff
                !     print *, 'jdiff: '          , jdiff
                !     print *, 'vel_cont: ', vel_cont
                !     print *, 'Nj: ', nj
                ! endif

                ! Decaying function       
                f_j = exp(-lssys%alpha*abs(agp_e_gr(GP_NR,other_grain)))

                ! Define vector vpvec  
                vp_e(:,GP_NR) = vp_e(:,GP_NR) + f_j*vel_cont*nj
            endif
    
        ! enddo

    enddo

    ! Deallocate
    deallocate(other_grains)

    return
end subroutine lvlset2D4_global_vp2




subroutine lvlset2D4vp_exp_brown(vpvec,coord,ed,P0, model_width, cu_height)
    ! Routine for computing vp for a 4-node element using experimental data related to radius of circle
    implicit none
    
    real(dp), intent(inout)        :: vpvec(:,:)
    
    double precision, intent(in)   :: coord(:,:), ed(:), model_width, cu_height
    real(dp), intent(in)           :: P0(2)
    
    real(dp)                       :: normal(2), gp_x, gp_y, threshold, r, vp_length, p
    
    double precision               :: JT(8,2), JTinv(2,2), DNX(2,4), B(2,4), N(4), x0, y0
    double precision               :: DetJ
    integer                        :: GP_NR, indx(2)
    integer, parameter             :: NGP=4
    
    ! Center of circle
    x0 = P0(1)
    y0 = P0(2)
    
    JT = MATMUL(DNR4,transpose(coord))
    do GP_NR=1,NGP
    
        indx = (/ 2*gp_nr-1, 2*gp_nr /)
        detJ = det2(JT(indx,:))
        call inv2(JTinv,JT(indx,:))
        dNX = matmul(JTinv,DNR4(indx,:))
        B = dNX
    
        ! Normal direction in gp (gradient of a_ls)
        normal = matmul(B,ed)
        normal = normal/norm2(normal)
    
        ! N
        N = NR4(GP_NR,:)
    
        ! Coordinates of gauss point (in global domain)
        gp_x = dot_product(N,coord(1,:))
        gp_y = dot_product(N,coord(2,:))
    
        ! Radius of circle at given gauss point
        r = sqrt((gp_x-x0)**2 + (gp_y-y0)**2)
    
        ! Lenght of vp vector (from experiments)
        p         = sqrt(2d0*model_width*0.014d-3/pi)
        vp_length = 0.3d0*p**(10d0/3d0)*r**(-7d0/3d0)
    
        ! Exclude high velocity terms
        threshold = 10d-6
        if (vp_length.gt.threshold) then
        vp_length = threshold
        endif

        if (gp_y.lt.cu_height) then
            vp_length = 0d0
        endif



    
        ! Define vector vp
        vpvec(:,gp_nr) = vp_length*normal
    
    enddo
    
    return
    end subroutine lvlset2D4vp_exp_brown




end module ls_vp_routines
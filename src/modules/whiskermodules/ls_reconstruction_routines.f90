module ls_reconstruction_routines
    ! Last modified 
    ! E. Jacobsson 2023-09-06
      
    ! mf_datatypes
    use mf_datatypes

    ! Somelib
    use mesh_module

    ! whiskerlib
    use ls_types

    implicit none

contains

subroutine gb_reconstruction(lssys,mesh)
    ! --- Routine for fixing grain boundaries s.t they are compatible (check elements containing two ls isocontours) ---
    implicit none

    ! Intent inout
    type(ls_system), intent(inout) :: lssys

    ! Intent in
    type(mesh_system), intent(in)  :: mesh

    ! Subroutine variables
    integer, allocatable           :: gb_elms(:)
    integer                        :: grain_idx(lssys%ngrains), line_elm_idx(lssys%nseg_alloc)
    integer                        :: ngb_elms, gb_ie, ls_at_gb(2), ga, gb, g_rowsa(2), g_rowsb(2), iseg_ga(1), iseg_gb(1)
    integer                        :: i, k, ierr
    real(dp)                       :: ex_e(4), ey_e(4), ed_e_ga(4), ed_e_gb(4), ed_e(4), line_exeg(2), line_eyeg(2)

    ! Determine gb elements
    ngb_elms = count(count(lssys%int_elms,2).eq.2)    
    allocate(gb_elms(ngb_elms),stat=ierr)
    gb_elms = pack(mesh%elmidx,count(lssys%int_elms,2).eq.2)  

    ! Initialize grain_idx and line_elm_idx arrays
    grain_idx    = [(i, i = 1, lssys%ngrains)]
    line_elm_idx = [(i, i = 1, mesh%nelm)]
      

    ! Loop through all gb elms and fix
    do k = 1, ngb_elms

        ! Grain boundary element
        gb_ie = gb_elms(k)
    
        ! Identify grains to be matched at grain boundary        
        ls_at_gb = pack(grain_idx, lssys%int_elms(gb_ie, :))
        ga = ls_at_gb(1)
        gb = ls_at_gb(2)
        g_rowsa = [(2 * (ga - 1) + 1), (2 * (ga - 1) + 2)]
        g_rowsb = [(2 * (gb - 1) + 1), (2 * (gb - 1) + 2)]
    
        ! Make interface reconstruction of grain boundary element
        ex_e    = mesh%ex(:,gb_ie)
        ey_e    = mesh%ey(:,gb_ie)
        ed_e_ga = lssys%ed(:, gb_ie, ga)
        ed_e_gb = lssys%ed(:, gb_ie, gb)
        ed_e    = 0.5d0*(ed_e_ga - ed_e_gb) ! New ed to be used to interpolate new line
        
        ! Call the get_single_line function to compute line_exeg and line_eyeg for ed_e
        call get_single_line(line_exeg,line_eyeg,ex_e,ey_e,ed_e)
        
        ! Update line_ex_gr and line_ey_gr arrays
        iseg_ga = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(ga), ga) - gb_ie .eq. 0)
        iseg_gb = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gb), gb) - gb_ie .eq. 0)
        lssys%line_ex(iseg_ga(1), g_rowsa(1):g_rowsa(2)) = line_exeg
        lssys%line_ey(iseg_ga(1), g_rowsa(1):g_rowsa(2)) = line_eyeg
        lssys%line_ex(iseg_gb(1), g_rowsb(1):g_rowsb(2)) = line_exeg
        lssys%line_ey(iseg_gb(1), g_rowsb(1):g_rowsb(2)) = line_eyeg        
    end do
      
    ! Deallocate
    deallocate(gb_elms)

    return
end subroutine gb_reconstruction


subroutine gb_reconstruction_spatial(lssys,mesh)
    ! --- Routine for fixing grain boundaries s.t they are compatible (check elements containing two ls isocontours) ---
    implicit none

    ! Intent inout
    type(ls_system), intent(inout) :: lssys

    ! Intent in
    type(mesh_system), intent(in)  :: mesh

    ! Subroutine variables
    integer, allocatable           :: gb_elms(:)
    integer                        :: grain_idx(lssys%ngrains), line_elm_idx(lssys%nseg_alloc)
    integer                        :: ngb_elms, gb_ie, ls_at_gb(2), ga, gb, g_rowsa(2), g_rowsb(2), iseg_ga(1), iseg_gb(1)
    integer                        :: i, k, ierr
    real(dp)                       :: ex_e(4), ey_e(4), ed_e_ga(4), ed_e_gb(4), ed_e(4), line_exeg(2), line_eyeg(2)

    ! Determine gb elements
    ngb_elms = count(count(lssys%int_elms,2).eq.2)    
    allocate(gb_elms(ngb_elms),stat=ierr)
    gb_elms = pack(mesh%elmidx,count(lssys%int_elms,2).eq.2)  

    ! Initialize grain_idx and line_elm_idx arrays
    grain_idx    = [(i, i = 1, lssys%ngrains)]
    line_elm_idx = [(i, i = 1, mesh%nelm)]

    ! Loop through all gb elms and fix
    do k = 1, ngb_elms

        ! Grain boundary element
        gb_ie = gb_elms(k)
    
        ! Identify grains to be matched at grain boundary
        ls_at_gb = pack(grain_idx, lssys%int_elms(gb_ie, :))
        ga = ls_at_gb(1)
        gb = ls_at_gb(2)
        g_rowsa = [(2 * (ga - 1) + 1), (2 * (ga - 1) + 2)]
        g_rowsb = [(2 * (gb - 1) + 1), (2 * (gb - 1) + 2)]
    
        ! Make interface reconstruction of grain boundary element
        ex_e    = mesh%newex(:,gb_ie)
        ey_e    = mesh%newey(:,gb_ie)
        ed_e_ga = lssys%ed(:, gb_ie, ga)
        ed_e_gb = lssys%ed(:, gb_ie, gb)
        ed_e    = 0.5d0*(ed_e_ga - ed_e_gb) ! New ed to be used to interpolate new line
        
        ! Call the get_single_line function to compute line_exeg and line_eyeg for ed_e
        call get_single_line(line_exeg,line_eyeg,ex_e,ey_e,ed_e)
        
        ! Update line_ex_gr and line_ey_gr arrays
        iseg_ga = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(ga), ga) - gb_ie .eq. 0)
        iseg_gb = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gb), gb) - gb_ie .eq. 0)
        lssys%line_ex(iseg_ga(1), g_rowsa(1):g_rowsa(2)) = line_exeg
        lssys%line_ey(iseg_ga(1), g_rowsa(1):g_rowsa(2)) = line_eyeg
        lssys%line_ex(iseg_gb(1), g_rowsb(1):g_rowsb(2)) = line_exeg
        lssys%line_ey(iseg_gb(1), g_rowsb(1):g_rowsb(2)) = line_eyeg        
    end do
      
    ! Deallocate
    deallocate(gb_elms)

    return
end subroutine gb_reconstruction_spatial

subroutine get_single_line(line_exeg,line_eyeg,ex_e,ey_e,ed_e)
    ! --- Routine for computing a single line of an element based on ed_e
    implicit none

    ! Intent out
    real(dp), intent(out) :: line_exeg(2),line_eyeg(2)

    ! Intent in
    real(dp), intent(in) :: ex_e(4), ey_e(4), ed_e(4)

    ! Subroutine variables
    integer                        :: counter, inod
    real(dp)                       :: isov=0d0, intersection_points(4,2)=0d0, edges_intersected(4) = 0d0
    real(dp)                       :: x1,y1,x2,y2,t1,t2,alpha,intersection_x,intersection_y


    ! Loop through all edges of the element and find intersection points
    counter = 0
    do inod=1,4
        x1 = ex_e(inod)
        y1 = ey_e(inod)
        x2 = ex_e(mod(inod, 4) + 1)
        y2 = ey_e(mod(inod, 4) + 1)

        ! Nodal values along element edge
        t1 = ed_e(inod)
        t2 = ed_e(mod(inod, 4) + 1)

        ! Check if the contour crosses the current edge
        if ((t1-isov) * (t2-isov) .lt. 0d0) then
            
            ! Interpolate the intersection point
            alpha = (isov - t1) / (t2 - t1)
            intersection_x = x1 + alpha * (x2 - x1)
            intersection_y = y1 + alpha * (y2 - y1)
            
            ! Store the segment
            counter = counter + 1
            intersection_points(inod,:) = [intersection_x, intersection_y]
            edges_intersected(counter)  = inod
        endif

    enddo

    ! Find nonzero intersection_points. OBS either 2 or 4 intersected edges
    if (counter.eq.2) then
            
        ! One line segment in element        
        line_exeg     = intersection_points(edges_intersected(1:counter),1)
        line_eyeg     = intersection_points(edges_intersected(1:counter),2)

    elseif (counter.eq.4) then
        print *, 'Error: A single level set function has two line segments in a grain boundary element'
    else
        print *, 'Error in get_single_line in gb_reconstruction'
    endif

    return
end subroutine get_single_line


subroutine tp_reconstruction(lssys,mesh,critical_length)
    ! --- Routine for fixing triple junctions s.t they are compatible (check elements containing three ls isocontours) ---
    implicit none

    ! Intent inout
    type(ls_system), intent(inout) :: lssys

    ! Intent in
    type(mesh_system), intent(in)  :: mesh
    real(dp)                       :: critical_length

    ! Subroutine variables    
    integer, allocatable   :: tp_elms(:), fp_elms(:), tp_elms_reordered(:), tmp(:), ls_at_tp(:)
    integer                :: grain_idx(lssys%ngrains), line_elm_idx(lssys%nseg_alloc), nls_at_tp
    integer                :: ntp_elms, tp_ie, ga, gb, gc, gd, g_rowsa(2), g_rowsb(2), g_rowsc(2), g_rowsd(2)
    integer                :: nfp_elms, rm_idx, junc
    integer                :: i, ierr
    real(dp)               :: ex_e(4), ey_e(4), points(6,2)
    real(dp)               :: points3(12,2), points4(4,2), points4b(4,2)
    logical                :: skip_tp

    ! Short 
    integer                :: nsega, nsegb, nsegc, nsegd, nunique
    integer, allocatable   :: lsa_conn_lines(:), lsb_conn_lines(:), lsc_conn_lines(:), lsd_conn_lines(:),unique_points_idx(:)
    integer, allocatable   :: gg_tplines(:), incompatible_edges(:), duplicate_points_idx(:), keep_lines(:)
    real(dp), allocatable  :: line_exa(:,:), line_eya(:,:), line_exb(:,:), line_eyb(:,:), line_exc(:,:), line_eyc(:,:)
    real(dp), allocatable  :: line_exd(:,:), line_eyd(:,:)
    logical, allocatable   :: lsa_log_conn(:,:), lsb_log_conn(:,:), lsc_log_conn(:,:), lsd_log_conn(:,:), pe(:,:)
    real(dp), dimension(2) :: lsa_connA,lsa_iconnA,lsa_connB,lsa_iconnB,lsa_connC,lsa_iconnC,lsa_connD,lsa_iconnD
    real(dp), dimension(2) :: lsb_connA,lsb_iconnA,lsb_connB,lsb_iconnB,lsb_connC,lsb_iconnC,lsb_connD,lsb_iconnD
    real(dp), dimension(2) :: lsc_connA,lsc_iconnA,lsc_connB,lsc_iconnB,lsc_connC,lsc_iconnC,lsc_connD,lsc_iconnD
    real(dp), dimension(2) :: lsd_connA,lsd_iconnA,lsd_connB,lsd_iconnB,lsd_connC,lsd_iconnC,lsd_connD,lsd_iconnD
    real(dp), dimension(2) :: A, B, C, D, E, F, PP, S, compatible_point, Aprev, Bprev, Cprev, M
    real(dp)               :: minx,maxx,miny,maxy
    logical                :: step_junc, is_inside, fp_elm
    integer                :: idx3(3) = [1,2,3], gi(1), gidx3(3), gg, g_rows(2), gg_rows(2), gother2(2), p1(1), jj, g, minidx(1)
    integer                :: idx4(4) = [1,2,3,4], nedges, npoints, incompatible_edge, idx6(6)=[1,2,3,4,5,6], steps, ival(1), ie
    integer                :: nkeep, lsa_rm_lines(2), lsb_rm_lines(2), lsc_rm_lines(2)
    logical                :: short_tpline
        


    ! Zero out
    lssys%rm_lines     = 0
    lssys%rm_lines_tmp = 0
    lssys%tplines      = .false.

    ! Determine tp elements
    ntp_elms = count(count(lssys%int_elms,2).ge.3)    
    allocate(tp_elms(ntp_elms),stat=ierr)
    tp_elms = pack(mesh%elmidx,count(lssys%int_elms,2).ge.3)

    ! Check if any elements contain four ls and reorder tpelms s.t. these fp_elms are placed last in tp_elms
    nfp_elms = count(count(lssys%int_elms,2).eq.4)  
    
    ! Reorder tpelms
    if (nfp_elms .gt. 0) then
        allocate(fp_elms(nfp_elms),stat=ierr)
        fp_elms = pack(mesh%elmidx,count(lssys%int_elms,2).eq.4)
        allocate(tp_elms_reordered(ntp_elms),stat=ierr)
        call setdiff(tmp,tp_elms,fp_elms)
        tp_elms_reordered(1:ntp_elms-nfp_elms)          = tmp
        tp_elms_reordered(ntp_elms-nfp_elms+1:ntp_elms) = fp_elms
        tp_elms                                         = tp_elms_reordered
        print *, 'OBS 4 ls active in at least one element!'
    end if

    ! Initialize grain_idx and line_elm_idx arrays
    grain_idx    = [(i, i = 1, lssys%ngrains)]
    line_elm_idx = [(i, i = 1, lssys%nseg_alloc)]

    ! Allocate rm_lines
    rm_idx = 0

    do junc=1,ntp_elms

        ! Triple junction element
        tp_ie = tp_elms(junc)

        ! Identify grains to be matched at triple junction
        nls_at_tp = count(lssys%int_elms(tp_ie, :))
        allocate(ls_at_tp(nls_at_tp),stat=ierr)  
        ls_at_tp = pack(grain_idx, lssys%int_elms(tp_ie, :))
        ga       = ls_at_tp(1)
        gb       = ls_at_tp(2)
        gc       = ls_at_tp(3)
        g_rowsa  = [(2 * (ga - 1) + 1), (2 * (ga - 1) + 2)]
        g_rowsb  = [(2 * (gb - 1) + 1), (2 * (gb - 1) + 2)]
        g_rowsc  = [(2 * (gc - 1) + 1), (2 * (gc - 1) + 2)]
        skip_tp  = .false.
        short_tpline=.false.

        ! Check if four point junction element
        fp_elm = count(lssys%int_elms(tp_ie,:)).eq.4

        ! Short form
        nsega     = lssys%line_seg(ga)
        nsegb     = lssys%line_seg(gb)
        nsegc     = lssys%line_seg(gc)
        allocate(line_exa(nsega,2),stat=ierr)
        allocate(line_eya(nsega,2),stat=ierr)
        allocate(line_exb(nsegb,2),stat=ierr)
        allocate(line_eyb(nsegb,2),stat=ierr)
        allocate(line_exc(nsegc,2),stat=ierr)
        allocate(line_eyc(nsegc,2),stat=ierr)
        line_exa  = lssys%line_ex(1:nsega,g_rowsa)
        line_eya  = lssys%line_ey(1:nsega,g_rowsa)
        line_exb  = lssys%line_ex(1:nsegb,g_rowsb)
        line_eyb  = lssys%line_ey(1:nsegb,g_rowsb)
        line_exc  = lssys%line_ex(1:nsegc,g_rowsc)
        line_eyc  = lssys%line_ey(1:nsegc,g_rowsc)


        ! Ls a
        call get_connecting_points_ie(mesh%ex, mesh%ey, tp_ie, line_elm_idx, line_exa, line_eya, nsega, lssys%line_elms(:,ga), &
        lsa_connA, lsa_iconnA, lsa_connB, lsa_iconnB, lsa_connC, lsa_iconnC, lsa_connD, lsa_iconnD, lsa_conn_lines, lsa_log_conn)

        ! Ls b
        call get_connecting_points_ie(mesh%ex, mesh%ey, tp_ie, line_elm_idx, line_exb, line_eyb, nsegb, lssys%line_elms(:,gb), &
        lsb_connA, lsb_iconnA, lsb_connB, lsb_iconnB, lsb_connC, lsb_iconnC, lsb_connD, lsb_iconnD, lsb_conn_lines, lsb_log_conn)

        ! Ls c
        call get_connecting_points_ie(mesh%ex, mesh%ey, tp_ie, line_elm_idx, line_exc, line_eyc, nsegc, lssys%line_elms(:,gc), &
        lsc_connA, lsc_iconnA, lsc_connB, lsc_iconnB, lsc_connC, lsc_iconnC, lsc_connD, lsc_iconnD, lsc_conn_lines, lsc_log_conn)

                
        ! Determine unique connecting points
        points(1,:) = lsa_connA
        points(2,:) = lsa_connB
        points(3,:) = lsb_connA
        points(4,:) = lsb_connB
        points(5,:) = lsc_connA
        points(6,:) = lsc_connB
        
        
        call get_unique_idx(unique_points_idx, points)
        nunique = size(unique_points_idx)
        
        
        ! Find unique points A, B, and C
        if (nunique.eq.3) then
            A = points(unique_points_idx(1),:)
            B = points(unique_points_idx(2),:)
            C = points(unique_points_idx(3),:)

            
            if (size(lsa_conn_lines).eq.1 .or. size(lsb_conn_lines).eq.1 .or. size(lsc_conn_lines).eq.1) then
                ! Tp at edge of domain. Dont enter stepping loop
                F         = (A + B + C)/3d0
                step_junc = .false.
            else
                ! Compute isogonic point
                call get_isogonic(F,step_junc,A,B,C)

                ! if ((step_junc.eqv..false.) .and. ((norm2(A-F).lt.critical_length) .or. (norm2(B-F).lt.critical_length) .or. &
                ! (norm2(C-F).lt.critical_length))) then
                !     print *, 'Too short line in tp elm'
                !     step_junc    = .true.
                !     short_tpline = .true.
                ! endif

            endif
            
        elseif(fp_elm) then
            ! 4 point junction - take special care

            ! Dont consider as a triple junction element. And don't step.
            skip_tp   = .true.
            step_junc = .false.

            ! Identify grains to be matched at 4p junction
            gd        = ls_at_tp(4)
            g_rowsd   = [2*(gd-1) + 1,2*(gd-1) + 2]
            nsegd     = lssys%line_seg(gd)
            allocate(line_exc(nsega,2),stat=ierr)
            allocate(line_eyc(nsega,2),stat=ierr)
            line_exd  = lssys%line_ex(1:nsegd,g_rowsd)
            line_eyd  = lssys%line_ey(1:nsegd,g_rowsd)

            ! Ls d
            call get_connecting_points_ie(mesh%ex, mesh%ey, tp_ie, line_elm_idx, line_exd, line_eyd, nsegd,lssys%line_elms(:,gd), &
            lsd_connA, lsd_iconnA, lsd_connB, lsd_iconnB, lsd_connC, lsd_iconnC, lsd_connD, lsd_iconnD, lsd_conn_lines,lsd_log_conn)

            ! Add lines in triple junction element to rm_lines
            rm_idx = rm_idx + 1
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(ga),ga) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,ga) = ival(1)
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gb),gb) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gb) = ival(1)
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gc),gc) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gc) = ival(1)
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gd),gd) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gd) = ival(1)

            ! Set isogonic point as mean of all input points
            F = (lsa_connA + lsa_connB + lsb_connA + lsb_connB + lsc_connA + lsc_connB + lsd_connA + lsd_connB)/8d0
            
            ! Add isogonic point to level set function a
            call add_isogonic(lssys%line_ex(:,g_rowsa(1):g_rowsa(2)),lssys%line_ey(:,g_rowsa(1):g_rowsa(2)),lssys%tplines(:,ga), &
            lsa_connA,lsa_connB,F,nsega,lssys%int_elms(:,ga),mesh)

            ! Add isogonic point to level set function b
            call add_isogonic(lssys%line_ex(:,g_rowsb(1):g_rowsb(2)),lssys%line_ey(:,g_rowsb(1):g_rowsb(2)),lssys%tplines(:,gb), &
            lsb_connA,lsb_connB,F,nsegb,lssys%int_elms(:,gb),mesh)

            ! Add isogonic point to level set function c
            call add_isogonic(lssys%line_ex(:,g_rowsc(1):g_rowsc(2)),lssys%line_ey(:,g_rowsc(1):g_rowsc(2)),lssys%tplines(:,gc), &
            lsc_connA,lsc_connB,F,nsegc,lssys%int_elms(:,gc),mesh)

            ! Add isogonic point to level set function d
            call add_isogonic(lssys%line_ex(:,g_rowsd(1):g_rowsd(2)),lssys%line_ey(:,g_rowsd(1):g_rowsd(2)),lssys%tplines(:,gd), &
            lsd_connA,lsd_connB,F,nsegd,lssys%int_elms(:,gd),mesh)

            ! Update line_seg counter
            lssys%line_seg(ga) = nsega + 2
            lssys%line_seg(gb) = nsegb + 2
            lssys%line_seg(gc) = nsegc + 2
            lssys%line_seg(gd) = nsegd + 2
        
            ! Update rm idx
            rm_idx = rm_idx + 1

        elseif (nunique.gt.3) then
            ! More than 3 connecting points

            if (size(lsa_conn_lines).eq.4 .or. size(lsb_conn_lines).eq.4 .or. size(lsc_conn_lines).eq.4) then                
                ! One ls has two line segments in the element. Don't consider as a triple junction.
                points3(1,:)  = lsa_connA
                points3(2,:)  = lsa_connB
                points3(3,:)  = lsa_connC
                points3(4,:)  = lsa_connD
                points3(5,:)  = lsb_connA
                points3(6,:)  = lsb_connB
                points3(7,:)  = lsb_connC
                points3(8,:)  = lsb_connD
                points3(9,:)  = lsc_connA
                points3(10,:) = lsc_connB
                points3(11,:) = lsc_connC
                points3(12,:) = lsc_connD

                ! Determine what level set (assume it can only be one) that has two active lines at the element
                gi = pack(idx3, [size(lsa_conn_lines).eq.4, size(lsb_conn_lines).eq.4, size(lsc_conn_lines).eq.4])
                
                ! Connecting points
                A = points3((gi(1)-1)*4+1,:)
                B = points3((gi(1)-1)*4+2,:)
                C = points3((gi(1)-1)*4+3,:)
                D = points3((gi(1)-1)*4+4,:)

                points4(1,:) = A
                points4(2,:) = B
                points4(3,:) = C
                points4(4,:) = D

                ! Determine gg (level set with two line segments in tp_ie)
                gidx3    = [ga, gb, gc]
                gg       = gidx3(gi(1))
                gg_rows  = [2*(gg-1) + 1,2*(gg-1) + 2]

                ! Find other 2 grains at tp_ie
                gother2  = pack(gidx3,.not.[size(lsa_conn_lines).eq.4, size(lsb_conn_lines).eq.4, size(lsc_conn_lines).eq.4])

                ! Find which line segments of gg is in tp_ie
                gg_tplines = pack(line_elm_idx(1:lssys%line_seg(gg)),lssys%line_elms(1:lssys%line_seg(gg),gg).eq.tp_ie)
                

                ! Change line segments in gother(1) and gg to connecting points
                do jj=1,2
                    g        = gother2(jj)
                    g_rows   = [2*(g-1) + 1,2*(g-1) + 2]; 
                    p1       = pack(line_elm_idx(1:lssys%line_seg(g)),lssys%line_elms(1:lssys%line_seg(g),g).eq.tp_ie) ! Size 1
                    E        = [lssys%line_ex(p1(1),g_rows(1)), lssys%line_ey(p1(1),g_rows(1))]
                    F        = [lssys%line_ex(p1(1),g_rows(2)), lssys%line_ey(p1(1),g_rows(2))]                                
                    points4b = points4

                    ! Determine what point E is closest to
                    PP = E
                    minidx = minloc([norm2(PP-points4b(1,:)), norm2(PP-points4b(2,:)), &
                                     norm2(PP-points4b(3,:)), norm2(PP-points4b(4,:))])                                    
                    S = points4b(minidx(1),:)
                    points4b(minidx,:) = 1d8

                    ! Change lines
                    lssys%line_ex(p1(1),g_rows(1))           = S(1)
                    lssys%line_ey(p1(1),g_rows(1))           = S(2)
                    lssys%line_ex(gg_tplines(jj),gg_rows(1)) = S(1)
                    lssys%line_ey(gg_tplines(jj),gg_rows(1)) = S(2)
            
                    ! Determine what point E is closest to
                    PP = F
                    minidx = minloc([norm2(PP-points4b(1,:)), norm2(PP-points4b(2,:)), &
                                     norm2(PP-points4b(3,:)), norm2(PP-points4b(4,:))])                                    
                    S = points4b(minidx(1),:)

                    ! Change lines
                    lssys%line_ex(p1(1),g_rows(2))           = S(1)
                    lssys%line_ey(p1(1),g_rows(2))           = S(2)
                    lssys%line_ex(gg_tplines(jj),gg_rows(2)) = S(1)
                    lssys%line_ey(gg_tplines(jj),gg_rows(2)) = S(2)

                enddo

                ! Lines in tp_ie changed. Dont make any further changes. 
                step_junc = .false.  
                skip_tp   = .true.    
                
            else
                ! Two neighbouring tp elms. All three active ls has a single line segment in element.
            
                ! Determine which two points are on the same element edge                
                ex_e = mesh%ex(:,tp_ie); ey_e = mesh%ey(:,tp_ie)
                minx = minval(ex_e)
                maxx = maxval(ex_e)
                miny = minval(ey_e)
                maxy = maxval(ey_e)

                ! Determine unique connecting points
                npoints = 6

                ! Points on the same edge
                nedges = 4
                allocate(pe(nedges,npoints),stat=ierr)
                pe = .false.
                pe(1,unique_points_idx) = abs(points(unique_points_idx,1)-minx).lt.1d-12 ! left
                pe(2,unique_points_idx) = abs(points(unique_points_idx,1)-maxx).lt.1d-12 ! right
                pe(3,unique_points_idx) = abs(points(unique_points_idx,2)-miny).lt.1d-12 ! lower
                pe(4,unique_points_idx) = abs(points(unique_points_idx,2)-maxy).lt.1d-12 ! upper
                incompatible_edges      = pack(idx4,count(pe,2).eq.2)                

                if (size(incompatible_edges).eq.0) then
                    skip_tp = .true.
                else

                    do jj=1,size(incompatible_edges)
                        
                        incompatible_edge = incompatible_edges(jj);
        
                        ! Find duplicates of points along incompatible edge
                        duplicate_points_idx = pack(idx6,pe(incompatible_edge,:))
                                        
                        ! Determine new compatible point
                        compatible_point = 0d0
                        do i = 1, 6
                            if (pe(incompatible_edge,i)) then
                                compatible_point = compatible_point + points(i, :)
                            end if
                        end do
                        compatible_point = compatible_point/2d0

                        ! Update points
                        do i=1,size(duplicate_points_idx)
                            points(duplicate_points_idx(i),:) = compatible_point
                        enddo

                        ! Update connecting line points
                        lsa_connA = points(1,:)
                        lsa_connB = points(2,:)
                        lsb_connA = points(3,:)
                        lsb_connB = points(4,:)
                        lsc_connA = points(5,:)
                        lsc_connB = points(6,:)

                        where (lsa_log_conn(1,:)) line_exa(lsa_conn_lines(1),:) = lsa_connA(1)
                        where (lsa_log_conn(1,:)) line_eya(lsa_conn_lines(1),:) = lsa_connA(2)
                        where (lsa_log_conn(2,:)) line_exa(lsa_conn_lines(2),:) = lsa_connB(1)
                        where (lsa_log_conn(2,:)) line_eya(lsa_conn_lines(2),:) = lsa_connB(2)

                        where (lsb_log_conn(1,:)) line_exb(lsb_conn_lines(1),:) = lsb_connA(1)
                        where (lsb_log_conn(1,:)) line_eyb(lsb_conn_lines(1),:) = lsb_connA(2)
                        where (lsb_log_conn(2,:)) line_exb(lsb_conn_lines(2),:) = lsb_connB(1)
                        where (lsb_log_conn(2,:)) line_eyb(lsb_conn_lines(2),:) = lsb_connB(2)

                        where (lsc_log_conn(1,:)) line_exc(lsc_conn_lines(1),:) = lsc_connA(1)
                        where (lsc_log_conn(1,:)) line_eyc(lsc_conn_lines(1),:) = lsc_connA(2)
                        where (lsc_log_conn(2,:)) line_exc(lsc_conn_lines(2),:) = lsc_connB(1)
                        where (lsc_log_conn(2,:)) line_eyc(lsc_conn_lines(2),:) = lsc_connB(2)
                        

                        ! Determine which of the connecting points are unique
                        call get_unique_idx(unique_points_idx, points)
                        nunique = size(unique_points_idx)                        
                        A = points(unique_points_idx(1),:)
                        B = points(unique_points_idx(2),:)
                        C = points(unique_points_idx(3),:)

                        if (nunique>3) then
                            print *, 'Error. Not finding unique points for neighb. tp'
                        endif

                        ! If tp at edge of domain, dont enter stepping tp loop
                        if (size(lsa_conn_lines).eq.1 .or. size(lsb_conn_lines).eq.1 .or. size(lsc_conn_lines).eq.1) then
                            F         = (A + B + C)/3d0
                            step_junc = .false.
                        else
                            ! Compute isogonic point
                            call get_isogonic(F,step_junc,A,B,C)
                        endif

                    enddo
                endif
            endif
        endif


        ! --- Enter stepping loop if F not found (internal angle > 120 degrees)

        Aprev = 0d0
        Bprev = 0d0
        Cprev = 0d0
        steps = 1
        do while (step_junc)

            ! Find next point along line in lsa
            call get_connecting_lines(line_elm_idx, tp_ie, lssys%int_elms, lsa_rm_lines, mesh, line_exa, line_eya, &
            nsega, lsa_conn_lines, lsa_connA, lsa_iconnA, lsa_connB, lsa_iconnB)

            ! Find next point along line in lsb
            call get_connecting_lines(line_elm_idx, tp_ie, lssys%int_elms, lsb_rm_lines, mesh, line_exb, line_eyb, &
            nsegb, lsb_conn_lines, lsb_connA, lsb_iconnA, lsb_connB, lsb_iconnB)

            ! Find next point along line in lsc
            call get_connecting_lines(line_elm_idx, tp_ie, lssys%int_elms, lsc_rm_lines, mesh, line_exc, line_eyc, &
            nsegc, lsc_conn_lines, lsc_connA, lsc_iconnA, lsc_connB, lsc_iconnB)


            ! Update lines to be removed
            rm_idx = rm_idx + 1
            lssys%rm_lines(rm_idx*2+1:rm_idx*2+2,ga)     = lsa_rm_lines
            lssys%rm_lines(rm_idx*2+1:rm_idx*2+2,gb)     = lsb_rm_lines
            lssys%rm_lines(rm_idx*2+1:rm_idx*2+2,gc)     = lsc_rm_lines

            lssys%rm_lines_tmp(rm_idx*2+1:rm_idx*2+2,ga) = lsa_rm_lines
            lssys%rm_lines_tmp(rm_idx*2+1:rm_idx*2+2,gb) = lsb_rm_lines
            lssys%rm_lines_tmp(rm_idx*2+1:rm_idx*2+2,gc) = lsc_rm_lines


            ! Determine unique connecting points
            points(1,:) = lsa_connA
            points(2,:) = lsa_connB
            points(3,:) = lsb_connA
            points(4,:) = lsb_connB
            points(5,:) = lsc_connA
            points(6,:) = lsc_connB
            
            
            call get_unique_idx(unique_points_idx, points)
            nunique = size(unique_points_idx)
                
            ! Find unique points A, B, and C
            if (nunique.eq.3) then
                A = points(unique_points_idx(1),:)
                B = points(unique_points_idx(2),:)
                C = points(unique_points_idx(3),:)

                
                if (short_tpline) then
                    step_junc    = .false.
                    short_tpline = .false.
                else
                    if ((norm2(A-Aprev).lt.1d-14 .and. norm2(B-Bprev).lt.1d-14 .and. norm2(C-Cprev).lt.1d-14) .or. steps.ge.1) then
                        ! Same unique points or stepp_max reached in tp elm
                        F         = (A + B + C)/3d0
                        step_junc = .false.                    
                    else
                        ! Compute isogonic point
                        call get_isogonic(F,step_junc,A,B,C)
                        Aprev = A
                        Bprev = B
                        Cprev = C
                    endif
                endif

            else
                print *, 'Error. More than three connecting points inside step tp loop'
            endif

            steps = steps + 1


        enddo

        ! --- Done with stepping. Isogonic point found or tpelm skipped ---

        if (.not. skip_tp) then
            ! OBS DONT zero out lines if skip_tp (all lines already adjusted in this case)
            
            ! Zero out rm_lines in int_elms
            do g=1,lssys%ngrains                                
                lssys%int_elms(lssys%line_elms(pack(lssys%rm_lines_tmp(:,g),lssys%rm_lines_tmp(:,g).gt.0),g),g) = 0
            enddo
            lssys%rm_lines_tmp = 0

            ! Zero out old lines in triple junction element
            rm_idx = rm_idx + 1
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(ga),ga) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,ga) = ival(1)            
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gb),gb) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gb) = ival(1)            
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gc),gc) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gc) = ival(1)
            
            ! Add isogonic point to level set functions
            call add_isogonic(lssys%line_ex(:,g_rowsa(1):g_rowsa(2)),lssys%line_ey(:,g_rowsa(1):g_rowsa(2)),lssys%tplines(:,ga), &
            lsa_connA,lsa_connB,F,nsega,lssys%int_elms(:,ga),mesh)

            ! Add isogonic point to level set functions
            call add_isogonic(lssys%line_ex(:,g_rowsb(1):g_rowsb(2)),lssys%line_ey(:,g_rowsb(1):g_rowsb(2)),lssys%tplines(:,gb), &
            lsb_connA,lsb_connB,F,nsegb,lssys%int_elms(:,gb),mesh)

            ! Add isogonic point to level set functions
            call add_isogonic(lssys%line_ex(:,g_rowsc(1):g_rowsc(2)),lssys%line_ey(:,g_rowsc(1):g_rowsc(2)),lssys%tplines(:,gc), &
            lsc_connA,lsc_connB,F,nsegc,lssys%int_elms(:,gc),mesh)


            ! Add elm containing isogonic point to int_elms_gr
            do ie = 1, mesh%nelm
                ! Check if point F is inside the current element
                call point_in_rectangle(is_inside, F(1), F(2), mesh%ex(:,ie), mesh%ey(:,ie))
                if (is_inside) then
                    ! Stop iterating once a containing element is found
                    exit
                endif
            enddo

            ! Update line_seg counter
            lssys%line_seg(ga) = nsega + 2
            lssys%line_seg(gb) = nsegb + 2
            lssys%line_seg(gc) = nsegc + 2
        
            ! Update rm idx
            rm_idx = rm_idx + 1
                
        endif

        ! Zero out rm lines in line_ex and line_ey
        do g=1,lssys%ngrains
            g_rows = [2*(g-1)+1, 2*(g-1)+2]
            lssys%line_ex(pack(lssys%rm_lines(:,g),lssys%rm_lines(:,g).gt.0),g_rows) = 0
            lssys%line_ey(pack(lssys%rm_lines(:,g),lssys%rm_lines(:,g).gt.0),g_rows) = 0
        enddo

    enddo


    ! Reorder lines such that only those not listed in rm_lines kept
    do g=1,lssys%ngrains        
        g_rows                                         = [2*(g-1)+1, 2*(g-1)+2];
        ! call setdiff(keep_lines, [(i, i = 1, lssys%line_seg(g))], pack(lssys%rm_lines(1:lssys%line_seg(g),g), &
        ! lssys%rm_lines(1:lssys%line_seg(g),g).gt.0))
        call setdiff(keep_lines, [(i, i = 1, lssys%line_seg(g))], pack(lssys%rm_lines(:,g), &
        lssys%rm_lines(:,g).gt.0))      
        nkeep                                          = size(keep_lines)
        lssys%line_ex(1:nkeep,g_rows)                  = lssys%line_ex(keep_lines,g_rows)
        lssys%line_ex(nkeep+1:lssys%nseg_alloc,g_rows) = 0
        lssys%line_ey(1:nkeep,g_rows)                  = lssys%line_ey(keep_lines,g_rows)
        lssys%line_ey(nkeep+1:lssys%nseg_alloc,g_rows) = 0
        lssys%tplines(1:nkeep,g)                       = lssys%tplines(keep_lines,g)
        lssys%tplines(nkeep+1:lssys%nseg_alloc,g)      = .false.
        lssys%line_seg(g)                              = nkeep
    enddo    

    return
end subroutine tp_reconstruction

subroutine tp_reconstruction_spatial(lssys,mesh)
    ! --- Routine for fixing triple junctions s.t they are compatible (check elements containing three ls isocontours) ---
    implicit none

    ! Intent inout
    type(ls_system), intent(inout) :: lssys

    ! Intent in
    type(mesh_system), intent(in)  :: mesh

    ! Subroutine variables    
    integer, allocatable   :: tp_elms(:), fp_elms(:), tp_elms_reordered(:), tmp(:), ls_at_tp(:)
    integer                :: grain_idx(lssys%ngrains), line_elm_idx(lssys%nseg_alloc), nls_at_tp
    integer                :: ntp_elms, tp_ie, ga, gb, gc, gd, g_rowsa(2), g_rowsb(2), g_rowsc(2), g_rowsd(2)
    integer                :: nfp_elms, rm_idx, junc
    integer                :: i, ierr
    real(dp)               :: ex_e(4), ey_e(4), points(6,2)
    real(dp)               :: points3(12,2), points4(4,2), points4b(4,2)
    logical                :: skip_tp

    ! Short 
    integer                :: nsega, nsegb, nsegc, nsegd, nunique
    integer, allocatable   :: lsa_conn_lines(:), lsb_conn_lines(:), lsc_conn_lines(:), lsd_conn_lines(:),unique_points_idx(:)
    integer, allocatable   :: gg_tplines(:), incompatible_edges(:), duplicate_points_idx(:), keep_lines(:)
    real(dp), allocatable  :: line_exa(:,:), line_eya(:,:), line_exb(:,:), line_eyb(:,:), line_exc(:,:), line_eyc(:,:)
    real(dp), allocatable  :: line_exd(:,:), line_eyd(:,:)
    logical, allocatable   :: lsa_log_conn(:,:), lsb_log_conn(:,:), lsc_log_conn(:,:), lsd_log_conn(:,:), pe(:,:)
    real(dp), dimension(2) :: lsa_connA,lsa_iconnA,lsa_connB,lsa_iconnB,lsa_connC,lsa_iconnC,lsa_connD,lsa_iconnD
    real(dp), dimension(2) :: lsb_connA,lsb_iconnA,lsb_connB,lsb_iconnB,lsb_connC,lsb_iconnC,lsb_connD,lsb_iconnD
    real(dp), dimension(2) :: lsc_connA,lsc_iconnA,lsc_connB,lsc_iconnB,lsc_connC,lsc_iconnC,lsc_connD,lsc_iconnD
    real(dp), dimension(2) :: lsd_connA,lsd_iconnA,lsd_connB,lsd_iconnB,lsd_connC,lsd_iconnC,lsd_connD,lsd_iconnD
    real(dp), dimension(2) :: A, B, C, D, E, F, PP, S, compatible_point, Aprev, Bprev, Cprev, M
    real(dp)               :: minx,maxx,miny,maxy
    logical                :: step_junc, is_inside, fp_elm
    integer                :: idx3(3) = [1,2,3], gi(1), gidx3(3), gg, g_rows(2), gg_rows(2), gother2(2), p1(1), jj, g, minidx(1)
    integer                :: idx4(4) = [1,2,3,4], nedges, npoints, incompatible_edge, idx6(6)=[1,2,3,4,5,6], steps, ival(1), ie
    integer                :: nkeep, lsa_rm_lines(2), lsb_rm_lines(2), lsc_rm_lines(2)
    logical                :: short_tpline
        
    ! Zero out
    lssys%rm_lines     = 0
    lssys%rm_lines_tmp = 0
    lssys%tplines      = .false.

    ! Determine tp elements
    ntp_elms = count(count(lssys%int_elms,2).ge.3)    
    allocate(tp_elms(ntp_elms),stat=ierr)
    tp_elms = pack(mesh%elmidx,count(lssys%int_elms,2).ge.3)

    ! Check if any elements contain four ls and reorder tpelms s.t. these fp_elms are placed last in tp_elms
    nfp_elms = count(count(lssys%int_elms,2).eq.4)  
    
    ! Reorder tpelms
    if (nfp_elms .gt. 0) then
        allocate(fp_elms(nfp_elms),stat=ierr)
        fp_elms = pack(mesh%elmidx,count(lssys%int_elms,2).eq.4)
        allocate(tp_elms_reordered(ntp_elms),stat=ierr)
        call setdiff(tmp,tp_elms,fp_elms)
        tp_elms_reordered(1:ntp_elms-nfp_elms)          = tmp
        tp_elms_reordered(ntp_elms-nfp_elms+1:ntp_elms) = fp_elms
        tp_elms                                         = tp_elms_reordered
        print *, 'OBS 4 ls active in at least one element!'
    end if

    ! Initialize grain_idx and line_elm_idx arrays
    grain_idx    = [(i, i = 1, lssys%ngrains)]
    line_elm_idx = [(i, i = 1, lssys%nseg_alloc)]

    ! Allocate rm_lines
    rm_idx = 0

    do junc=1,ntp_elms

        ! Triple junction element
        tp_ie = tp_elms(junc)

        ! Identify grains to be matched at triple junction
        nls_at_tp    = count(lssys%int_elms(tp_ie, :))
        allocate(ls_at_tp(nls_at_tp),stat=ierr)  
        ls_at_tp     = pack(grain_idx, lssys%int_elms(tp_ie, :))
        ga           = ls_at_tp(1)
        gb           = ls_at_tp(2)
        gc           = ls_at_tp(3)
        g_rowsa      = [(2 * (ga - 1) + 1), (2 * (ga - 1) + 2)]
        g_rowsb      = [(2 * (gb - 1) + 1), (2 * (gb - 1) + 2)]
        g_rowsc      = [(2 * (gc - 1) + 1), (2 * (gc - 1) + 2)]
        skip_tp      = .false.
        short_tpline = .false.

        ! Check if tp_ie is a four point junction element
        fp_elm = count(lssys%int_elms(tp_ie,:)).eq.4

        ! Short form
        nsega     = lssys%line_seg(ga)
        nsegb     = lssys%line_seg(gb)
        nsegc     = lssys%line_seg(gc)
        allocate(line_exa(nsega,2),stat=ierr)
        allocate(line_eya(nsega,2),stat=ierr)
        allocate(line_exb(nsegb,2),stat=ierr)
        allocate(line_eyb(nsegb,2),stat=ierr)
        allocate(line_exc(nsegc,2),stat=ierr)
        allocate(line_eyc(nsegc,2),stat=ierr)
        line_exa  = lssys%line_ex(1:nsega,g_rowsa)
        line_eya  = lssys%line_ey(1:nsega,g_rowsa)
        line_exb  = lssys%line_ex(1:nsegb,g_rowsb)
        line_eyb  = lssys%line_ey(1:nsegb,g_rowsb)
        line_exc  = lssys%line_ex(1:nsegc,g_rowsc)
        line_eyc  = lssys%line_ey(1:nsegc,g_rowsc)

        ! Ls a
        call get_connecting_points_ie_spatial(mesh%newex, mesh%newey, tp_ie, line_elm_idx, line_exa, line_eya, nsega, &
        lssys%line_elms(:,ga), lsa_connA, lsa_iconnA, lsa_connB, lsa_iconnB, lsa_connC, lsa_iconnC, lsa_connD, lsa_iconnD, &
        lsa_conn_lines, lsa_log_conn)

        ! Ls b
        call get_connecting_points_ie_spatial(mesh%newex, mesh%newey, tp_ie, line_elm_idx, line_exb, line_eyb, nsegb, &
        lssys%line_elms(:,gb), lsb_connA, lsb_iconnA, lsb_connB, lsb_iconnB, lsb_connC, lsb_iconnC, lsb_connD, lsb_iconnD, &
        lsb_conn_lines, lsb_log_conn)

        ! Ls c
        call get_connecting_points_ie_spatial(mesh%newex, mesh%newey, tp_ie, line_elm_idx, line_exc, line_eyc, nsegc, &
        lssys%line_elms(:,gc), lsc_connA, lsc_iconnA, lsc_connB, lsc_iconnB, lsc_connC, lsc_iconnC, lsc_connD, lsc_iconnD, &
        lsc_conn_lines, lsc_log_conn)

        
        ! Determine unique connecting points
        points(1,:) = lsa_connA
        points(2,:) = lsa_connB
        points(3,:) = lsb_connA
        points(4,:) = lsb_connB
        points(5,:) = lsc_connA
        points(6,:) = lsc_connB
        
        
        call get_unique_idx(unique_points_idx, points)
        nunique = size(unique_points_idx)
        
        ! Find unique points A, B, and C
        if (nunique.eq.3) then
            A = points(unique_points_idx(1),:)
            B = points(unique_points_idx(2),:)
            C = points(unique_points_idx(3),:)

            
            if (size(lsa_conn_lines).eq.1 .or. size(lsb_conn_lines).eq.1 .or. size(lsc_conn_lines).eq.1) then
                ! Tp at edge of domain. Dont enter stepping loop
                F         = (A + B + C)/3d0
                step_junc = .false.
            else
                ! Compute isogonic point
                call get_isogonic(F,step_junc,A,B,C)

                ! if ((step_junc.eqv..false.) .and. ((norm2(A-F).lt.critical_length) .or. (norm2(B-F).lt.critical_length) .or. &
                ! (norm2(C-F).lt.critical_length))) then
                !     print *, 'Too short line in tp elm'
                !     step_junc    = .true.
                !     short_tpline = .true.
                ! endif

            endif
            
        elseif(fp_elm) then
            ! 4 point junction - take special care

            ! Dont consider as a triple junction element. And don't step.
            skip_tp   = .true.
            step_junc = .false.

            ! Identify grains to be matched at 4p junction
            gd        = ls_at_tp(4)
            g_rowsd   = [2*(gd-1) + 1,2*(gd-1) + 2]
            nsegd     = lssys%line_seg(gd)
            allocate(line_exc(nsega,2),stat=ierr)
            allocate(line_eyc(nsega,2),stat=ierr)
            line_exd  = lssys%line_ex(1:nsegd,g_rowsd)
            line_eyd  = lssys%line_ey(1:nsegd,g_rowsd)

            ! Ls d
            call get_connecting_points_ie_spatial(mesh%newex, mesh%newey, tp_ie, line_elm_idx, line_exd, line_eyd, nsegd, &
            lssys%line_elms(:,gd), lsd_connA, lsd_iconnA, lsd_connB, lsd_iconnB, lsd_connC, lsd_iconnC, lsd_connD, lsd_iconnD, &
            lsd_conn_lines,lsd_log_conn)

            ! Add lines in triple junction element to rm_lines
            rm_idx = rm_idx + 1
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(ga),ga) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,ga) = ival(1)
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gb),gb) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gb) = ival(1)
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gc),gc) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gc) = ival(1)
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gd),gd) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gd) = ival(1)

            ! Set isogonic point as mean of all input points
            F = (lsa_connA + lsa_connB + lsb_connA + lsb_connB + lsc_connA + lsc_connB + lsd_connA + lsd_connB)/8d0

            ! Add isogonic point to level set function a
            call add_isogonic(lssys%line_ex(:,g_rowsa(1):g_rowsa(2)),lssys%line_ey(:,g_rowsa(1):g_rowsa(2)), &
            lssys%tplines(:,ga), lsa_connA,lsa_connB,F,nsega,lssys%int_elms(:,ga),mesh)

            ! Add isogonic point to level set function b
            call add_isogonic(lssys%line_ex(:,g_rowsb(1):g_rowsb(2)),lssys%line_ey(:,g_rowsb(1):g_rowsb(2)), &
            lssys%tplines(:,gb), lsb_connA,lsb_connB,F,nsegb,lssys%int_elms(:,gb),mesh)

            ! Add isogonic point to level set function c
            call add_isogonic(lssys%line_ex(:,g_rowsc(1):g_rowsc(2)),lssys%line_ey(:,g_rowsc(1):g_rowsc(2)), &
            lssys%tplines(:,gc), lsc_connA,lsc_connB,F,nsegc,lssys%int_elms(:,gc),mesh)

            ! Add isogonic point to level set function d
            call add_isogonic(lssys%line_ex(:,g_rowsd(1):g_rowsd(2)),lssys%line_ey(:,g_rowsd(1):g_rowsd(2)), &
            lssys%tplines(:,gd), lsd_connA,lsd_connB,F,nsegd,lssys%int_elms(:,gd),mesh)

            ! Add isogonic to tppoints
            lssys%ntp_points = lssys%ntp_points  + 1
            lssys%tp_points(lssys%ntp_points ,:) = F

            ! Update line_seg counter
            lssys%line_seg(ga) = nsega + 2
            lssys%line_seg(gb) = nsegb + 2
            lssys%line_seg(gc) = nsegc + 2
            lssys%line_seg(gd) = nsegd + 2
        
            ! Update rm idx
            rm_idx = rm_idx + 1

        elseif (nunique.gt.3) then
            ! More than 3 connecting points

            if (size(lsa_conn_lines).eq.4 .or. size(lsb_conn_lines).eq.4 .or. size(lsc_conn_lines).eq.4) then                
                ! One ls has two line segments in the element. Don't consider as a triple junction.
                points3(1,:)  = lsa_connA
                points3(2,:)  = lsa_connB
                points3(3,:)  = lsa_connC
                points3(4,:)  = lsa_connD
                points3(5,:)  = lsb_connA
                points3(6,:)  = lsb_connB
                points3(7,:)  = lsb_connC
                points3(8,:)  = lsb_connD
                points3(9,:)  = lsc_connA
                points3(10,:) = lsc_connB
                points3(11,:) = lsc_connC
                points3(12,:) = lsc_connD

                ! Determine what level set (assume it can only be one) that has two active lines at the element
                gi = pack(idx3, [size(lsa_conn_lines).eq.4, size(lsb_conn_lines).eq.4, size(lsc_conn_lines).eq.4])
                
                ! Connecting points
                A = points3((gi(1)-1)*4+1,:)
                B = points3((gi(1)-1)*4+2,:)
                C = points3((gi(1)-1)*4+3,:)
                D = points3((gi(1)-1)*4+4,:)

                points4(1,:) = A
                points4(2,:) = B
                points4(3,:) = C
                points4(4,:) = D

                ! Determine gg (level set with two line segments in tp_ie)
                gidx3    = [ga, gb, gc]
                gg       = gidx3(gi(1))
                gg_rows  = [2*(gg-1) + 1,2*(gg-1) + 2]

                ! Find other 2 grains at tp_ie
                gother2  = pack(gidx3,.not.[size(lsa_conn_lines).eq.4, size(lsb_conn_lines).eq.4, size(lsc_conn_lines).eq.4])

                ! Find which line segments of gg is in tp_ie
                gg_tplines = pack(line_elm_idx(1:lssys%line_seg(gg)),lssys%line_elms(1:lssys%line_seg(gg),gg).eq.tp_ie)
                

                ! Change line segments in gother(1) and gg to connecting points
                do jj=1,2
                    g        = gother2(jj)
                    g_rows   = [2*(g-1) + 1,2*(g-1) + 2]; 
                    p1       = pack(line_elm_idx(1:lssys%line_seg(g)),lssys%line_elms(1:lssys%line_seg(g),g).eq.tp_ie) ! Size 1
                    E        = [lssys%line_ex(p1(1),g_rows(1)), lssys%line_ey(p1(1),g_rows(1))]
                    F        = [lssys%line_ex(p1(1),g_rows(2)), lssys%line_ey(p1(1),g_rows(2))]                                
                    points4b = points4

                    ! Determine what point E is closest to
                    PP = E
                    minidx = minloc([norm2(PP-points4b(1,:)), norm2(PP-points4b(2,:)), &
                                     norm2(PP-points4b(3,:)), norm2(PP-points4b(4,:))])                                    
                    S = points4b(minidx(1),:)
                    points4b(minidx,:) = 1d8

                    ! Change lines
                    lssys%line_ex(p1(1),g_rows(1))           = S(1)
                    lssys%line_ey(p1(1),g_rows(1))           = S(2)
                    lssys%line_ex(gg_tplines(jj),gg_rows(1)) = S(1)
                    lssys%line_ey(gg_tplines(jj),gg_rows(1)) = S(2)
            
                    ! Determine what point E is closest to
                    PP = F
                    minidx = minloc([norm2(PP-points4b(1,:)), norm2(PP-points4b(2,:)), &
                                     norm2(PP-points4b(3,:)), norm2(PP-points4b(4,:))])                                    
                    S = points4b(minidx(1),:)

                    ! Change lines
                    lssys%line_ex(p1(1),g_rows(2))           = S(1)
                    lssys%line_ey(p1(1),g_rows(2))           = S(2)
                    lssys%line_ex(gg_tplines(jj),gg_rows(2)) = S(1)
                    lssys%line_ey(gg_tplines(jj),gg_rows(2)) = S(2)

                enddo

                ! Lines in tp_ie changed. Dont make any further changes. 
                step_junc = .false.  
                skip_tp   = .true.    
                
            else
                ! Two neighbouring tp elms. All three active ls has a single line segment in element.
            
                ! FIX
                print *, 'Two neighbouring tp elms - Not fixed in spatial setting'

                ! Determine which two points are on the same element edge                
                ex_e = mesh%newex(:,tp_ie); ey_e = mesh%newey(:,tp_ie)
                minx = minval(ex_e)
                maxx = maxval(ex_e)
                miny = minval(ey_e)
                maxy = maxval(ey_e)

                ! Determine unique connecting points
                npoints = 6

                ! Points on the same edge
                nedges = 4
                allocate(pe(nedges,npoints),stat=ierr)
                pe = .false.
                pe(1,unique_points_idx) = abs(points(unique_points_idx,1)-minx).lt.1d-12 ! left
                pe(2,unique_points_idx) = abs(points(unique_points_idx,1)-maxx).lt.1d-12 ! right
                pe(3,unique_points_idx) = abs(points(unique_points_idx,2)-miny).lt.1d-12 ! lower
                pe(4,unique_points_idx) = abs(points(unique_points_idx,2)-maxy).lt.1d-12 ! upper
                incompatible_edges      = pack(idx4,count(pe,2).eq.2)

                if (size(incompatible_edges).eq.0) then
                    skip_tp = .true.
                else

                    do jj=1,size(incompatible_edges)
                        
                        incompatible_edge = incompatible_edges(jj);
        
                        ! Find duplicates of points along incompatible edge
                        duplicate_points_idx = pack(idx6,pe(incompatible_edge,:))
                                        
                        ! Determine new compatible point
                        compatible_point = 0d0
                        do i = 1, 6
                            if (pe(incompatible_edge,i)) then
                                compatible_point = compatible_point + points(i, :)
                            end if
                        end do
                        compatible_point = compatible_point/2d0

                        ! Update points
                        do i=1,size(duplicate_points_idx)
                            points(duplicate_points_idx(i),:) = compatible_point
                        enddo

                        ! Update connecting line points
                        lsa_connA = points(1,:)
                        lsa_connB = points(2,:)
                        lsb_connA = points(3,:)
                        lsb_connB = points(4,:)
                        lsc_connA = points(5,:)
                        lsc_connB = points(6,:)

                        where (lsa_log_conn(1,:)) line_exa(lsa_conn_lines(1),:) = lsa_connA(1)
                        where (lsa_log_conn(1,:)) line_eya(lsa_conn_lines(1),:) = lsa_connA(2)
                        where (lsa_log_conn(2,:)) line_exa(lsa_conn_lines(2),:) = lsa_connB(1)
                        where (lsa_log_conn(2,:)) line_eya(lsa_conn_lines(2),:) = lsa_connB(2)

                        where (lsb_log_conn(1,:)) line_exb(lsb_conn_lines(1),:) = lsb_connA(1)
                        where (lsb_log_conn(1,:)) line_eyb(lsb_conn_lines(1),:) = lsb_connA(2)
                        where (lsb_log_conn(2,:)) line_exb(lsb_conn_lines(2),:) = lsb_connB(1)
                        where (lsb_log_conn(2,:)) line_eyb(lsb_conn_lines(2),:) = lsb_connB(2)

                        where (lsc_log_conn(1,:)) line_exc(lsc_conn_lines(1),:) = lsc_connA(1)
                        where (lsc_log_conn(1,:)) line_eyc(lsc_conn_lines(1),:) = lsc_connA(2)
                        where (lsc_log_conn(2,:)) line_exc(lsc_conn_lines(2),:) = lsc_connB(1)
                        where (lsc_log_conn(2,:)) line_eyc(lsc_conn_lines(2),:) = lsc_connB(2)
                        

                        ! Determine which of the connecting points are unique
                        call get_unique_idx(unique_points_idx, points)
                        nunique = size(unique_points_idx)
                        A = points(unique_points_idx(1),:)
                        B = points(unique_points_idx(2),:)
                        C = points(unique_points_idx(3),:)

                        if (nunique>3) then
                            print *, 'Error. Not finding unique points for neighb. tp'
                        endif

                        ! If tp at edge of domain, dont enter stepping tp loop
                        if (size(lsa_conn_lines).eq.1 .or. size(lsb_conn_lines).eq.1 .or. size(lsc_conn_lines).eq.1) then
                            F         = (A + B + C)/3d0
                            step_junc = .false.
                        else
                            ! Compute isogonic point
                            call get_isogonic(F,step_junc,A,B,C)
                        endif

                    enddo
                endif
            endif
        endif


        ! --- Enter stepping loop if F not found (internal angle > 120 degrees)

        Aprev = 0d0
        Bprev = 0d0
        Cprev = 0d0
        steps = 1
        do while (step_junc)

            ! Find next point along line in lsa
            call get_connecting_lines(line_elm_idx, tp_ie, lssys%int_elms, lsa_rm_lines, mesh, line_exa, line_eya, &
            nsega, lsa_conn_lines, lsa_connA, lsa_iconnA, lsa_connB, lsa_iconnB)

            ! Find next point along line in lsb
            call get_connecting_lines(line_elm_idx, tp_ie, lssys%int_elms, lsb_rm_lines, mesh, line_exb, line_eyb, &
            nsegb, lsb_conn_lines, lsb_connA, lsb_iconnA, lsb_connB, lsb_iconnB)

            ! Find next point along line in lsc
            call get_connecting_lines(line_elm_idx, tp_ie, lssys%int_elms, lsc_rm_lines, mesh, line_exc, line_eyc, &
            nsegc, lsc_conn_lines, lsc_connA, lsc_iconnA, lsc_connB, lsc_iconnB)


            ! Update lines to be removed
            rm_idx = rm_idx + 1
            lssys%rm_lines(rm_idx*2+1:rm_idx*2+2,ga)     = lsa_rm_lines
            lssys%rm_lines(rm_idx*2+1:rm_idx*2+2,gb)     = lsb_rm_lines
            lssys%rm_lines(rm_idx*2+1:rm_idx*2+2,gc)     = lsc_rm_lines

            lssys%rm_lines_tmp(rm_idx*2+1:rm_idx*2+2,ga) = lsa_rm_lines
            lssys%rm_lines_tmp(rm_idx*2+1:rm_idx*2+2,gb) = lsb_rm_lines
            lssys%rm_lines_tmp(rm_idx*2+1:rm_idx*2+2,gc) = lsc_rm_lines


            ! Determine unique connecting points
            points(1,:) = lsa_connA
            points(2,:) = lsa_connB
            points(3,:) = lsb_connA
            points(4,:) = lsb_connB
            points(5,:) = lsc_connA
            points(6,:) = lsc_connB
            
            
            call get_unique_idx(unique_points_idx, points)
            nunique = size(unique_points_idx)
                
            ! Find unique points A, B, and C
            if (nunique.eq.3) then
                A = points(unique_points_idx(1),:)
                B = points(unique_points_idx(2),:)
                C = points(unique_points_idx(3),:)

                
                if (short_tpline) then
                    step_junc    = .false.
                    short_tpline = .false.
                else
                    if ((norm2(A-Aprev).lt.1d-14 .and. norm2(B-Bprev).lt.1d-14 .and. norm2(C-Cprev).lt.1d-14) .or. steps.ge.1) then
                        ! Same unique points or stepp_max reached in tp elm
                        F         = (A + B + C)/3d0
                        step_junc = .false.                    
                    else
                        ! Compute isogonic point
                        call get_isogonic(F,step_junc,A,B,C)
                        Aprev = A
                        Bprev = B
                        Cprev = C
                    endif
                endif

            else
                print *, 'Error. More than three connecting points inside step tp loop'
            endif

            steps = steps + 1


        enddo

        ! --- Done with stepping. Isogonic point found or tpelm skipped ---

        if (.not. skip_tp) then
            ! OBS DONT zero out lines if skip_tp (all lines already adjusted in this case)
            
            ! Zero out rm_lines in int_elms
            do g=1,lssys%ngrains                                
                lssys%int_elms(lssys%line_elms(pack(lssys%rm_lines_tmp(:,g),lssys%rm_lines_tmp(:,g).gt.0),g),g) = 0
            enddo
            lssys%rm_lines_tmp = 0

            ! Zero out old lines in triple junction element
            rm_idx = rm_idx + 1
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(ga),ga) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,ga) = ival(1)            
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gb),gb) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gb) = ival(1)            
            ival = pack(line_elm_idx, lssys%line_elms(1:lssys%line_seg(gc),gc) - tp_ie .eq. 0)
            lssys%rm_lines(rm_idx*2+1,gc) = ival(1)
            
            ! Add isogonic point to level set functions
            call add_isogonic(lssys%line_ex(:,g_rowsa(1):g_rowsa(2)),lssys%line_ey(:,g_rowsa(1):g_rowsa(2)), & 
            lssys%tplines(:,ga), lsa_connA,lsa_connB,F,nsega,lssys%int_elms(:,ga),mesh)

            ! Add isogonic point to level set functions
            call add_isogonic(lssys%line_ex(:,g_rowsb(1):g_rowsb(2)),lssys%line_ey(:,g_rowsb(1):g_rowsb(2)), & 
            lssys%tplines(:,gb), lsb_connA,lsb_connB,F,nsegb,lssys%int_elms(:,gb),mesh)

            ! Add isogonic point to level set functions
            call add_isogonic(lssys%line_ex(:,g_rowsc(1):g_rowsc(2)),lssys%line_ey(:,g_rowsc(1):g_rowsc(2)), &
            lssys%tplines(:,gc), lsc_connA,lsc_connB,F,nsegc,lssys%int_elms(:,gc),mesh)

            ! Add isogonic to tppoints
            lssys%ntp_points = lssys%ntp_points  + 1
            lssys%tp_points(lssys%ntp_points ,:) = F

            ! Add elm containing isogonic point to int_elms_gr
            do ie = 1, mesh%nelm
                ! Check if point F is inside the current element
                call point_in_rectangle(is_inside, F(1), F(2), mesh%newex(:,ie), mesh%newey(:,ie))
                if (is_inside) then
                    ! Stop iterating once a containing element is found
                    exit
                endif
            enddo

            ! Update line_seg counter
            lssys%line_seg(ga) = nsega + 2
            lssys%line_seg(gb) = nsegb + 2
            lssys%line_seg(gc) = nsegc + 2
        
            ! Update rm idx
            rm_idx = rm_idx + 1
                
        endif

        ! Zero out rm lines in line_ex and line_ey
        do g=1,lssys%ngrains
            g_rows = [2*(g-1)+1, 2*(g-1)+2]
            lssys%line_ex(pack(lssys%rm_lines(:,g),lssys%rm_lines(:,g).gt.0),g_rows) = 0
            lssys%line_ey(pack(lssys%rm_lines(:,g),lssys%rm_lines(:,g).gt.0),g_rows) = 0
        enddo

    enddo


    ! Reorder lines such that only those not listed in rm_lines kept
    do g=1,lssys%ngrains        
        g_rows                                         = [2*(g-1)+1, 2*(g-1)+2]
        ! call setdiff(keep_lines, [(i, i = 1, lssys%line_seg(g))], pack(lssys%rm_lines(1:lssys%line_seg(g),g), &
        ! lssys%rm_lines(1:lssys%line_seg(g),g).gt.0))
        call setdiff(keep_lines, [(i, i = 1, lssys%line_seg(g))], pack(lssys%rm_lines(:,g), &
        lssys%rm_lines(:,g).gt.0))
        nkeep                                            = size(keep_lines)
        lssys%line_ex(1:nkeep,g_rows)                  = lssys%line_ex(keep_lines,g_rows)
        lssys%line_ex(nkeep+1:lssys%nseg_alloc,g_rows) = 0
        lssys%line_ey(1:nkeep,g_rows)                  = lssys%line_ey(keep_lines,g_rows)
        lssys%line_ey(nkeep+1:lssys%nseg_alloc,g_rows) = 0
        lssys%tplines(1:nkeep,g)                       = lssys%tplines(keep_lines,g)
        lssys%tplines(nkeep+1:lssys%nseg_alloc,g)      = .false.
        lssys%line_seg(g)                              = nkeep
    enddo    

    return
end subroutine tp_reconstruction_spatial


subroutine add_isogonic(line_ex,line_ey,tplines,connA,connB,F,nseg,int_elms,mesh)
    ! --- Adding isogonic point to lines --- 
    implicit none

    ! Intent inout
    real(dp), intent(inout)       :: line_ex(:,:), line_ey(:,:)
    logical,intent(inout)         :: tplines(:)
    logical, intent(inout)        :: int_elms(:)
    
    ! Intent in
    real(dp), intent(in)          :: connA(2), connB(2), F(2)   
    integer, intent(in)           :: nseg    
    type(mesh_system), intent(in) :: mesh

    ! Subroutine variables
    integer  :: k,ntest=10, j, ie
    real(dp) :: A(2), B(2), v(2), xvec(10), yvec(10), M(2)
    logical :: is_inside
    

    ! Add line segemnt to go from connA to F
    line_ex(nseg+1,:) = [connA(1), F(1)]
    line_ey(nseg+1,:) = [connA(2), F(2)]

    ! Add extra line segemnt to go from F to connB
    line_ex(nseg+2,:) = [F(1), connB(1)]
    line_ey(nseg+2,:) = [F(2), connB(2)]

    ! Mark which lines are triple junction lines
    tplines(nseg+1)   = .true.
    tplines(nseg+2)   = .true.


    ! Update bool int_elms
    do k = 1, 2
        A      = [line_ex(nseg+k,1), line_ey(nseg+k,1)]
        B      = [line_ex(nseg+k,2), line_ey(nseg+k,2)]
        v      = A - B    
        call linspace_real(A(1), B(1), ntest, xvec)
        yvec   = v(2)/v(1)*xvec + A(2) - A(1)*v(2)/v(1)
        do j=1,ntest
            M   = [xvec(j),yvec(j)]
            do ie = 1,mesh%nelm
                ! Check if point F is inside the current element                  
                call point_in_rectangle(is_inside, M(1), M(2), mesh%ex(:,ie), mesh%ey(:,ie))
                if (is_inside) then
                    ! Stop iterating once a containing element is found
                    exit
                endif
            enddo
        
            ! Update int_elms
            int_elms(ie) = .true.
        enddo
    enddo


    return
end subroutine add_isogonic


subroutine get_connecting_lines(line_elm_idx, tp_ie, int_elms, rm_lines, mesh, line_ex, line_ey, nseg, conn_lines, &
    connA, iconnA, connB, iconnB)
    ! --- Subroutine for computing next connecting points along line segment ---
    implicit none

    ! Intent out
    integer, dimension(2), intent(inout) :: conn_lines, rm_lines
    real(dp), dimension(2)               :: connA, iconnA, connB, iconnB

    ! Intent in
    integer, intent(in)                  :: line_elm_idx(:), tp_ie, nseg
    real(dp), intent(in)                 :: line_ex(:,:), line_ey(:,:)
    logical, intent(in)                  :: int_elms(:,:)
    type(mesh_system), intent(in)        :: mesh

    ! Subroutine variables
    real(dp), dimension(2)               :: connA_in, iconnA_in, connB_in, iconnB_in, A1, A2
    integer                              :: k, ie, njunctions, iconnA_line, iconnB_line, rm_lineA, rm_lineB
    integer, dimension(:), allocatable   :: current_tpelms, other_tp_elms, iconnA_lines, iconnB_lines, iconnA_line1, iconnB_line1
    logical                              :: is_inside

    ! Define input
    connA_in  = connA
    iconnA_in = iconnA
    connB_in  = connB
    iconnB_in = iconnB

    ! Define rm lines as the connecting lines
    rm_lineA = conn_lines(1)
    rm_lineB = conn_lines(2)

    ! Define current_tpelms
    current_tpelms = pack(mesh%elmidx, count(int_elms, dim=2) >= 3)

    ! Define other_tp_elms
    call setdiff(other_tp_elms, current_tpelms, [tp_ie])

    ! Get the number of other triple junctions
    njunctions = size(other_tp_elms)

    ! Loop through other junction elements
    do k = 1, njunctions
        ie = other_tp_elms(k)
        
        ! Check if iconnA_in is present in other triple junction element ie
        call point_in_rectangle(is_inside, iconnA_in(1), iconnA_in(2), mesh%ex(:,ie), mesh%ey(:,ie))
        if (is_inside) then
            iconnA_line = 0
            rm_lineA    = 0
            connA       = connA_in
            iconnA      = iconnA_in
        end if

        ! Check if iconnB_in is present in other triple junction element ie
        call point_in_rectangle(is_inside, iconnB_in(1), iconnB_in(2), mesh%ex(:,ie), mesh%ey(:,ie))
        if (is_inside) then
            iconnB_line = 0
            rm_lineB    = 0
            connB       = connB_in
            iconnB      = iconnB_in
        end if

    end do

    ! Find next line segments connecting to iconnA
    if (rm_lineA .ne. 0) then
        iconnA_lines = pack(line_elm_idx(1:nseg), count(sqrt((line_ex-iconnA_in(1))**2 + (line_ey-iconnA_in(2))**2) .lt. &
        1.0d-13,2).eq.1)        
        if (size(iconnA_lines) .eq. 1) then
            ! There are no further lines connecting to iconnA -> Choose same points for connA and iconnA as before
            iconnA_line = iconnA_lines(1)
            connA       = connA_in
            iconnA      = iconnA_in
            rm_lineA    = 0
        else
            call setdiff(iconnA_line1, iconnA_lines, conn_lines)
            iconnA_line = iconnA_line1(1)
            ! Find new end point A
            A1 = [line_ex(iconnA_line,1), line_ey(iconnA_line,1)]
            A2 = [line_ex(iconnA_line,2), line_ey(iconnA_line,2)]

            if (norm2(A1-iconnA_in).lt.1d-14) then
                connA  = A1
                iconnA = A2
            else
                connA  = A2
                iconnA = A1
            endif
        end if
    end if

    ! Find next line segments connecting to iconnA
    if (rm_lineB .ne. 0) then
        iconnB_lines = pack(line_elm_idx(1:nseg), count(sqrt((line_ex-iconnB_in(1))**2 + (line_ey-iconnB_in(2))**2) .lt. &
        1.0d-13,2).eq.1)        
        if (size(iconnB_lines) .eq. 1) then
            ! There are no further lines connecting to iconnA -> Choose same points for connA and iconnA as before
            iconnB_line = iconnB_lines(1)
            connB       = connB_in
            iconnB      = iconnB_in
            rm_lineB    = 0
        else
            call setdiff(iconnB_line1, iconnB_lines, conn_lines)
            iconnB_line = iconnB_line1(1)
            ! Find new end point A
            A1 = [line_ex(iconnB_line,1), line_ey(iconnB_line,1)]
            A2 = [line_ex(iconnB_line,2), line_ey(iconnB_line,2)]

            if (norm2(A1-iconnB_in).lt.1d-14) then
                connB  = A1
                iconnB = A2
            else
                connB  = A2
                iconnB = A1
            endif
        end if
    end if

    ! Define new conn_lines
    conn_lines = [iconnA_line,iconnB_line]

    ! Define lines to be removed
    rm_lines = [rm_lineA, rm_lineB]

    ! Deallocate     
    deallocate(current_tpelms)
    deallocate(other_tp_elms)
    if (allocated(iconnA_lines)) then
        deallocate(iconnA_lines)
    endif
    if (allocated(iconnB_lines)) then
        deallocate(iconnB_lines)
    endif
    if (allocated(iconnA_line1)) then
        deallocate(iconnA_line1)
    endif
    if (allocated(iconnB_line1)) then
        deallocate(iconnB_line1)
    endif
    return
end subroutine get_connecting_lines


subroutine get_connecting_points_ie(ex, ey, tp_ie, line_elm_idx, line_ex, line_ey, nseg, line_elms, &
    connA, iconnA, connB, iconnB, connC, iconnC, connD, iconnD, conn_lines, log_conn)
    ! --- Routine for finding connecting lines and points to tp_ie ---
    
    implicit none

    ! Intent out
    real(dp), dimension(2), intent(out)  :: connA, iconnA, connB, iconnB, connC, iconnC, connD, iconnD
    integer, allocatable, intent(out)    :: conn_lines(:)
    logical, allocatable, intent(out)    :: log_conn(:,:)

    ! Intent in
    real(dp), dimension(:,:), intent(in) :: ex, ey, line_ex, line_ey
    integer, dimension(:), intent(in)    :: line_elm_idx
    integer, intent(in)                  :: nseg, tp_ie, line_elms(:)

    ! Subroutine variables
    real(dp), allocatable                :: xpoints(:,:), ypoints(:,:)
    real(dp)                             :: ex_e(4), ey_e(4), A(2), B(2)
    real(dp)                             :: minx, maxx, miny, maxy
    logical                              :: l1(size(line_ex,1),size(line_ex,2)), l2(size(line_ex,1),size(line_ex,2))
    logical                              :: l3(size(line_ex,1),size(line_ex,2)), l4(size(line_ex,1),size(line_ex,2))
    logical                              :: logg(size(line_elm_idx,1),2)
    integer                              :: nconn_lines, tp_line(1), ierr

    ! Find line segment points in ls a and ls b connecting to triple junction (ls c redundant)
    ex_e = ex(:,tp_ie)
    ey_e = ey(:,tp_ie)
    minx = minval(ex_e)
    maxx = maxval(ex_e)
    miny = minval(ey_e)
    maxy = maxval(ey_e)

    ! Find which line segments connect to the triple junction element
    l1 = .false.
    l2 = .false.
    l3 = .false.
    l4 = .false.
    
    where ((line_ex - minx) .ge. -1.0d-14) l1 = .true.
    where ((maxx - line_ex) .ge. -1.0d-14) l2 = .true.
    where ((line_ey - miny) .ge. -1.0d-14) l3 = .true.
    where ((maxy - line_ey) .ge. -1.0d-14) l4 = .true.
    logg = .false.
    logg(1:nseg,:) = l1 .and. l2 .and. l3 .and. l4    
    conn_lines  = pack(line_elm_idx, count(logg,2).eq.1)
    nconn_lines = size(conn_lines)


    allocate(xpoints(nconn_lines,2),stat=ierr)
    allocate(ypoints(nconn_lines,2),stat=ierr)
    xpoints  = line_ex(conn_lines,:)
    ypoints  = line_ey(conn_lines,:)
    allocate(log_conn(nconn_lines,2),stat=ierr)
    log_conn = logg(conn_lines,:)

    if (nconn_lines.eq.2) then
        connA  = [pack(xpoints(1,:),     log_conn(1,:)), pack(ypoints(1,:),     log_conn(1,:))]
        iconnA = [pack(xpoints(1,:),.not.log_conn(1,:)), pack(ypoints(1,:),.not.log_conn(1,:))]
        connB  = [pack(xpoints(2,:),     log_conn(2,:)), pack(ypoints(2,:),     log_conn(2,:))]
        iconnB = [pack(xpoints(2,:),.not.log_conn(2,:)), pack(ypoints(2,:),.not.log_conn(2,:))]
        connC  = 0d0
        iconnC = 0d0
        connD  = 0d0
        iconnD = 0d0
    elseif (nconn_lines.eq.1) then
        ! At domain border
        connA   = [pack(xpoints,log_conn), pack(ypoints,log_conn)]        
        tp_line = pack(line_elm_idx(1:nseg),line_elms(1:nseg).eq.tp_ie)
        A       = [line_ex((tp_line(1)),1), line_ey(tp_line(1),1)]
        B       = [line_ex(tp_line(1),2), line_ey(tp_line(1),2)]
        if (norm2(connA-A).gt.norm2(connA-B)) then
            connB = A
        else
            connB = B
        endif
        iconnA = 0d0
        iconnB = 0d0
        connC  = 0d0
        iconnC = 0d0
        connD  = 0d0
        iconnD = 0d0
    elseif (nconn_lines.eq.4) then
        ! 4 connecting lines to tp elm
        connA  = [pack(xpoints(1,:),     log_conn(1,:)), pack(ypoints(1,:),     log_conn(1,:))]
        iconnA = [pack(xpoints(1,:),.not.log_conn(1,:)), pack(ypoints(1,:),.not.log_conn(1,:))]
        connB  = [pack(xpoints(2,:),     log_conn(2,:)), pack(ypoints(2,:),     log_conn(2,:))]
        iconnB = [pack(xpoints(2,:),.not.log_conn(2,:)), pack(ypoints(2,:),.not.log_conn(2,:))]
        connC  = [pack(xpoints(3,:),     log_conn(3,:)), pack(ypoints(3,:),     log_conn(3,:))]
        iconnC = [pack(xpoints(3,:),.not.log_conn(3,:)), pack(ypoints(3,:),.not.log_conn(3,:))]
        connD  = [pack(xpoints(4,:),     log_conn(4,:)), pack(ypoints(4,:),     log_conn(4,:))]
        iconnD = [pack(xpoints(4,:),.not.log_conn(4,:)), pack(ypoints(4,:),.not.log_conn(4,:))]
    endif

    ! Deallocate
    deallocate(xpoints)
    deallocate(ypoints)

    return
end subroutine get_connecting_points_ie

subroutine get_connecting_points_ie_spatial(ex, ey, tp_ie, line_elm_idx, line_ex, line_ey, nseg, line_elms, &
    connA, iconnA, connB, iconnB, connC, iconnC, connD, iconnD, conn_lines, log_conn)
    ! --- Routine for finding connecting lines and points to tp_ie ---
    
    implicit none

    ! Intent out
    real(dp), dimension(2), intent(out)  :: connA, iconnA, connB, iconnB, connC, iconnC, connD, iconnD
    integer, allocatable, intent(out)    :: conn_lines(:)
    logical, allocatable, intent(out)    :: log_conn(:,:)

    ! Intent in
    real(dp), dimension(:,:), intent(in) :: ex, ey, line_ex, line_ey
    integer, dimension(:), intent(in)    :: line_elm_idx
    integer, intent(in)                  :: nseg, tp_ie, line_elms(:)

    ! Subroutine variables
    real(dp), allocatable                :: xpoints(:,:), ypoints(:,:)
    real(dp)                             :: ex_e(4), ey_e(4), A(2), B(2)
    logical                              :: logg(size(line_elm_idx,1),2)
    logical                              :: A_is_inside, B_is_inside
    integer                              :: nconn_lines, tp_line(1), ierr, iseg


    ! Node coordinates of tp elm
    ex_e = ex(:,tp_ie)
    ey_e = ey(:,tp_ie)

    ! Find line segment points in ls a and ls b connecting to triple junction (ls c redundant)
    logg = .false.
    do iseg = 1, nseg

        ! Line segment points A and B
        A =  [line_ex(iseg, 1), line_ey(iseg, 1)]
        B =  [line_ex(iseg, 2), line_ey(iseg, 2)]

        ! Check if point A is inside tp element (arbitrary qudarilateral)
        call point_in_arbitrary_rectangle(A_is_inside, A(1), A(2), ex_e, ey_e)

        ! Check if point B is inside tp element (arbitrary qudarilateral)
        call point_in_arbitrary_rectangle(B_is_inside, B(1), B(2), ex_e, ey_e)

        ! Assign to logg
        logg(iseg,:) = [A_is_inside, B_is_inside]

     enddo
    conn_lines  = pack(line_elm_idx, count(logg,2).eq.1)
    nconn_lines = size(conn_lines)


    allocate(xpoints(nconn_lines,2),stat=ierr)
    allocate(ypoints(nconn_lines,2),stat=ierr)
    xpoints  = line_ex(conn_lines,:)
    ypoints  = line_ey(conn_lines,:)
    allocate(log_conn(nconn_lines,2),stat=ierr)
    log_conn = logg(conn_lines,:)

    if (nconn_lines.eq.2) then
        connA  = [pack(xpoints(1,:),     log_conn(1,:)), pack(ypoints(1,:),     log_conn(1,:))]
        iconnA = [pack(xpoints(1,:),.not.log_conn(1,:)), pack(ypoints(1,:),.not.log_conn(1,:))]
        connB  = [pack(xpoints(2,:),     log_conn(2,:)), pack(ypoints(2,:),     log_conn(2,:))]
        iconnB = [pack(xpoints(2,:),.not.log_conn(2,:)), pack(ypoints(2,:),.not.log_conn(2,:))]
        connC  = 0d0
        iconnC = 0d0
        connD  = 0d0
        iconnD = 0d0
    elseif (nconn_lines.eq.1) then
        ! At domain border
        connA   = [pack(xpoints,log_conn), pack(ypoints,log_conn)]        
        tp_line = pack(line_elm_idx(1:nseg),line_elms(1:nseg).eq.tp_ie)
        A       = [line_ex((tp_line(1)),1), line_ey(tp_line(1),1)]
        B       = [line_ex(tp_line(1),2), line_ey(tp_line(1),2)]
        if (norm2(connA-A).gt.norm2(connA-B)) then
            connB = A
        else
            connB = B
        endif
        iconnA = 0d0
        iconnB = 0d0
        connC  = 0d0
        iconnC = 0d0
        connD  = 0d0
        iconnD = 0d0
    elseif (nconn_lines.eq.4) then
        ! 4 connecting lines to tp elm
        connA  = [pack(xpoints(1,:),     log_conn(1,:)), pack(ypoints(1,:),     log_conn(1,:))]
        iconnA = [pack(xpoints(1,:),.not.log_conn(1,:)), pack(ypoints(1,:),.not.log_conn(1,:))]
        connB  = [pack(xpoints(2,:),     log_conn(2,:)), pack(ypoints(2,:),     log_conn(2,:))]
        iconnB = [pack(xpoints(2,:),.not.log_conn(2,:)), pack(ypoints(2,:),.not.log_conn(2,:))]
        connC  = [pack(xpoints(3,:),     log_conn(3,:)), pack(ypoints(3,:),     log_conn(3,:))]
        iconnC = [pack(xpoints(3,:),.not.log_conn(3,:)), pack(ypoints(3,:),.not.log_conn(3,:))]
        connD  = [pack(xpoints(4,:),     log_conn(4,:)), pack(ypoints(4,:),     log_conn(4,:))]
        iconnD = [pack(xpoints(4,:),.not.log_conn(4,:)), pack(ypoints(4,:),.not.log_conn(4,:))]
    endif

    ! Deallocate
    deallocate(xpoints)
    deallocate(ypoints)

    return
end subroutine get_connecting_points_ie_spatial


subroutine setdiff(result, arr1, arr2)
    ! --- Outputs elements in arr1 not present in arr2 in results ---
    implicit none
  
    ! Intent out    
    integer, allocatable, intent(out) :: result(:)
  
    ! Intent in
    integer, intent(in)               :: arr1(:), arr2(:)

    ! Subroutine variables
    logical, allocatable              :: mask(:)
    integer                           :: i, ierr

    ! Mask
    allocate(mask(size(arr1)),stat=ierr)
    mask = .false.    
    do i=1,size(arr1)
        if (all(arr1(i).ne.arr2)) then            
            mask(i) = .true.
        endif
    enddo
    
    ! Output
    result = pack(arr1,mask)

    return
end subroutine setdiff
  

subroutine get_unique_idx(unique_points_idx_out,points)
    ! --- ---
    implicit none

    ! Intent out
    integer, allocatable :: unique_points_idx_out(:)

    ! Intent in
    real(dp), intent(in) :: points(:,:)

    ! Subrioutine variables
    integer              :: unique_points_idx(size(points,1))
    integer              :: i, j , nunique, ierr
    logical              :: is_unique

    ! Initialize unique_points_idx
    unique_points_idx = 0

    ! Loop through each point
    nunique = 0
    do i = 1,size(points, 1)
        is_unique = .true.
        
        ! Compare the current point with the previously checked points
        do j = 1,nunique
            if (norm2(points(i, :) - points(unique_points_idx(j), :)).lt.1d-14) then
                is_unique = .false.                
                exit ! No need to continue checking if a match is found
            endif
        enddo
        
        ! If the point is unique, add it to the unique_points array
        if (is_unique) then
            nunique = nunique + 1
            unique_points_idx(nunique) = i
        endif
    enddo

    ! unique_points_idx_out
    allocate(unique_points_idx_out(nunique),stat=ierr)
    unique_points_idx_out = unique_points_idx(1:nunique)
    
    return
end subroutine get_unique_idx


subroutine get_isogonic(F,step_junc,AA,BB,CC)
    ! --- Computing isogonic point for triangle AA-BB-CC. If any angle>120deg, ouput F=0 and step_junc=.true.

    ! Intent inout
    real(dp)               :: F(2)
    logical                :: step_junc

    ! Intent in
    real(dp), dimension(2) :: AA, BB, CC

    ! Subroutine variables
    real(dp)               :: a, b, c, Ad, Bd, Cd
    real(dp)               :: x, y, z
    real(dp), dimension(2) :: Avec, Bvec


    ! Side lengths
    a = norm2(BB - CC)
    b = norm2(AA - CC)
    c = norm2(AA - BB)

    if (a .eq. 0.0d0 .or. b .eq. 0.0d0 .or. c .eq. 0.0d0) then
        print *, 'Error in computing isogonic point'
        print *, 'One or more side lengths equal zero'
    end if

    ! Angles at vertices
    Ad = acosd((b**2 + c**2 - a**2) / (2.0d0 * b * c))
    Bd = acosd((a**2 + c**2 - b**2) / (2.0d0 * a * c))
    Cd = acosd((a**2 + b**2 - c**2) / (2.0d0 * a * b))

    if (Ad .ge. 120.0d0 .or. Bd .ge. 120.0d0 .or. Cd .ge. 120.0d0) then
        step_junc = .true.
        F         = [0.0d0, 0.0d0]
    else
        step_junc = .false.

        ! Triangular coordinates of isogonic point
        x = 1.0d0 / sind(Ad + 60.0d0)
        y = 1.0d0 / sind(Bd + 60.0d0)
        z = 1.0d0 / sind(Cd + 60.0d0)

        ! Directions given vertex C as origin
        Avec = AA - CC
        Bvec = BB - CC

        F(1) = CC(1) + a * x / (a * x + b * y + c * z) * Avec(1) + b * y / (a * x + b * y + c * z) * Bvec(1)
        F(2) = CC(2) + a * x / (a * x + b * y + c * z) * Avec(2) + b * y / (a * x + b * y + c * z) * Bvec(2)
    end if

    return
end subroutine get_isogonic


subroutine linspace_real(start, end, num, result)
    implicit none
    real(dp), intent(in)  :: start, end
    integer, intent(in)   :: num
    real(dp), intent(out) :: result(num)
    real(dp)              :: step
    integer               :: i
    
    step = (end - start) / (num - 1)
    
    do i = 1, num
      result(i) = start + (i - 1) * step
    end do
    return
end subroutine linspace_real


subroutine point_in_rectangle(is_inside, x, y, ex_e, ey_e)
    ! --- Check if point inside rectangle. OBS rectangle can not be sheared ---
    implicit none

    ! Intent out
    logical, intent(out)  :: is_inside

    ! Intent in
    real(dp), intent(in)  :: x, y, ex_e(4), ey_e(4)

    ! Subroutine variables
    real(dp)              :: minx, maxx, miny, maxy, pad=1d-12

    minx = minval(ex_e)
    maxx = maxval(ex_e)
    miny = minval(ey_e)
    maxy = maxval(ey_e)

    ! Boolean. True if (x,y) inside polygon
    is_inside = (x.ge.minx-pad .and. x.le.maxx+pad .and. y.ge.miny-pad .and. y.le.maxy+pad)


    return
end subroutine point_in_rectangle


subroutine point_in_arbitrary_rectangle(is_inside, x, y, ex_e, ey_e)
    ! --- Check if point (x,y) is inside rectangle or at element edge. 
    !     OBS rectangle can have an arbitrary form (but with straight edges) ---
    implicit none

    ! Intent out
    logical, intent(out)  :: is_inside

    ! Intent in
    real(dp), intent(in)  :: x, y, ex_e(4), ey_e(4)

    ! Subroutine variables
    real(dp)              :: ex_e_loc(4), ey_e_loc(4), x1,x2,x3,y1,y2,y3
    
    ! Divide element in two triangles and use point in triangle function

    ! 4 --- 3      4 --- 3         3
    ! |     |      |    /       /  |
    ! |     |  ->  |(1)/       /(2)|
    ! |     |      |  /       /    |
    ! 1 --- 2      1         1 --- 2

    is_inside = .false.

    ! Add margins to ex_e and ey_e to find point also at element edge
    ex_e_loc = ex_e + [-1d0,  1d0, 1d0, -1d0]*1d-12
    ey_e_loc = ey_e + [-1d0, -1d0, 1d0,  1d0]*1d-12

    ! First triangle
    x1=ex_e_loc(1); x2=ex_e_loc(3); x3=ex_e_loc(4)
    y1=ey_e_loc(1); y2=ey_e_loc(3); y3=ey_e_loc(4)
    call point_within_triangle(is_inside,x,y,x1,x2,x3,y1,y2,y3)
    if (is_inside) then
        return
    endif  
    
    ! Second triangle
    x1=ex_e_loc(1); x2=ex_e_loc(2); x3=ex_e_loc(3)
    y1=ey_e_loc(1); y2=ey_e_loc(2); y3=ey_e_loc(3)
    call point_within_triangle(is_inside,x,y,x1,x2,x3,y1,y2,y3)
    if (is_inside) then
        return
    endif 

    return
end subroutine point_in_arbitrary_rectangle


subroutine point_within_triangle(condition_satisfied,x,y,x1,x2,x3,y1,y2,y3)

    ! Subroutine for determining if point (x,y) is within triangle with vertices xte and yte
    implicit none
  
    ! Intent inout
    logical, intent(inout) :: condition_satisfied  
  
    ! Intent in
    real(dp), intent(in)   :: x, y, x1, x2, x3, y1, y2, y3
  
    ! Subroutine variables
    real(dp)               :: detT, c1, c2, c3
  
    ! Barycentric coordinates of (x,y) in the triangle
    detT = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3)
    c1   = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3))/detT
    c2   = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3))/detT
    c3   = 1d0 - c1 - c2
  
    if ((0d0.le.c1 .and. c1.le.1d0) .and. (0d0.le.c2 .and. c2.le.1d0) .and. (0d0.le.c3 .and. c3.le.1d0) &
    .and. ((c1+c2+c3).le.1d0)) then
      condition_satisfied = .true.    
    endif
  
    return
  end subroutine point_within_triangle


end module ls_reconstruction_routines
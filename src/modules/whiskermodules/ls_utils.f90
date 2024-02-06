module ls_utils
    ! Last modified 
    ! E. Jacobsson 2023-09-04
      
    ! mf_datatypes
    use mf_datatypes

    ! Somelib
    use matrix_util

    ! Whiskerlib
    use ls_sorting_routines
    use ls_reconstruction_routines
    use ls_vp_routines
    use ls_types    
    use mesh_module
    use diffusion_types
    use diffusion

    implicit none

contains


subroutine init_level_set_function(a,x,y,r,coord)
    ! --- Init level set function as a circle with mid point (x,y) and radius r ---
    implicit none

    ! Intent inout 
    real(dp), intent(inout) :: a(:)

    ! Intent in
    real(dp), intent(in)    :: x, y, r, coord(:,:)  

    ! Level set function
    a = sqrt((coord(1,:) - x)**2d0 + (coord(2,:) - y)**2d0) - r 

    return 
end subroutine init_level_set_function


subroutine init_level_set_function_ellips(a,x,y,h1,h2,coord)
    ! --- Init level set function as a circle with mid point (x,y) and radius r ---
    implicit none

    ! Intent inout 
    real(dp), intent(inout) :: a(:)

    ! Intent in
    real(dp), intent(in)    :: x, y, h1, h2, coord(:,:)

    ! Level set function
    a = ((((coord(1,:)-x)/h1)**2 + ((coord(2,:)-y)/h2)**2)) - 1d0

    return 
end subroutine init_level_set_function_ellips


subroutine interaction_correction(a,ngrains)
    ! --- Interaction-correction step to make level set functions compatible ---
    implicit none
    
    ! Intent inout
    real(dp),intent(inout) :: a(:,:)
    
    ! Intent in
    integer                :: ngrains
    
    ! Subroutine variables
    integer                :: grain_idx(ngrains), k, ierr, g
    real(dp)               :: sign=-1d0, other_grains(ngrains-1)
    real(dp), allocatable  :: a_copy(:,:)
    
    ! Allocate a_copy    
    allocate(a_copy(size(a,1),size(a,2)),stat=ierr)
    a_copy = a
    
    ! Define grain idx
    do k=1,ngrains
        grain_idx(k) = k
    enddo
    
    ! Interaction-correction node wise
    do g = 1, ngrains
    
        ! Find what other grains are present at the grain g
        other_grains = pack(grain_idx, grain_idx .ne. g)
    
        ! Correct a(:,g)
        a(:,g) = 0.5d0 * (a_copy(:,g) - minval(a_copy(:,other_grains), dim=2))
    
        ! Ensure that no ls intersect exactly at a node    
        if (any(abs(a(:,g)) .lt. 1d-13)) then
            print *, 'Adjusting LS since it is intersecting exactly at a node'
            where (abs(a(:,g)) .lt. 1d-13) a(:,g) = sign * 1d-13
            sign = 1d0
        end if
    
    enddo
    
    ! Deallocate subroutine variables
    deallocate(a_copy)
    
    return
end subroutine interaction_correction


subroutine get_interface(line_ex, line_ey, line_elms, line_seg, ed, mesh)
    ! --- Extracting zero isocontour interface for a single level set function ---
    implicit none

    ! Intent inout
    real(dp),intent(out)           :: line_ex(:,:), line_ey(:,:)
    integer,intent(out)            :: line_elms(:), line_seg

    ! Intent in
    type(mesh_system), intent(in)  :: mesh
    real(dp), intent(in)           :: ed(:,:)

    ! Subroutine variables 
    integer                        :: nint_elms, iseg, ie, k, counter, inod, ierr, n1(2), n2(2)
    real(dp)                       :: isov, ed_e(4), ex_e(4), ey_e(4), intersection_points(4,2), edges_intersected(4)
    real(dp)                       :: x1,y1,x2,y2,t1,t2,alpha,intersection_x,intersection_y
    integer, allocatable           :: interface_elms(:)
    


    ! 1) Find elements containing interface (the level set zero isocontour)
    nint_elms = count(maxval(ed,1)*minval(ed,1).lt.0d0) ! Number of interface elements
    allocate(interface_elms(nint_elms),stat=ierr)
    interface_elms = pack(mesh%elmidx,maxval(ed,1)*minval(ed,1).lt.0d0)

    ! 2) Find the line segments in each interface element
    iseg = 0
    isov = 0d0  ! isocontour to be found
    
    do k=1,nint_elms

        ! Interface element
        ie = interface_elms(k)

        ! Short form of local coordinates and nodal values
        ed_e = ed(:,ie); ex_e = mesh%ex(:,ie); ey_e = mesh%ey(:,ie)

        ! Loop through all edges of the element and find intersection points
        intersection_points = 0d0
        edges_intersected   = 0d0
        counter             = 0
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
            iseg            = iseg + 1
            line_ex(iseg,:) = intersection_points(edges_intersected(1:counter),1)
            line_ey(iseg,:) = intersection_points(edges_intersected(1:counter),2)
            line_elms(iseg) = ie

        elseif (counter.eq.4) then

            ! Two line segments in element            
            if (norm2(sign([1d0,1d0,1d0,1d0],ed_e) - [-1d0, 1d0, -1d0, 1d0]) .lt. 1d-13) then
                n1 = [1,2]
                n2 = [3,4]
            elseif (norm2(sign([1d0,1d0,1d0,1d0],ed_e) - [1d0, -1d0, 1d0, -1d0]) .lt. 1d-13) then
                n1 = [1,4]
                n2 = [2,3]
            endif   

            iseg            = iseg + 1
            line_ex(iseg,:) = intersection_points(n1,1)
            line_ey(iseg,:) = intersection_points(n1,2)
            line_elms(iseg) = ie

            iseg            = iseg + 1
            line_ex(iseg,:) = intersection_points(n2,1)
            line_ey(iseg,:) = intersection_points(n2,2)
            line_elms(iseg) = ie
        else
            print *, 'Error: Level set cutting neither 2 or 4 edges'
        endif
    enddo

    ! Number of line segments
    line_seg = iseg

    ! Deallocate subroutine variables
    deallocate(interface_elms)

    return
end subroutine get_interface


subroutine get_interface_spatial(line_ex_s, line_ey_s, line_elms, line_seg, ed, mesh)
    ! --- Extracting zero isocontour interface for a single level set function ---
    implicit none

    ! Intent inout
    real(dp), intent(out)          :: line_ex_s(:,:), line_ey_s(:,:)
    integer, intent(out)           :: line_seg, line_elms(:)

    ! Intent in
    type(mesh_system), intent(in)  :: mesh
    real(dp), intent(in)           :: ed(:,:)

    ! Subroutine variables 
    integer                        :: nint_elms, iseg, ie, k, counter, inod, ierr, n1(2), n2(2), edges_intersected(4)
    real(dp)                       :: isov, ed_e(4), ex_e(4), ey_e(4), intersection_points(4,2)
    real(dp)                       :: x1,y1,x2,y2,t1,t2,alpha,intersection_x,intersection_y
    integer, allocatable           :: interface_elms(:)

    ! 1) Find elements containing the interface (the level set zero isocontour)
    nint_elms      = count(maxval(ed,1)*minval(ed,1).lt.0d0) ! Number of interface elements
    allocate(interface_elms(nint_elms),stat=ierr)
    interface_elms = pack(mesh%elmidx,maxval(ed,1)*minval(ed,1).lt.0d0)

    ! 2) Find the line segments in each interface element
    iseg = 0
    isov = 0d0 ! isocontour to be found
    
    do k=1,nint_elms

        ! Interface element
        ie = interface_elms(k)

        ! Short form of local coordinates and nodal values
        ed_e = ed(:,ie); ex_e = mesh%newex(:,ie); ey_e = mesh%newey(:,ie)

        ! Loop through all edges of the element and find intersection points
        intersection_points = 0d0
        edges_intersected   = 0
        counter             = 0
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
            iseg              = iseg + 1
            line_ex_s(iseg,:) = intersection_points(edges_intersected(1:counter),1)
            line_ey_s(iseg,:) = intersection_points(edges_intersected(1:counter),2)
            line_elms(iseg)   = ie

        elseif (counter.eq.4) then

            ! Two line segments in element            
            if (norm2(sign([1d0,1d0,1d0,1d0],ed_e) - [-1d0, 1d0, -1d0, 1d0]) .lt. 1d-13) then
                n1 = [1,2]
                n2 = [3,4]
            elseif (norm2(sign([1d0,1d0,1d0,1d0],ed_e) - [1d0, -1d0, 1d0, -1d0]) .lt. 1d-13) then
                n1 = [1,4]
                n2 = [2,3]
            endif   

            iseg              = iseg + 1
            line_ex_s(iseg,:) = intersection_points(n1,1)
            line_ey_s(iseg,:) = intersection_points(n1,2)
            line_elms(iseg)   = ie

            iseg            = iseg + 1
            line_ex_s(iseg,:) = intersection_points(n2,1)
            line_ey_s(iseg,:) = intersection_points(n2,2)
            line_elms(iseg) = ie
        else
            print *, 'Error in ls interpolate interface: Level set cutting neither 2 nor 4 edges'
        endif
    enddo

    ! Number of line segments
    line_seg = iseg

    ! Deallocate subroutine variables
    deallocate(interface_elms)

    return
end subroutine get_interface_spatial


subroutine interface_reconstruction(lssys,mesh)
    ! --- Routine for making interface reconstruction, a post-processing step to ensure that interfaces compatibel. 
    !     Note that void regions exist even after the interaction-correction step ---
    implicit none

    ! Intent inout
    type(ls_system), intent(inout) :: lssys
    
    ! Intent in
    type(mesh_system), intent(in)  :: mesh

    ! Subroutine variables    
    real(dp)                       :: a_copy(size(lssys%a,1)), critical_length
    integer                        :: g, ie, nods(mesh%nodel)    

    ! Critical length
    critical_length = 0.4d0*mesh%min_elm_width

    ! Fix grain boundaries
    call gb_reconstruction(lssys,mesh)


    ! Fix triple junctions    
    call tp_reconstruction(lssys,mesh,critical_length)

    
    ! Make sure that all nods in the non interface elements have the same sign, i.e. make sign of a_gr more accurate        
    do g=1,lssys%ngrains
        ! Find elements that dont contain a line segment
        a_copy = lssys%a(:,g)      
        do ie=1,mesh%nelm
            if (lssys%int_elms(ie,g).eq.0) then
                ! Non-interface element -> all nods should have same sign            
                nods = mesh%enod(:,ie)
                lssys%a(nods,g) = sign(a_copy(nods), sum(lssys%ed(:,ie,g)))
            endif
        enddo
    enddo
    

    return
end subroutine interface_reconstruction


subroutine interface_reconstruction_spatial(lssys,mesh)
    ! --- Routine for making an interface reconstruction, a post-processing step to ensure that interfaces are compatibel. 
    !     Note that void regions exist even after the interaction-correction step. In this routine, all void regions are
    !     eliminated ---
    implicit none

    ! Intent inout
    type(ls_system), intent(inout) :: lssys
    
    ! Intent in
    type(mesh_system), intent(in)  :: mesh

    ! Subroutine variables    
    real(dp)                       :: a_copy(size(lssys%a,1)), critical_length
    integer                        :: g, ie, nods(mesh%nodel)    

    ! Critical length
    critical_length = 0.4d0*mesh%min_elm_width

    ! Fix grain boundaries
    call gb_reconstruction_spatial(lssys,mesh)

    ! Fix triple junctions
    call tp_reconstruction_spatial(lssys,mesh)

    
    ! Make sure that all nods in the non interface elements have the same sign, i.e. make sign of a_gr more accurate        
    do g=1,lssys%ngrains
        ! Find elements that dont contain a line segment
        a_copy = lssys%a(:,g)      
        do ie=1,mesh%nelm
            if (lssys%int_elms(ie,g).eq.0) then
                ! Non-interface element -> all nods should have same sign            
                nods = mesh%enod(:,ie)
                lssys%a(nods,g) = sign(a_copy(nods), sum(lssys%ed(:,ie,g)))
            endif
        enddo
    enddo
    

    return
end subroutine interface_reconstruction_spatial


subroutine sort_lines(line_ex,line_ey,tplines,line_seg,sep_lines,nsep_lines, mesh)
    ! --- Routine for sorting line segments, removing too short lines and adding lines to too long line segments ---
    implicit none

    ! Intent inout
    real(dp),intent(inout)              :: line_ex(:,:), line_ey(:,:)
    integer, intent(inout)              :: line_seg
    logical, intent(inout)              :: tplines(:)
    integer,intent(inout)               :: sep_lines(:,:), nsep_lines

    ! Intent in
    type(mesh_system), intent(in)       :: mesh

    ! Subroutine variables
    integer, allocatable                :: separate_lines_idx(:), sep_l_idx_ext(:,:)
    real(dp)                            :: critical_length
    
    ! 1) Make sure line segments are connected
    call connect_line_segments(line_ex,line_ey,tplines,separate_lines_idx,line_seg,mesh, sep_l_idx_ext)

    ! 2) Remove line segments shorter than threshold
    sep_l_idx_ext = 0
    critical_length = 0.6d0*mesh%min_elm_width
    call remove_line_segments(line_ex,line_ey,tplines,critical_length,separate_lines_idx,line_seg,sep_l_idx_ext)

    ! Sep lines
    nsep_lines                = size(sep_l_idx_ext,2)
    sep_lines(1:nsep_lines,:) = transpose(sep_l_idx_ext)
    ! OBS this does not work with sep_lines!!

    ! 3) Add line segements longer than threshold
    critical_length = 1.2d0*mesh%min_elm_width
    call add_line_segments(line_ex, line_ey, critical_length, line_seg, sep_lines, nsep_lines)


    ! Deallocate
    deallocate(separate_lines_idx)
    deallocate(sep_l_idx_ext)


    return
end subroutine sort_lines


subroutine reinit_level_set(a_g,line_ex,line_ey,closest_line,nseg,mesh)
    ! --- Routine for reinitializing level set function. Valid for a quad mesh ---
    implicit none

    ! Intent inout
    real(dp),intent(inout)         :: a_g(:), line_ex(:,:), line_ey(:,:)
    integer, intent(inout)         :: closest_line(:), nseg

    ! Intent in
    type(mesh_system), intent(in)  :: mesh    

    ! Subroutine variables
    integer                        :: inod, iIntElm, ierr
    real(dp)                       :: P(2), A(2), B(2), E(2),v(2), u(2), t
    real(dp), allocatable          :: seg_dists(:)


    ! Loop through all nodes and for each node compute distance to all line segments
    do inod = 1,mesh%nnod

        ! Point from which distances to line segments are computed
        P = mesh%coord(:,inod)

        ! Compute distance to all line segments for given node point P
        allocate(seg_dists(nseg),stat=ierr)
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

        ! Store closest line
        closest_line(inod) = minloc(seg_dists, 1)
    enddo

    ! Deallocate    
    deallocate(seg_dists)

    return
end subroutine reinit_level_set


subroutine reinit_level_set_spatial(a_g,line_ex,line_ey,closest_line,nseg,mesh)
    ! --- Routine for reinitializing level set function in deformed configuration ---
    implicit none

    ! Intent inout
    real(dp),intent(inout)         :: a_g(:), line_ex(:,:), line_ey(:,:)
    integer, intent(inout)         :: closest_line(:), nseg

    ! Intent in
    type(mesh_system), intent(in)  :: mesh    

    ! Subroutine variables
    integer                        :: inod, iIntElm, ierr
    real(dp)                       :: P(2), A(2), B(2), E(2),v(2), u(2), t
    real(dp), allocatable          :: seg_dists(:)


    ! Loop through all nodes. For each node -> compute distance to all line segments
    do inod = 1,mesh%nnod

        ! Coordinates of node
        P = mesh%newcoord(:,inod)

        ! Compute distance to all line segments for given node P
        allocate(seg_dists(nseg),stat=ierr)
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

        ! Store closest line
        closest_line(inod) = minloc(seg_dists, 1)
    enddo

    ! Deallocate    
    deallocate(seg_dists)

    return
end subroutine reinit_level_set_spatial


subroutine compute_common_vp(lssys, mesh, diffsys)
     ! --- Routine for computing a common velocity field ---
    implicit none

    ! Intent inout
    type(diffusion_system), intent(inout) :: diffsys
    type(ls_system), intent(inout)        :: lssys
    
    ! Intent in
    type(mesh_system), intent(in)  :: mesh

    ! Subrouitine variables
    integer  :: g, ie, Igrain_gp(mesh%nrgp), inod
    real(dp) :: ed_e_gr(mesh%nodel,lssys%ngrains), jed_e_gr(mesh%nodel,lssys%ngrains), bcvaled_e_gr(mesh%nodel,lssys%ngrains)
    real(dp) :: agp_e_gr(mesh%nrgp,lssys%ngrains), jgp_e_gr(mesh%nrgp,lssys%ngrains), bcvalgp_e_gr(mesh%nrgp,lssys%ngrains)
    integer  :: g_cols(2), minIdx(1), i, snsn_closest(1)
    real(dp) :: xnod, ynod, xint, yint, grain_nod
    ! integer, allocatable :: indices(:)
    ! real(dp) :: tp_points_possible(400,2), tp_points(400,2)
    real(dp) :: mean_x, model_edge_x_distance, mgb1pos, mgb2pos, mapply, mzetagbapply
    real(dp) :: right_distance, left_distance, mgbpos(2), numsnsngb, snsn_gb_distance

    ! 1) Extract magnitude of normal flux at all interface points for all level sets. (Flux outward from the grain).
    ! Should be obtained from diffusion problem
    ! lssys%jint = 0d0
    ! do g=1,lssys%ngrains    
    !     lssys%jint(1:lssys%line_seg(g),g) = lssys%EEvec(g)
    ! enddo



    ! New closest interface point
    do g=1,lssys%ngrains
        do inod=1,mesh%nnod

            ! g_cols
            g_cols = [2*(g-1) + 1, 2*(g-1) + 2]

            ! Coordinates of nod
            xnod = mesh%coord(1,inod)
            ynod = mesh%coord(2,inod)            

            ! Find coordinates of closest point on interface
            minIdx =  minloc(sqrt((xnod - lssys%line_ex(1:lssys%line_seg(g),g_cols(1)))**2 + & 
            (ynod - lssys%line_ey(1:lssys%line_seg(g),g_cols(1)))**2))
            xint = lssys%line_ex(minIdx(1),g_cols(1))
            yint = lssys%line_ey(minIdx(1),g_cols(1))

            ! Find corresponding nod in grain_mesh
            minIdx = minloc(sqrt((diffsys%grain_meshes(g)%coord(1,:) - xint)**2 + (diffsys%grain_meshes(g)%coord(2,:) - yint)**2))
            grain_nod = minIdx(1)

            ! if (abs(mesh%coord(1,inod)-0.00036d0).lt.1d-5 .and. abs(mesh%coord(2,inod)-0.00044d0).lt.1d-5) then
            !     print *, 'nod coord', mesh%coord(:,inod)
            !     print *, 'clnod coord', lssys%grain_meshes(g)%coord(:,grain_nod)                    
            !     print *, 'as'
            ! endif            

            ! See if closest interface node is not present in bcnod
            if (all(abs(diffsys%grain_meshes(g)%bcnod(:,1) - grain_nod).gt.0)) then
                ! Closest interface node does not exist in bcnod

                ! if (abs(mesh%coord(1,inod)-0.00104d0).lt.1d-4 .and. abs(mesh%coord(2,inod)-0.00040d0).lt.1d-4) then
                !     print *, 'nod coord', mesh%coord(:,inod)
                !     print *, 'clnod coord', lssys%grain_meshes(g)%coord(:,grain_nod)                    
                !     print *, 'as'
                ! endif

                lssys%jnod(inod,g)  = 0d0
                ! lssys%bcval(inod,g) = 0.0111d0

                ! Closest interface node exist in bcnod     
                ! Find local position of nod in bcnod
                minIdx = minloc(abs(diffsys%grain_meshes(g)%bcnod(:,1) - grain_nod))                
                lssys%bcval(inod,g) = diffsys%grain_meshes(g)%bcval(minIdx(1))

            else
                ! Closest interface node exist in bcnod     
                ! Find local position of nod in bcnod
                minIdx = minloc(abs(diffsys%grain_meshes(g)%bcnod(:,1) - grain_nod))
                lssys%jnod(inod,g)  = diffsys%grain_meshes(g)%jint(minIdx(1))
                lssys%bcval(inod,g) = diffsys%grain_meshes(g)%bcval(minIdx(1))
            endif
        enddo
    enddo    


    ! ! 2) Transfer interface point normal velocities to all nods by closest line segment
    ! do g=1,lssys%ngrains
    !     ! lssys%jnod(:,g) = lssys%jint(lssys%closest_line(:,g),g)
    !     lssys%jnod(:,g)  = lssys%grain_meshes(g)%jint(lssys%closest_line(:,g)+1)  ! Must be changed
    !     lssys%bcval(:,g) = lssys%grain_meshes(g)%bcval(lssys%closest_line(:,g)+1)
    ! enddo
    

    ! 3) Loop through all elements and compute velocity field in each gauss point
    do ie=1,mesh%nelm

        ! ed_e_gr  : [elmvals_g1 , elmvals_g2 , ....., elmvals_gn]
        ! jed_e_gr : [jelmvals_g1, jelmvals_g2, ....., jelmvals_gn]


        do g=1,lssys%ngrains
            ed_e_gr(:,g) = lssys%ed(:,ie,g)
        enddo

        do g=1,lssys%ngrains
            jed_e_gr(:,g) = lssys%jnod(mesh%enod(:,ie),g)
        enddo

        do g=1,lssys%ngrains
            bcvaled_e_gr(:,g) = lssys%bcval(mesh%enod(:,ie),g)
        enddo

        ! 
        where (abs(jed_e_gr(:,g)).lt.1d-22) jed_e_gr(:,g) = 0d0

        ! Convert from nodvals to gauss point vals
        call elm2D4_nodmat_to_gpmat(agp_e_gr, ed_e_gr, lssys%ngrains)
        call elm2D4_nodmat_to_gpmat(jgp_e_gr, jed_e_gr, lssys%ngrains)
        call elm2D4_nodmat_to_gpmat(bcvalgp_e_gr, bcvaled_e_gr, lssys%ngrains)


        ! Find which grain each gauss point in the element belong to
        Igrain_gp = minloc(agp_e_gr,2)

        ! Compute velocity in all gauss points of element
        call lvlset2D4_global_vp(lssys%vp(:,:,ie),mesh%coord(:,mesh%enod(:,ie)),ed_e_gr,jgp_e_gr,agp_e_gr,&
        bcvalgp_e_gr,Igrain_gp, lssys, diffsys, ie)    
        
        ! Apply mobility based on xdistance to Sn gb's
        mean_x  = sum(mesh%coord(1,mesh%enod(:,ie)))/size(mesh%enod(:,ie))
                

        ! Sn/Sn gb positions
        numsnsngb = 4
        mgbpos(1) = 0d0

        ! g = 6
        ! g_cols = [2*(g-1) + 1, 2*(g-1) + 2]
        ! mgbpos(2) = maxval(lssys%line_ex(:,g_cols))

        ! g = 7
        ! g_cols = [2*(g-1) + 1, 2*(g-1) + 2]
        ! mgbpos(3) = maxval(lssys%line_ex(:,g_cols))

        ! mgbpos(4) = mesh%model_width      
        mgbpos(2) = mesh%model_width
        

        
        ! Closest distance to Sn/Sn gb
        snsn_gb_distance = minval(abs(mgbpos - mean_x))
        snsn_closest     = minloc(abs(mgbpos - mean_x))        

        ! Sn/Sn gb mobility
        ! mzetagbapply = snsn_closest(1)/2d0 * lssys%mzetagb
        mzetagbapply = lssys%mzetagb

        mapply = lssys%mzeta + (mzetagbapply - lssys%mzeta)*exp(-lssys%malpha*snsn_gb_distance)




        ! ! Different mobilities
        ! right_distance = abs(mesh%model_width - mean_x)
        ! left_distance  = abs(0d0 - mean_x)
        ! if (right_distance.lt.left_distance) then
        !     mapply = lssys%mzeta + (lssys%mzetagb*2d0 - lssys%mzeta)*exp(-lssys%malpha*right_distance)            
        ! else
        !     mapply = lssys%mzeta + (lssys%mzetagb/2d0 - lssys%mzeta)*exp(-lssys%malpha*left_distance)
        ! endif

        ! If distance to Sn less than threshold -> apply Sn/IMC mobility (4 IMC grains)
        ! if (abs(sum(lssys%ed(:,ie,6))/4d0).lt.2.5d-5*2d0 .or. abs(sum(lssys%ed(:,ie,7))/4d0).lt.2.5d-5*2d0 &
        ! .or. abs(sum(lssys%ed(:,ie,8))/4d0).lt.2.5d-5*2d0) then
        !     lssys%vp(:,:,ie) = lssys%vp(:,:,ie)*mapply
        ! else
        !     lssys%vp(:,:,ie) = lssys%vp(:,:,ie)*lssys%mzeta
        ! endif

        ! If distance to Sn less than threshold -> apply Sn/IMC mobility (3 IMC grains)
        if (abs(sum(lssys%ed(:,ie,5))/4d0).lt.2.5d-5*2d0) then
            lssys%vp(:,:,ie) = lssys%vp(:,:,ie)*mapply
        else
            lssys%vp(:,:,ie) = lssys%vp(:,:,ie)*lssys%mzeta
        endif

    enddo

    return
end subroutine compute_common_vp


subroutine compute_common_vp_spatial(lssys, mesh, diffsys)
    ! --- Routine for computing a common velocity field in deformed configuration ---
   implicit none

   ! Intent inout
   type(diffusion_system), intent(inout) :: diffsys
   type(ls_system), intent(inout)        :: lssys
   
   ! Intent in
   type(mesh_system), intent(in)  :: mesh

   ! Subrouitine variables
   integer  :: g, ie, Igrain_gp(mesh%nrgp), inod
   real(dp) :: ed_e_gr(mesh%nodel,lssys%ngrains), jed_e_gr(mesh%nodel,lssys%ngrains), bcvaled_e_gr(mesh%nodel,lssys%ngrains)
   real(dp) :: agp_e_gr(mesh%nrgp,lssys%ngrains), jgp_e_gr(mesh%nrgp,lssys%ngrains), bcvalgp_e_gr(mesh%nrgp,lssys%ngrains)
   integer  :: g_cols(2), minIdx(1), snsn_closest(1)
   real(dp) :: xnod, ynod, xint, yint, grain_nod
   real(dp) :: mean_x, mapply, mzetagbapply
   real(dp) :: mgbpos(2), numsnsngb, snsn_gb_distance

   ! New closest interface point
   do g=1,lssys%ngrains
       do inod=1,mesh%nnod

           ! g_cols
           g_cols = [2*(g-1) + 1, 2*(g-1) + 2]

           ! Coordinates of nod
           xnod = mesh%newcoord(1,inod)
           ynod = mesh%newcoord(2,inod)            

           ! Find coordinates of closest point on interface
           minIdx =  minloc(sqrt((xnod - lssys%line_ex(1:lssys%line_seg(g),g_cols(1)))**2 + & 
           (ynod - lssys%line_ey(1:lssys%line_seg(g),g_cols(1)))**2))
           xint = lssys%line_ex(minIdx(1),g_cols(1))
           yint = lssys%line_ey(minIdx(1),g_cols(1))

           ! Find corresponding nod in grain_mesh
           minIdx = minloc(sqrt((diffsys%grain_meshes(g)%coord(1,:) - xint)**2 + (diffsys%grain_meshes(g)%coord(2,:) - yint)**2))
           grain_nod = minIdx(1)

           ! See if closest interface node is not present in bcnod
           if (all(abs(diffsys%grain_meshes(g)%bcnod(:,1) - grain_nod).gt.0)) then
               ! Closest interface node does not exist in bcnod

               lssys%jnod(inod,g) = 0d0
               ! lssys%bcval(inod,g) = 0.0111d0

               ! Closest interface node exist in bcnod     
               ! Find local position of nod in bcnod
               minIdx = minloc(abs(diffsys%grain_meshes(g)%bcnod(:,1) - grain_nod))                
               lssys%bcval(inod,g) = diffsys%grain_meshes(g)%bcval(minIdx(1))

           else
               ! Closest interface node exist in bcnod     
               ! Find local position of nod in bcnod
               minIdx = minloc(abs(diffsys%grain_meshes(g)%bcnod(:,1) - grain_nod))
               lssys%jnod(inod,g)  = diffsys%grain_meshes(g)%jint(minIdx(1))
               lssys%bcval(inod,g) = diffsys%grain_meshes(g)%bcval(minIdx(1))
           endif
       enddo
   enddo    

   ! Loop through all elements and compute velocity field in each gauss point
   do ie=1,mesh%nelm

       ! ed_e_gr  : [elmvals_g1 , elmvals_g2 , ....., elmvals_gn]
       ! jed_e_gr : [jelmvals_g1, jelmvals_g2, ....., jelmvals_gn]


       do g=1,lssys%ngrains
           ed_e_gr(:,g) = lssys%ed(:,ie,g)
       enddo

       do g=1,lssys%ngrains
           jed_e_gr(:,g) = lssys%jnod(mesh%enod(:,ie),g)
       enddo

       do g=1,lssys%ngrains
           bcvaled_e_gr(:,g) = lssys%bcval(mesh%enod(:,ie),g)
       enddo

       ! 
       where (abs(jed_e_gr(:,g)).lt.1d-22) jed_e_gr(:,g) = 0d0

       ! Convert from nodvals to gauss point vals
       call elm2D4_nodmat_to_gpmat(agp_e_gr, ed_e_gr, lssys%ngrains)
       call elm2D4_nodmat_to_gpmat(jgp_e_gr, jed_e_gr, lssys%ngrains)
       call elm2D4_nodmat_to_gpmat(bcvalgp_e_gr, bcvaled_e_gr, lssys%ngrains)


       ! Find which grain each gauss point in the element belong to
       Igrain_gp = minloc(agp_e_gr,2)

       ! Compute velocity in all gauss points of element
       call lvlset2D4_global_vp(lssys%vp(:,:,ie),mesh%newcoord(:,mesh%enod(:,ie)),ed_e_gr,jgp_e_gr,agp_e_gr,&
       bcvalgp_e_gr,Igrain_gp, lssys, diffsys, ie)
       
       ! Apply mobility based on xdistance to Sn gb's
       mean_x  = sum(mesh%newcoord(1,mesh%enod(:,ie)))/size(mesh%enod(:,ie))
               
       ! Sn/Sn gb positions
       numsnsngb = 4
       mgbpos(1) = 0d0

       ! g = 6
       ! g_cols = [2*(g-1) + 1, 2*(g-1) + 2]
       ! mgbpos(2) = maxval(lssys%line_ex(:,g_cols))

       ! g = 7
       ! g_cols = [2*(g-1) + 1, 2*(g-1) + 2]
       ! mgbpos(3) = maxval(lssys%line_ex(:,g_cols))

       ! mgbpos(4) = mesh%model_width      
       mgbpos(2) = mesh%model_width
       
       ! Closest distance to Sn/Sn gb
       snsn_gb_distance = minval(abs(mgbpos - mean_x))
       snsn_closest     = minloc(abs(mgbpos - mean_x))        

       ! Sn/Sn gb mobility
       mapply = lssys%mzeta + (lssys%mzetagb - lssys%mzeta)*exp(-lssys%malpha*snsn_gb_distance)

       ! If distance to Sn less than threshold -> apply Sn/IMC mobility
       if (abs(sum(lssys%ed(:,ie,5))/4d0).lt.2.5d-5*2d0) then
           lssys%vp(:,:,ie) = lssys%vp(:,:,ie)*mapply
       else
           lssys%vp(:,:,ie) = lssys%vp(:,:,ie)*lssys%mzeta
       endif

   enddo

   return
end subroutine compute_common_vp_spatial


subroutine remove_ls(lssys)
    ! --- Routine for removing level set if it is larger than zero everywhere ---
    implicit none

    ! Intent inout
    type(ls_system), intent(inout) :: lssys    


    ! Subroutine variables
    integer              :: i, g, ngrains, nkeep, counter, ierr
    integer, allocatable :: gidx(:), keep_g(:), keep_lines(:)

    ! ngrains
    ngrains = lssys%ngrains

    ! Grain idx
    gidx = pack([(i, i = 1, ngrains)], [(.true., i = 1, ngrains)])



    ! Determine what grains to be removed
    allocate(keep_g(ngrains),stat=ierr)
    counter = 1
    do g=1,ngrains
        if (.not. all(lssys%a(:,g).gt.0)) then
            keep_g(counter) = g
            counter = counter + 1
        endif
    enddo 

    ! Reduce keep_g
    nkeep  = count(keep_g.gt.0)
    keep_g = keep_g(1:nkeep)

    if (ngrains-nkeep.gt.0) then
        print *, 'Removing ', ngrains-nkeep, 'grains'
    

        ! Keep lines
        allocate(keep_lines(nkeep*2),stat=ierr)
        keep_lines(1:2*nkeep-1:2) = (keep_g-1)*2 + 1
        keep_lines(2:2*nkeep:2)   = (keep_g-1)*2 + 2


        ! Reshape matrices
        lssys%a            = lssys%a(:,keep_g)
        lssys%ed           = lssys%ed(:,:,keep_g)
        lssys%jnod         = lssys%jnod(:,keep_g)        
        lssys%jed          = lssys%jed(:,:,keep_g)
        lssys%line_ex      = lssys%line_ex(:,keep_lines)
        lssys%line_ey      = lssys%line_ey(:,keep_lines)
        lssys%line_seg     = lssys%line_seg(keep_g)
        lssys%closest_line = lssys%closest_line(:,keep_g)
        lssys%int_elms     = lssys%int_elms(:,keep_g)    
        lssys%line_elms    = lssys%line_elms(:,keep_g)
        lssys%tplines      = lssys%tplines(:,keep_g)
        lssys%sep_lines    = lssys%sep_lines(:,keep_lines)
        lssys%nsep_lines   = lssys%nsep_lines(keep_g)
        lssys%rm_lines     = lssys%rm_lines(:,keep_g)
        lssys%rm_lines_tmp = lssys%rm_lines_tmp(:,keep_g)
        lssys%ngrains      = nkeep
        lssys%material     = lssys%material(keep_g)

    endif


    return
end subroutine remove_ls


subroutine init_random_seed
    implicit none

    integer :: values(1:8), k
    integer, dimension(:), allocatable :: seed

    call date_and_time(values=values)

    call random_seed(size=k)
    allocate(seed(1:k))
    seed(:) = values(8)
    call random_seed(put=seed)

    return
end subroutine init_random_seed


subroutine rand_vec(a, randMin, randMax)
    implicit none

    real(dp), intent(inout) :: a(:)
    real(dp), intent(in) :: randMin, randMax

    call random_number(a)

    !    min ---- max

    a = randMin + (randMax-randMin) * a

    return
end subroutine rand_vec


subroutine GetUniqueIndices(inputArray, uniqueIndices_out)
    implicit none
    integer, intent(in)               :: inputArray(:)
    integer, intent(out), allocatable :: uniqueIndices_out(:)
    integer :: i, j, k
    logical, dimension(size(inputArray)) :: isUnique
    integer, allocatable :: uniqueIndices(:)
  
    !
    allocate(uniqueIndices(size(inputArray)))


    ! Initialize the isUnique array
    isUnique = .TRUE.
  
    ! Find unique indices
    k = 0
    do i = 1, size(inputArray)
      if (isUnique(i)) then
        k = k + 1
        uniqueIndices(k) = i
        do j = i + 1, size(inputArray)
          if (inputArray(i) == inputArray(j)) then
            isUnique(j) = .FALSE.
          end if
        end do
      end if
    end do
  
    ! Resize the uniqueIndices array to the actual number of unique indices
    allocate(uniqueIndices_out(k))
    uniqueIndices_out = uniqueIndices(1:k)

end subroutine GetUniqueIndices

end module ls_utils
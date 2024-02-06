module ls_sorting_routines

    ! Last modified 
    ! E. Jacobsson 2023-09-05
      
    ! mf_datatypes
    use mf_datatypes

    ! whiskerlib
    use mesh_module


    ! --- Contains routines used to sort line segments in level set ---

    implicit none

contains


subroutine connect_line_segments(line_ex, line_ey, tplines, separate_lines_idx, nseg, mesh, separate_lines_idx_extended)
    ! Routine for sorting line segments such that they are connected
    implicit none

    ! Intent inout
    real(dp), intent(inout)             :: line_ex(:,:), line_ey(:,:)
    integer, allocatable, intent(inout) :: separate_lines_idx(:)
    logical, intent(inout)              :: tplines(:)
    integer, allocatable, intent(inout) :: separate_lines_idx_extended(:,:)

    ! Intent in
    integer, intent(in)                 :: nseg
    type(mesh_system), intent(in)       :: mesh    

    ! Subroutine variables
    real(dp), allocatable               :: line_ex_tmp(:,:), line_ey_tmp(:,:)
    real(dp), allocatable               :: line_ex_red(:,:), line_ey_red(:,:), line_ex_remaining(:,:), line_ey_remaining(:,:)
    real(dp)                            :: start_point(2), search_point(2), leftLinePoint(2), rightLinePoint(2)
    integer, allocatable                :: visited_lines(:), separate_lines_idx2(:)
    integer                             :: i, j, ierr, ij_leftmost_point(2), nvisited, nvisited2, line_idx, k, jk(1)
    integer                             :: n_seperate_sets, nseparate_lines
    logical                             :: skip_line, line_not_found, separate_lines
    logical, allocatable                :: tplines_sorted(:)

    
    ! Nonzero lines stored as line_eX_tmp
    allocate(line_ex_tmp(2,nseg),stat=ierr)
    allocate(line_ey_tmp(2,nseg),stat=ierr)
    line_ex_tmp = transpose(line_ex(1:nseg,:))
    line_ey_tmp = transpose(line_ey(1:nseg,:))

    ! Tplines
    allocate(tplines_sorted(size(tplines)),stat=ierr)
    tplines_sorted = tplines

    ! Allocate line_ex_red, line_ey_red and visited_lines
    allocate(line_ex_red(2,nseg),stat=ierr)
    allocate(line_ey_red(2,nseg),stat=ierr)
    allocate(line_ex_remaining(2,nseg),stat=ierr)
    allocate(line_ey_remaining(2,nseg),stat=ierr)
    line_ex_remaining = line_ex_tmp
    line_ey_remaining = line_ey_tmp
    allocate(visited_lines(nseg),stat=ierr)
    line_ex_red   = 0d0
    line_ey_red   = 0d0
    visited_lines = 0
    allocate(separate_lines_idx2(nseg),stat=ierr)
    separate_lines_idx2 = 0

    k         = 0
    nvisited  = 0
    nvisited2 = 0
    do while (nvisited.lt.nseg)
        
        ! Find start point
        jk = minloc([minval(abs(line_ex_remaining - mesh%axisbc(1))), minval(abs(line_ex_remaining - mesh%axisbc(2))), &
        minval(abs(line_ey_remaining - mesh%axisbc(3))), minval(abs(line_ey_remaining - mesh%axisbc(4)))])

        if (jk(1).eq.1) then
            ij_leftmost_point = minloc(abs(line_ex_remaining - mesh%axisbc(1)))
        elseif (jk(1).eq.2) then
            ij_leftmost_point = minloc(abs(line_ex_remaining - mesh%axisbc(2)))
        elseif (jk(1).eq.3) then
            ij_leftmost_point = minloc(abs(line_ey_remaining - mesh%axisbc(3)))
        elseif (jk(1).eq.4) then
            ij_leftmost_point = minloc(abs(line_ey_remaining - mesh%axisbc(4)))
        endif

        ! ij_leftmost_point  = minloc(line_ex_remaining)
        start_point        = [line_ex_remaining(ij_leftmost_point(1),ij_leftmost_point(2)), &
                              line_ey_remaining(ij_leftmost_point(1),ij_leftmost_point(2))]

        ! Search for next line
        do i=nvisited+1,nseg

            if (i.eq.nvisited+1) then
                search_point = start_point
            else
                search_point = [line_ex_red(2,i-1), line_ey_red(2,i-1)]
            endif

            ! Loop through all non-visited line segments
            do line_idx=1,nseg

                ! Check if line_idx has already been visited
                skip_line = .false.
                do j = 1, nvisited2
                    if (line_idx .eq. visited_lines(j)) then
                        skip_line = .true.
                        exit
                    endif
                enddo

                if (.not. skip_line) then
                    leftLinePoint  = [line_ex_tmp(1,line_idx), line_ey_tmp(1,line_idx)]
                    rightLinePoint = [line_ex_tmp(2,line_idx), line_ey_tmp(2,line_idx)]
                    line_not_found = .true.
                    if (norm2(leftLinePoint-search_point).lt.1d-14) then
                        line_ex_red(:,i)  = line_ex_tmp(:,line_idx)
                        line_ey_red(:,i)  = line_ey_tmp(:,line_idx)
                        tplines_sorted(i) = tplines(line_idx)
                        visited_lines(i)  = line_idx
                        line_not_found    = .false.
                        exit
                    elseif (norm2(rightLinePoint-search_point).lt.1d-14) then
                        line_ex_red(:,i)  = line_ex_tmp([2,1],line_idx)
                        line_ey_red(:,i)  = line_ey_tmp([2,1],line_idx)
                        tplines_sorted(i) = tplines(line_idx)
                        visited_lines(i)  = line_idx
                        line_not_found    = .false.      
                        exit         
                    endif
                endif
            enddo

            if (line_not_found) then
                exit
            endif

            nvisited2 = count(visited_lines.ne.0)
        enddo

        nvisited               = nvisited2
        k                      = k + 1
        separate_lines_idx2(k) = nvisited

        ! Set values in line_xx_remaining at visited_lines to large number s.t. not deteceted in finiding start_point for next
        ! discontinous set of line segments
        line_ex_remaining(:,visited_lines(1:nvisited)) = 1d15
        line_ey_remaining(:,visited_lines(1:nvisited)) = 1d15
    
    enddo

    ! Save separate_lines_idx    
    allocate(separate_lines_idx(k-1),stat=ierr)
    separate_lines_idx = separate_lines_idx2(1:k-1)

    ! Assign line_ex and line_ey to line_ex_red and line_ey_red    
    line_ex(1:nseg,:) = transpose(line_ex_red)
    line_ey(1:nseg,:) = transpose(line_ey_red)

    ! Assign tplines
    tplines = tplines_sorted

    ! Determine if lines seperated in different continous sets or not
    nseparate_lines = size(separate_lines_idx)
    separate_lines  = nseparate_lines.gt.0


    if (separate_lines) then
        ! Multiple sets of continous line segments

        ! Find sequences of connecting segments
        n_seperate_sets                                    = nseparate_lines + 1
        allocate(separate_lines_idx_extended(2,n_seperate_sets),stat=ierr)
        separate_lines_idx_extended(1,1)                   = 1
        separate_lines_idx_extended(1,2:n_seperate_sets)   = separate_lines_idx + 1
        separate_lines_idx_extended(2,1:n_seperate_sets-1) = separate_lines_idx
        separate_lines_idx_extended(2,n_seperate_sets)     = nseg

    else
        n_seperate_sets = 1
        allocate(separate_lines_idx_extended(2,1),stat=ierr)
        separate_lines_idx_extended(:,1) = [1, nseg]
    endif

    ! Deallocate
    deallocate(line_ex_tmp)
    deallocate(line_ey_tmp)
    deallocate(line_ex_red)
    deallocate(line_ey_red)
    deallocate(visited_lines)
    deallocate(separate_lines_idx2)
    deallocate(tplines_sorted)

    return
end subroutine connect_line_segments


subroutine remove_line_segments(line_ex, line_ey, tplines, critical_length, separate_lines_idx, line_seg, &
    separate_lines_idx_extended)
    ! Routine for removing line segments shorter than threshold
    implicit none

    ! Intent inout
    real(dp), intent(inout) :: line_ex(:,:), line_ey(:,:)
    integer, intent(inout)  :: line_seg, separate_lines_idx(:)
    integer, allocatable, intent(inout) :: separate_lines_idx_extended(:,:)

    ! Intent in    
    real(dp), intent(in)    :: critical_length
    logical, intent(in)     :: tplines(:)

    ! Subroutine variables
    real(dp), allocatable   :: line_ex_new(:,:), line_ey_new(:,:), line_lengths(:)    
    
    logical, allocatable    :: keep_lines(:), closed_sets(:)
    real(dp)                :: new_x, new_y, A(2), B(2)
    integer                 :: nkeep, i, ierr, next_line, nseg, nseparate_lines, n_seperate_sets, iseg_start, iseg_end, k
    logical                 :: separate_lines
    real(dp)                :: Lseg1, Lseglast, rmedge

    ! Number of line segments in level set function
    nseg = line_seg

    ! Length of line segments
    allocate(line_lengths(nseg),stat=ierr)
    line_lengths = sqrt((line_ex(:,1)-line_ex(:,2))**2 + (line_ey(:,1)-line_ey(:,2))**2)

    ! Lines to be kept (logical)
    allocate(keep_lines(nseg), stat=ierr)
    keep_lines = .false.
    keep_lines = line_lengths.gt.critical_length .or. tplines(1:nseg)

    ! Determine if lines seperated in different continous sets or not
    nseparate_lines = size(separate_lines_idx)
    separate_lines  = nseparate_lines.gt.0


    if (separate_lines) then
        ! Multiple sets of continous line segments

        ! Find sequences of connecting segments
        n_seperate_sets                                    = nseparate_lines + 1
        allocate(separate_lines_idx_extended(2,n_seperate_sets),stat=ierr)
        separate_lines_idx_extended(1,1)                   = 1
        separate_lines_idx_extended(1,2:n_seperate_sets)   = separate_lines_idx + 1
        separate_lines_idx_extended(2,1:n_seperate_sets-1) = separate_lines_idx
        separate_lines_idx_extended(2,n_seperate_sets)     = nseg

    else
        n_seperate_sets = 1
        allocate(separate_lines_idx_extended(2,1),stat=ierr)
        separate_lines_idx_extended(:,1) = [1, nseg]
    endif

    ! Allocate closed curve sets (logical)
    allocate(closed_sets(n_seperate_sets),stat=ierr)

    ! Loop through each continous set and check if closed curve
    do k = 1,n_seperate_sets
        iseg_start = separate_lines_idx_extended(1,k)
        iseg_end   = separate_lines_idx_extended(2,k)

        ! Check if set is a closed curve
        closed_sets(k) = (abs(line_ex(iseg_end,2) - line_ex(iseg_start,1)).lt.1d-14) .and. &
                         (abs(line_ey(iseg_end,2) - line_ey(iseg_start,1)).lt.1d-14)

        ! If segment set is not closed, make sure first and last line segment kept
        if (.not. closed_sets(k)) then
            keep_lines(iseg_start) = 1
            keep_lines(iseg_end)   = 1
        endif

        if (k.eq.1) then
            separate_lines_idx(k) = count(keep_lines(iseg_start:iseg_end))
        elseif (k .lt. n_seperate_sets) then
            separate_lines_idx(k) = separate_lines_idx(k-1) + count(keep_lines(iseg_start:iseg_end))
        endif
    enddo        

    ! Obtain a new, reduced, set of lines where the short lines have been removed. I.e. only keep lines listed in keep_lines
    nkeep = count(keep_lines)    
    allocate(line_ex_new(2,nkeep), stat=ierr)
    allocate(line_ey_new(2,nkeep), stat=ierr)
    line_ex_new(1,:) = pack(line_ex(:,1),keep_lines)
    line_ex_new(2,:) = pack(line_ex(:,2),keep_lines)
    line_ey_new(1,:) = pack(line_ey(:,1),keep_lines)
    line_ey_new(2,:) = pack(line_ey(:,2),keep_lines)
    print *, 'Removing', nseg-nkeep, 'lines'



    ! --- Loop through reduced line sets and make them connect ---

    if (separate_lines) then
        ! Multiple sets of continous line segments

        ! Find sequences of connecting segments. Observe that separate_lines_idx has changed after reducing lines
        separate_lines_idx_extended(1,1)                   = 1
        separate_lines_idx_extended(1,2:n_seperate_sets)   = separate_lines_idx + 1
        separate_lines_idx_extended(2,1:n_seperate_sets-1) = separate_lines_idx
        separate_lines_idx_extended(2,n_seperate_sets)     = nkeep
    else
        separate_lines_idx_extended(:,1) = [1, nkeep]

    endif

    do k=1,n_seperate_sets

        ! Start and fininsh line segment
        iseg_start = separate_lines_idx_extended(1,k)
        iseg_end   = separate_lines_idx_extended(2,k)

        ! If set not supposed to be closed, dont check first and last
        if (.not. closed_sets(k)) then
            iseg_end = iseg_end - 1
        endif

        ! Loop through segments and make sure they connect  
        do i=iseg_start,iseg_end
            
            next_line = i + 1
            if (i.eq.iseg_end .and. closed_sets(k)) then
                next_line = iseg_start
            endif
        
            ! If endpoint not coincide with next line
            A = [line_ex_new(2,i), line_ey_new(2,i)]
            B = [line_ex_new(1,next_line), line_ey_new(1,next_line)]

            if (norm2(A-B).gt.1d-14) then 
                
                ! New x and y position
                new_x                    = (line_ex_new(2,i) + line_ex_new(1,next_line))/2d0
                new_y                    = (line_ey_new(2,i) + line_ey_new(1,next_line))/2d0

                ! Update lines
                line_ex_new(2,i)         = new_x
                line_ex_new(1,next_line) = new_x
                line_ey_new(2,i)         = new_y
                line_ey_new(1,next_line) = new_y
            endif
        enddo
    enddo    

    ! Adjust line_ex and line_ey and assign line_ex_new and line_ey_new to them
    line_ex               = 0d0
    line_ey               = 0d0
    line_seg              = nkeep
    line_ex(1:line_seg,:) = transpose(line_ex_new)
    line_ey(1:line_seg,:) = transpose(line_ey_new)    

    

    ! Deallocate
    deallocate(line_ex_new)
    deallocate(line_ey_new)
    deallocate(line_lengths)
    deallocate(keep_lines)


    ! Remove end segments if too short    
    do k = 1,n_seperate_sets
        if (.not. closed_sets(k)) then            

            rmedge = 0

            ! First line segment
            Lseg1 = sqrt((line_ex(1,2) - line_ex(1,1))**2 + (line_ey(1,2) - line_ey(1,1))**2)
            if (Lseg1.lt.critical_length) then

                ! Remove first line segment
                allocate(line_ex_new(line_seg-1,2), stat=ierr)
                allocate(line_ey_new(line_seg-1,2), stat=ierr)
                line_ex_new = 0d0
                line_ey_new = 0d0
                line_ex_new = line_ex(2:line_seg,:)
                line_ey_new = line_ey(2:line_seg,:)
                line_ex_new(1,1) = line_ex(1,1)
                line_ey_new(1,1) = line_ey(1,1)
                line_seg = line_seg - 1

                ! Adjust line_ex and line_ey and assign line_ex_new and line_ey_new to them
                line_ex               = 0d0
                line_ey               = 0d0
                line_ex(1:line_seg,:) = line_ex_new
                line_ey(1:line_seg,:) = line_ey_new
                separate_lines_idx_extended(2,k) = line_seg

                deallocate(line_ex_new)
                deallocate(line_ey_new)
            endif

            
            ! Last line segment
            Lseglast = sqrt((line_ex(line_seg,2) - line_ex(line_seg,1))**2 + (line_ey(line_seg,2) - line_ey(line_seg,1))**2)
            if (Lseglast.lt.critical_length) then

                ! Remove last line segment
                allocate(line_ex_new(line_seg-1,2), stat=ierr)
                allocate(line_ey_new(line_seg-1,2), stat=ierr)
                line_ex_new = 0d0
                line_ey_new = 0d0
                line_ex_new = line_ex(1:line_seg-1,:)
                line_ey_new = line_ey(1:line_seg-1,:)
                line_ex_new(line_seg-1,2) = line_ex(line_seg,2)
                line_ey_new(line_seg-1,2) = line_ey(line_seg,2)
                line_seg = line_seg - 1

                ! Adjust line_ex and line_ey and assign line_ex_new and line_ey_new to them
                line_ex               = 0d0
                line_ey               = 0d0
                line_ex(1:line_seg,:) = line_ex_new
                line_ey(1:line_seg,:) = line_ey_new
                separate_lines_idx_extended(2,k) = line_seg

                deallocate(line_ex_new)
                deallocate(line_ey_new)
            endif
            ! if (Lseglast.lt.critical_length) then
            !     ! Remove line 
            !     rm2 = line_seg
            !     line_ex(2,1) = line_ex(1,1)
            ! endif

        endif
    enddo    



    return
end subroutine remove_line_segments


subroutine add_line_segments(line_ex, line_ey, critical_length, line_seg, sep_lines, nsep_lines)
    ! Routine for adding node in middle of line segment if line longer than threshold

    ! Intent inout
    real(dp), intent(inout) :: line_ex(:,:), line_ey(:,:)
    integer, intent(inout)  :: line_seg, sep_lines(:,:)

    ! Intent in
    real(dp), intent(in)    :: critical_length
    integer, intent(in)     :: nsep_lines

    ! Subroutine variables
    real(dp), allocatable   :: line_ex_new(:,:), line_ey_new(:,:)
    real(dp), allocatable   :: line_ex_new_copy(:,:), line_ey_new_copy(:,:), line_lengths(:)
    real(dp)                :: new_x, new_y, A(2), B(2)
    integer                 :: nseg, nseg_extended, n_add_lines, k, ierr, line_idx, i, isegA, isegB
    real(dp), allocatable   :: sep_points(:,:), jk(:) 

    ! Number of line segments in level set function
    nseg = line_seg

    ! Number of extended line segments
    nseg_extended = 10*nseg

    ! Allocate sep_points
    allocate(sep_points(nsep_lines, 4),stat=ierr)

    ! Save seperated points
    do i=1,nsep_lines
        isegA = sep_lines(i,1)
        isegB = sep_lines(i,2)
        sep_points(i,1) = line_ex(isegA,1)
        sep_points(i,2) = line_ey(isegA,1)
        sep_points(i,3) = line_ex(isegB,2)
        sep_points(i,4) = line_ey(isegB,2)
    enddo

    

    ! Allocate new arrays for storing old lines as well as new added lines
    allocate(line_ex_new(2,nseg_extended),stat=ierr)
    allocate(line_ey_new(2,nseg_extended),stat=ierr)
    line_ex_new(:,1:nseg) = transpose(line_ex(1:nseg,:))
    line_ey_new(:,1:nseg) = transpose(line_ey(1:nseg,:))

    ! Allocate copy of new lines
    allocate(line_ex_new_copy(2,nseg_extended),stat=ierr)
    allocate(line_ey_new_copy(2,nseg_extended),stat=ierr)

    ! Length of line segments
    allocate(line_lengths(nseg_extended),stat=ierr)
    line_lengths(1:nseg) = sqrt((line_ex_new(1,1:nseg)-line_ex_new(2,1:nseg))**2 + (line_ey_new(1,1:nseg)-line_ey_new(2,1:nseg))**2)

    ! Number of lines to add
    n_add_lines = count(line_lengths(1:nseg).gt.critical_length)

    ! Add lines as until no line segemnt longer than threshold
    do while (n_add_lines.gt.0)

      ! Print info
      print *, 'Number of lines added: ', n_add_lines

      ! Copy lines
      line_ex_new_copy = line_ex_new
      line_ey_new_copy = line_ey_new

      ! Loop through lines and add point in middle
      k = 0
      do line_idx=1,nseg
        if (line_lengths(line_idx).gt.critical_length) then
          k = k + 2
          new_x              = (line_ex_new_copy(1,line_idx) + line_ex_new_copy(2,line_idx))/2d0
          new_y              = (line_ey_new_copy(1,line_idx) + line_ey_new_copy(2,line_idx))/2d0
          line_ex_new(:,k-1) = [line_ex_new_copy(1,line_idx), new_x]
          line_ey_new(:,k-1) = [line_ey_new_copy(1,line_idx), new_y]
          line_ex_new(:,k)   = [new_x,line_ex_new_copy(2,line_idx)]
          line_ey_new(:,k)   = [new_y,line_ey_new_copy(2,line_idx)]
        else 
          k = k + 1
          line_ex_new(:,k)   = line_ex_new_copy(:,line_idx)
          line_ey_new(:,k)   = line_ey_new_copy(:,line_idx)
        endif
      enddo

      ! New number of line segments
      nseg = nseg + n_add_lines

      ! Compute n_add_lines for new line segments
      line_lengths(1:nseg) = sqrt((line_ex_new(1,1:nseg)-line_ex_new(2,1:nseg))**2 + (line_ey_new(1,1:nseg)-line_ey_new(2,1:nseg))**2)      
      n_add_lines = count(line_lengths(1:nseg).gt.critical_length)

    enddo

    if (nseg.gt.nseg_extended) then
      print *, 'Error in add_lines - to few nseg_extended'
    endif

    ! Adjust line_ex and line_ey and assign line_ex_red and line_ey_red to them
    line_ex               = 0d0
    line_ey               = 0d0
    line_seg              = nseg
    line_ex(1:line_seg,:) = transpose(line_ex_new(:,1:nseg))
    line_ey(1:line_seg,:) = transpose(line_ey_new(:,1:nseg))


    ! Find saved seperated points
    do i=1,nsep_lines
        A = [sep_points(i,1), sep_points(i,2)]
        B = [sep_points(i,3), sep_points(i,4)]

        ! Find the index of the first true value of line containing A
        jk = pack([(k, k=1, nseg)],count(abs(line_ex(1:nseg,:) - A(1) + line_ey(1:nseg,:) - A(2)) .lt. 1d-14,2).eq.1)
        isegA = jk(1)

        ! Find the index of the first true value of line containing B
        jk = pack([(k, k=1, nseg)],count(abs(line_ex(1:nseg,:) - B(1) + line_ey(1:nseg,:) - B(2)) .lt. 1d-14,2).eq.1)        
        isegB = jk(size(jk))

        ! Edit sep_lines
         sep_lines(i,1) = isegA
         sep_lines(i,2) = isegB


    enddo
    
    ! Deallocate
    deallocate(line_ex_new)
    deallocate(line_ey_new)
    deallocate(line_ex_new_copy)
    deallocate(line_ey_new_copy)
    deallocate(line_lengths)

    return
end subroutine add_line_segments


end module ls_sorting_routines
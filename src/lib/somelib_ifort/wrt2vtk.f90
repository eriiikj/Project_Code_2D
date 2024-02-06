module wrt2vtk
 !============================================================================
!
!   Module wrt2vt contains subroutines for writing Fortran data to vtk files.
!   Only unstructured grid is implemented
!
!   Included routines:
!
!      VTKopen                    Open a VTK file and assign a pointer and a filename
!      VTKclose                   Close a VTK file
!      VTKPlotMesh                Creates an unstructured mesh (2D or 3D)
!      VTKPlotNodeVal             Write nodal values to an unstructured mesh
!      VTKPlotCellVal             Write cell values to an unstructured mesh
!     
!
!   Last rev.:
!      M. Ristinmaa, 2020-02-25
!      M. Ristinmaa, 2020-03-11
!          Changed to binary format, allowing several fields to be stored
!          Included plotting of gauss point data average over element
! ============================================================================

implicit none

character(len=1)                   :: lf=char(10)
character(len=80)                  :: cell_name
integer                            :: num_elem

private lf, cell_name

integer    :: ivtk=10, initate_point_data, initate_cell_data, lock_point_data, lock_cell_data
private ivtk, initate_point_data, initate_cell_data, lock_point_data, lock_cell_data

interface VTKOpen
   module procedure VTKopen1
   module procedure VTKopen2
end interface
private VTKopen1, VTKopen2

interface VTKPlotCellVal
   module procedure VTKPlotCellVal1
   module procedure VTKPlotCellVal2
   module procedure VTKPlotCellVal3
!   module procedure VTKPlotCellVal4
end interface
private VTKPlotCellVal1, VTKPlotCellVal2, VTKPlotCellVal3

interface VTKPlotCellVectorEL
   module procedure VTKPlotCellVectorEL1
   module procedure VTKPlotCellVectorEL2
end interface
private VTKPlotCellVectorEL1, VTKPlotCellVectorEL2

interface VTKPlotCellTensorEL
   module procedure VTKPlotCellTensorEL1
   module procedure VTKPlotCellTensorEL2
end interface
private VTKPlotCellTensorEL1, VTKPlotCellTensorEL2


contains

! ============================================================================

subroutine VTKopen1(filename)
   ! ===========================================================
   ! This routine opens VTK file.
   !
   !    Input:
   !          filename          Name of the new VTK file
   !
   ! ===========================================================

   implicit none

   character(len=*)         :: filename
   character                :: buffer*80

      open(unit=ivtk,file=filename,form='binary',convert='BIG_ENDIAN')
      ! File header
      ! ===========

      buffer = '# vtk DataFile Version 3.0'//lf                                             ;
      write(ivtk) trim(buffer)
      buffer = 'Model data'//lf                                                             ;
      write(ivtk) trim(buffer)
      buffer = 'BINARY'//lf                                                                 ;
      write(ivtk) trim(buffer)
      buffer = 'DATASET UNSTRUCTURED_GRID'//lf                                              ;
      write(ivtk) trim(buffer)
      
      initate_point_data=0
      initate_cell_data=0
      lock_point_data=0 
      lock_cell_data=0

   return

end subroutine VTKopen1

subroutine VTKopen2(filename,num)
   ! ===========================================================
   ! This routine opens VTK file.
   !
   !    Input:
   !          filename          Name of the new VTK file
   !          num               extension to VTK file filename.num
   !                            used for time series
   !
   ! ===========================================================

   implicit none

   character(len=*)         :: filename
   character(len=80)        :: filee, buffer
   character(len=30)        :: str
   integer                  :: num

!write(ci,'(I10)') step
!write(filename,'(A5,A10)') 'step_',adjustl(ci)
!filename = trim(filename) // '.vtk'
!open(unit=ivtk,file=filename,form='binary',convert='BIG_ENDIAN')

      write(str,*)num
      filee=trim(filename)//'.'//adjustl(str)
 
      open(unit=ivtk,file=filee,form='binary',convert='BIG_ENDIAN')
      ! File header
      ! ===========

      buffer = '# vtk DataFile Version 3.0'//lf                                             ;
      write(ivtk) trim(buffer)
      buffer = 'Model data'//lf                                                             ;
      write(ivtk) trim(buffer)
      buffer = 'BINARY'//lf                                                                 ;
      write(ivtk) trim(buffer)
      buffer = 'DATASET UNSTRUCTURED_GRID'//lf                                              ;
      write(ivtk) trim(buffer)
      
      initate_point_data=0
      initate_cell_data=0
      lock_point_data=0 
      lock_cell_data=0

   return

end subroutine VTKopen2

subroutine VTKclose()
   ! ===========================================================
   ! This routine closes VTK file.
   !
   ! ===========================================================

   implicit none

   close(unit=ivtk)
   cell_name=''

   return

end subroutine VTKclose

subroutine VTKplotMesh(eltype,enod,coord)

   implicit none

   character                           :: eltype*(*)
   integer, allocatable                :: enod(:,:)
   double precision, allocatable       :: coord(:,:)
   integer                             :: ie, im, nrelem, nrnode, nodeel, vtktype
!   character(len=10)                   :: str
   character                           :: buffer*80, str1*8, str2*8

   nrelem=size(enod,2)
   nodeel=size(enod,1)
   nrnode=size(coord,2)

! store internally used in cell plot
   num_elem=nrelem

! write coordinates for nodes
      write(str1(1:8),'(i8)') nrnode
      buffer = 'POINTS '//str1//'  double'//lf                                              ;
      write(ivtk) trim(buffer)

      if (size(coord,1).eq.3) then
        do im=1,nrnode
          write(ivtk)coord(1,im),coord(2,im),coord(3,im)
        enddo
      else
        do im=1,nrnode
          write(ivtk)coord(1,im),coord(2,im),0.d0
        enddo
      endif
   
! write cells or elements
      write(str1(1:8),'(i8)') nrelem
      write(str2(1:8),'(i8)') (nodeel+1)*nrelem
      buffer = 'CELLS '//str1//' '//str2//lf                                                ;
      write(ivtk) trim(buffer)
      do ie=1,nrelem
         write(ivtk)nodeel,(enod(im,ie)-1,im=1,nodeel)
     enddo

! write vtk cell type
     if (eltype.eq.'tr3') vtktype=5
     if (eltype.eq.'qu4') vtktype=9
     if (eltype.eq.'brick8') vtktype=12

      write(str1(1:8),'(i8)') nrelem
      buffer = lf//lf//'CELL_TYPES '//str1//lf                                              ;
      write(ivtk) trim(buffer)
      do ie=1,nrelem
        write(ivtk)vtktype
      enddo

      buffer = 'FIELD Elementdata 1'//lf                                                    ;
      write(ivtk) trim(buffer)
      write(str1(1:8),'(i8)') nrelem
      buffer = 'Elementnummer  1 '//str1//' int'//lf                                        ;
      write(ivtk) trim(buffer)
      do ie=1,nrelem
        write(ivtk)ie
      enddo

      
end subroutine VTKPlotMesh

subroutine VTKPlotNodeVal(val,name)
   double precision                    :: val(:)
   character                           :: name*(*)
   integer                             :: nrnode, im
   character                           :: buffer*80, str1*8

     if (initate_cell_data.eq.1) lock_cell_data=1

     if (lock_point_data.eq.1) then
        write(*,*)'Please write all node data after each other'
        write(*,*)'Not possible to switch between node and cell data'
        stop
     end if

     nrnode=size(val)

! Only written one time in file
     if (initate_point_data.eq.0) then
       write(str1(1:8),'(i8)') nrnode
       buffer = lf//lf//'POINT_DATA '//str1//lf                                              ;
       write(ivtk) trim(buffer)
       initate_point_data=1
     end if
     
! Data  part
      buffer = 'SCALARS '//name//'  double'//lf                                             ;
      write(ivtk) trim(buffer)
      buffer = 'LOOKUP_TABLE default'//lf                                                   ;
      write(ivtk) trim(buffer)
      write(ivtk) (val(im),im=1,nrnode)



end subroutine VTKPlotNodeVal

subroutine VTKPlotCellVal1(val,name)
   double precision                    :: val(:)
   character                           :: name*(*)
   integer                             :: nrelem, im
   character                           :: buffer*80, str1*8

     if (initate_point_data.eq.1) lock_point_data=1

      if (lock_cell_data.eq.1) then
        write(*,*)'Please write all cell data after each other'
        write(*,*)'Not possible to switch between node and cell data'
        stop
     end if

    nrelem=size(val)

! Only written one time in file
     if (initate_cell_data.eq.0) then
       write(str1(1:8),'(i8)') nrelem
       buffer = lf//lf//'CELL_DATA '//str1//lf                                              ;
       write(ivtk) trim(buffer)
       initate_cell_data=1
     end if
     
     
! Data  part
      buffer = 'SCALARS '//name//'  double'//lf                                             ;
      write(ivtk) trim(buffer)
      buffer = 'LOOKUP_TABLE default'//lf                                                   ;
      write(ivtk) trim(buffer)
      write(ivtk) (val(im),im=1,nrelem)


end subroutine VTKPlotCellVal1

subroutine VTKPlotCellVal2(val,name)
! stores the average gauss value in the element
! In paraview use filter/alphabetical 
! select !cell data to point data"
! to get a smooth plot
!
   double precision                    :: val(:,:)
   character                           :: name*(*)
   integer                             :: nrelem, nrgauss, im
   character                           :: buffer*80, str1*8

     if (initate_point_data.eq.1) lock_point_data=1

      if (lock_cell_data.eq.1) then
        write(*,*)'Please write all cell data after each other'
        write(*,*)'Not possible to switch between node and cell data'
        stop
     end if

     nrelem=size(val,2)
     nrgauss=size(val,1)

! Only written one time in file
     if (initate_cell_data.eq.0) then
       write(str1(1:8),'(i8)') nrelem
       buffer = lf//lf//'CELL_DATA '//str1//lf                                              ;
       write(ivtk) trim(buffer)
       initate_cell_data=1
     end if
     
! Data  part
      buffer = 'SCALARS '//name//'  double'//lf                                             ;
      write(ivtk) trim(buffer)
      buffer = 'LOOKUP_TABLE default'//lf                                                   ;
      write(ivtk) trim(buffer)
      write(ivtk) (sum(val(:,im))/dble(nrgauss),im=1,nrelem)


end subroutine VTKPlotCellVal2


subroutine VTKPlotCellVal3(val,name)
   integer                             :: val(:)
   character                           :: name*(*)
   integer                             :: nrelem, im
   character                           :: buffer*80, str1*8

     if (initate_point_data.eq.1) lock_point_data=1

      if (lock_cell_data.eq.1) then
        write(*,*)'Please write all cell data after each other'
        write(*,*)'Not possible to switch between node and cell data'
        stop
     end if

     nrelem=size(val)

! Only written one time in file
     if (initate_cell_data.eq.0) then
       write(str1(1:8),'(i8)') nrelem
       buffer = lf//lf//'CELL_DATA '//str1//lf                                              ;
       write(ivtk) trim(buffer)
       initate_cell_data=1
     end if
     
! Data  part
      buffer = 'SCALARS '//name//'  int 1'//lf                                              ;
      write(ivtk) trim(buffer)
      buffer = 'LOOKUP_TABLE default'//lf                                                   ;
      write(ivtk) trim(buffer)
      write(ivtk) (val(im),im=1,nrelem)


end subroutine VTKPlotCellVal3


! Routines for single element calls

subroutine VTKPlotCellValEl(val,name)
!subroutine VTKPlotCellValEl(val,name,nrelem)
!subroutine VTKPlotCellVal4(val,name,new,nrelem)
! stores the average gauss value in one element
! note must be called for all elements in 
! correct order
! In paraview use filter/alphabetical 
! select !cell data to point data"
! to get a smooth plot
!
   double precision                    :: val(:)
   character                           :: name*(*)
!   integer, optional                   :: nrelem
   integer                             :: nrgauss
   character                           :: buffer*80, str1*8

     if (initate_point_data.eq.1) lock_point_data=1

      if (lock_cell_data.eq.1) then
        write(*,*)'Please write all cell data after each other'
        write(*,*)'Not possible to switch between node and cell data'
        stop
     end if

     nrgauss=size(val)

! Only written one time in file
     if (initate_cell_data.eq.0) then
!       if(present(num_elem))then ! as we must call plotmesh first we have the number
         write(str1(1:8),'(i8)') num_elem   !  this is undefined must be given problem 
         buffer = lf//lf//'CELL_DATA '//str1//lf                                              ;
         write(ivtk) trim(buffer)
         initate_cell_data=1        
!       else
!         write(*,*)'number of elements must be give as this is first call to PlotCellVal'
!         stop
!       endif       
     end if
     
! Data  part initiate if new
!     if (new.eq.1) then
!       buffer = 'SCALARS '//name//'  double'//lf                                             ;
!       write(ivtk) trim(buffer)
!       buffer = 'LOOKUP_TABLE default'//lf                                                   ;
!       write(ivtk) trim(buffer)
!     end if

     if (0.eq.(cell_name.eq.name)) then
       buffer = 'SCALARS '//name//'  double'//lf                                             ;
       write(ivtk) trim(buffer)
       buffer = 'LOOKUP_TABLE default'//lf                                                   ;
       cell_name=trim(name)
     end if

! Data
     write(ivtk) sum(val)/dble(nrgauss)


end subroutine VTKPlotCellValEl


subroutine VTKPlotCellScalarEL(val,name)
!
   double precision                    :: val(:)
   character                           :: name*(*)
!   integer, optional                   :: nrelem
   character                           :: buffer*80, str1*8


     if (initate_point_data.eq.1) lock_point_data=1

      if (lock_cell_data.eq.1) then
        write(*,*)'Please write all cell data after each other'
        write(*,*)'Not possible to switch between node and cell data'
        stop
     end if

! Only written one time in file
     if (initate_cell_data.eq.0) then
!       if(present(num_elem))then
         write(str1(1:8),'(i8)') num_elem   
         buffer = lf//lf//'CELL_DATA '//str1//lf                                              ;
         write(ivtk) trim(buffer)
         initate_cell_data=1        
!       else
!         write(*,*)'number of elements must be give as this is first call to PlotCellVal'
!         stop
!       endif       
     end if
     
! Data  part initiate if new field
!     if (new.eq.1) then
     if (0.eq.(cell_name.eq.name)) then
       buffer = 'SCALARS '//name//'  double '// '1' //lf                                             ;
       write(ivtk) trim(buffer)
       buffer = 'LOOKUP_TABLE default'//lf
       write(ivtk) trim(buffer)
       cell_name=trim(name)
     end if

! Data
     write(ivtk) sum(val(:))/dble(size(val))

end subroutine VTKPlotCellScalarEL


subroutine VTKPlotCellVectorEL1(val,name)
!subroutine VTKPlotCellVectorEL(val,name,nrelem)
!subroutine VTKPlotCellVector(val,name,new,nrelem)
! stores the average gauss value in one element
! note must be called for all elements in 
! correct order
! In paraview use filter/alphabetical 
! select "cell data to point data"
! to get a smooth plot
!  val(components,gauss points) 
!
   double precision                    :: val(:,:)
   character                           :: name*(*)
!   integer, optional                   :: nrelem
   integer                             :: nrgauss, nrval, iv
   character                           :: buffer*80, str1*8

     if (initate_point_data.eq.1) lock_point_data=1

      if (lock_cell_data.eq.1) then
        write(*,*)'Please write all cell data after each other'
        write(*,*)'Not possible to switch between node and cell data'
        stop
     end if

     nrgauss=size(val,2)
     nrval=size(val,1)

! Only written one time in file
     if (initate_cell_data.eq.0) then
!       if(present(num_elem))then
         write(str1(1:8),'(i8)') num_elem   !  this is undefined must be given problem 
         buffer = lf//lf//'CELL_DATA '//str1//lf                                              ;
         write(ivtk) trim(buffer)
         initate_cell_data=1        
!       else
!         write(*,*)'number of elements must be give as this is first call to PlotCellVal'
!         stop
!       endif       
     end if
     
! Data  part initiate if new field
!     if (new.eq.1) then
     if (0.eq.(cell_name.eq.name)) then
       buffer = 'VECTORS '//name//'  double'//lf                                             ;
       write(ivtk) trim(buffer)
       cell_name=trim(name)
     end if

! Data
     if (nrval.eq.3) then
       write(ivtk) (sum(val(iv,:))/dble(nrgauss),iv=1,nrval)
     elseif (nrval.eq.2) then
       write(ivtk) (sum(val(iv,:))/dble(nrgauss),iv=1,nrval), 0d0
     else
       write(*,*)'VTKPlotCellVector: Wrong size first index first argument, 2 or 3'  
       stop
     end if

end subroutine VTKPlotCellVectorEL1

subroutine VTKPlotCellVectorEL2(val,name)
! Special for 1 gauss point
!
   double precision                    :: val(:), valm(size(val),1)
   character                           :: name*(*)

   valm(:,1)=val
   call VTKPlotCellVectorEL(valm,name)

end subroutine VTKPlotCellVectorEL2

subroutine VTKPlotCellTensorEL1(val,name)
! stores the average gauss value in one element
! note must be called for all elements in 
! correct order
! In paraview use filter/alphabetical 
! select "cell data to point data"
! to get a smooth plot
!  val(components,gauss points) 
!
   double precision                    :: val(:,:)
   character                           :: name*(*)
!   integer, optional                   :: nrelem
   integer                             :: nrgauss, nrval
   character                           :: buffer*80, str1*8

   double precision                    :: t11, t22, t33, t12, t13, t23

     if (initate_point_data.eq.1) lock_point_data=1

      if (lock_cell_data.eq.1) then
        write(*,*)'Please write all cell data after each other'
        write(*,*)'Not possible to switch between node and cell data'
        stop
     end if

     nrgauss=size(val,2)
     nrval=size(val,1)

! Only written one time in file
     if (initate_cell_data.eq.0) then
!       if(present(num_elem))then
         write(str1(1:8),'(i8)') num_elem   
         buffer = lf//lf//'CELL_DATA '//str1//lf                                              ;
         write(ivtk) trim(buffer)
         initate_cell_data=1        
!       else
!         write(*,*)'number of elements must be give as this is first call to PlotCellVal'
!         stop
!       endif       
     end if
     
! Data  part initiate if new field
!     if (new.eq.1) then
     if (0.eq.(cell_name.eq.name)) then
       buffer = 'TENSORS '//name//'  double'//lf                                             ;
       write(ivtk) trim(buffer)
       cell_name=trim(name)
     end if

! Data
     if (nrval.eq.6) then
       t11=sum(val(1,:))/dble(nrgauss)
       t22=sum(val(2,:))/dble(nrgauss)
       t33=sum(val(3,:))/dble(nrgauss)
       t12=sum(val(4,:))/dble(nrgauss)
       t13=sum(val(5,:))/dble(nrgauss)
       t23=sum(val(6,:))/dble(nrgauss)
       write(ivtk) t11, t12, t13
       write(ivtk) t12, t22, t23
       write(ivtk) t13, t13, t33
     elseif (nrval.eq.4) then
       t11=sum(val(1,:))/dble(nrgauss)
       t22=sum(val(2,:))/dble(nrgauss)
       t33=sum(val(3,:))/dble(nrgauss)
       t12=sum(val(4,:))/dble(nrgauss)
       t13=0
       t23=0
       write(ivtk) t11, t12, t13
       write(ivtk) t12, t22, t23
       write(ivtk) t13, t13, t33
     else
       write(*,*)'VTKPlotTensorVector: Wrong size first index first argument, 6 or 4'  
       stop
     end if

end subroutine VTKPlotCellTensorEL1

subroutine VTKPlotCellTensorEL2(val,name)
! Special for 1 gauss point
!
   double precision                    :: val(:), valm(size(val),1)
   character                           :: name*(*)

   valm(:,1)=val
   call VTKPlotCellTensorEL(valm,name)

end subroutine VTKPlotCellTensorEL2



end module wrt2vtk
